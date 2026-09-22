"""Guided temporal interpolation of the OPERA TROPO wet delay: hourly ERA5 and NEXRAD.

TROPO products exist every 6 hours, and a SAR acquisition usually falls between
two of them.  Linear interpolation assumes the wet delay changed along a straight
line over those 6 hours.  The guide corrects that straight line with two sources
that do see the weather in between:

* **ERA5**, hourly total column water vapour (TCWV).  Its departure from its own
  straight line between the two product times, as a fraction, scales the linearly
  interpolated wet delay::

      rho = E(t) / [(1 - w) E(t0) + w E(t1)] - 1

  Only ERA5's change in time is used, so its own biases cancel.

* **NEXRAD**, the national reflectivity composite (5 minutes, 0.005 deg), where
  available (contiguous US).  Radar sees rain, not vapour, so it guides in two
  indirect ways:

  - the change of the fraction of each cell covered by echoes >= 20 dBZ,
    smoothed over 40 km, relative to its own straight line; and
  - the motion of the echoes, tracked hour by hour through the interval and used
    to move the TROPO wet delay field (as in `interpolate_motion`, but with the
    displacement taken from radar instead of from the two TROPO fields).

The three terms are added to the linear wet delay with fixed weights::

    wet = wet_lin + a * rho * wet_lin + b * dcover * wet_lin / wet_lin(0 m)
                  + c * (wet_moved - wet_lin)

The weights were fitted against 5-minute GNSS zenith delays on ten dates over the
contiguous US (60,382 station-times) and then tested, unchanged, on four
Sentinel-1 frames: the guide lowered the error against GNSS by 4.6 % (Los Angeles),
7.4 % and 4.2 % (Puget Sound, two tracks) and 1.4 % (Houston, 23 minutes from a
product), and improved 58 of 59 Los Angeles DISP-S1 products (6.3 % lower
plane-removed residual).  Without radar, ERA5 alone gives roughly 60 % of the gain.

The guide corrects timing only.  It does nothing for moisture that the model puts
in the wrong place, which is the larger error in many regions.  The hydrostatic
delay is interpolated linearly.

"""

from __future__ import annotations

import logging
import urllib.request
from datetime import datetime, timedelta
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import rasterio
import xarray as xr
from rasterio.windows import Window, from_bounds
from scipy import ndimage

from ._motion import estimate_flow, warp

if TYPE_CHECKING:
    from collections.abc import Iterable, Sequence

__all__ = [
    "DEFAULT_GUIDE_WEIGHTS",
    "download_era5_tcwv",
    "download_guide_inputs",
    "download_nexrad_n0q",
    "guide_radar_times",
    "interp_in_time_guided",
    "load_era5_tcwv",
    "read_nexrad_echo_fraction",
]

logger = logging.getLogger(__name__)

# Weights fitted out of sample (see the module docstring).  "cover" is in metres of
# zenith wet delay per unit change of the echo-covered fraction at sea level.
DEFAULT_GUIDE_WEIGHTS: dict[str, float] = {
    "era5": 0.411,
    "cover": 0.0156,
    "motion": 0.157,
}
# Weight of the ERA5 term when no radar is available.
DEFAULT_ERA5_ONLY_WEIGHT = 0.522

ECHO_DBZ = 20.0
COVER_SMOOTH_KM = 40.0
TRACK_SMOOTH_KM = 10.0
# Radar echoes are tracked on pixels of about this size (degrees).
TRACK_PIXEL_DEG = 0.14
KM_PER_DEG = 111.2

N0Q_URL = (
    "https://mesonet.agron.iastate.edu/archive/data/{t:%Y/%m/%d}/GIS/uscomp/"
    "n0q_{t:%Y%m%d%H%M}.{ext}"
)
# n0q palette index -> dBZ; index 0 means no echo or no data.
N0Q_DBZ_OFFSET, N0Q_DBZ_STEP = -32.5, 0.5


# --------------------------------------------------------------------------- times
def guide_radar_times(
    t0: datetime, t1: datetime, t: datetime, step: timedelta = timedelta(hours=1)
) -> list[datetime]:
    """Radar times needed to guide the interval [t0, t1] at time `t`.

    Every `step` from `t0` to `t1` (to track the echoes), plus `t` rounded to the
    5-minute cadence of the composite.
    """
    out = []
    x = pd.Timestamp(t0).to_pydatetime()
    while x <= t1:
        out.append(x)
        x = x + step
    ts = pd.Timestamp(t).round("5min").to_pydatetime()
    if ts not in out:
        out.append(ts)
    return sorted(out)


def _bracketing_epochs(t: datetime, hours: int = 6) -> tuple[datetime, datetime]:
    t = pd.Timestamp(t).tz_localize(None)
    t0 = t.floor(f"{hours}h")
    if t0 == t:
        t0 = t0 - pd.Timedelta(hours=hours)
    return t0.to_pydatetime(), (t0 + pd.Timedelta(hours=hours)).to_pydatetime()


def nexrad_path(directory: Path | str, t: datetime) -> Path:
    """Local path of the n0q composite for time `t` in `directory`."""
    return Path(directory) / f"n0q_{pd.Timestamp(t):%Y%m%d%H%M}.png"


# --------------------------------------------------------------------------- NEXRAD
def download_nexrad_n0q(
    datetimes: Iterable[datetime],
    output_dir: Path | str,
    product_interval_hours: int = 6,
) -> list[Path]:
    """Download the NEXRAD n0q composites needed to guide each datetime.

    Source: Iowa Environmental Mesonet national composite (PNG + world file,
    0.005 deg, every 5 minutes since 2010), public, no credentials.  For each
    datetime the hourly composites from the preceding to the following TROPO
    product time are fetched, plus the one nearest the datetime.

    Returns the paths of the composites now on disk.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    needed: set[datetime] = set()
    for dt in datetimes:
        t0, t1 = _bracketing_epochs(dt, product_interval_hours)
        needed.update(guide_radar_times(t0, t1, dt))
    have = []
    for t in sorted(needed):
        png = nexrad_path(output_dir, t)
        for ext, dst in (("wld", png.with_suffix(".wld")), ("png", png)):
            if dst.exists():
                continue
            try:
                urllib.request.urlretrieve(N0Q_URL.format(t=t, ext=ext), dst)
            except OSError as e:
                logger.warning(f"NEXRAD {t:%Y-%m-%d %H:%M} not available: {e}")
                dst.unlink(missing_ok=True)
                break
        if png.exists() and png.with_suffix(".wld").exists():
            have.append(png)
    return have


def _echo_fraction_on_grid(
    index: np.ndarray,
    transform,
    latitude: np.ndarray,
    longitude: np.ndarray,
    dbz_threshold: float = ECHO_DBZ,
) -> np.ndarray:
    """Fraction of mosaic pixels >= `dbz_threshold` in each cell of a lat/lon grid.

    `index` holds n0q palette indices; `transform` maps (col, row) to (lon, lat)
    of pixel corners.  Cells the mosaic does not cover are NaN.
    """
    rows, cols = index.shape
    lon_px = transform.c + (np.arange(cols) + 0.5) * transform.a
    lat_px = transform.f + (np.arange(rows) + 0.5) * transform.e
    lat_desc = latitude[0] > latitude[-1]
    lat_sorted = latitude if not lat_desc else latitude[::-1]
    dlat = float(np.abs(np.diff(latitude)).mean())
    dlon = float(np.abs(np.diff(longitude)).mean())
    iy = np.floor((lat_px - (lat_sorted[0] - dlat / 2)) / dlat).astype(np.int64)
    ix = np.floor((lon_px - (longitude[0] - dlon / 2)) / dlon).astype(np.int64)
    ry = np.flatnonzero((iy >= 0) & (iy < latitude.size))
    rx = np.flatnonzero((ix >= 0) & (ix < longitude.size))
    out = np.full((latitude.size, longitude.size), np.nan)
    if ry.size == 0 or rx.size == 0:
        return out
    cell = (iy[ry][:, None] * longitude.size + ix[rx][None, :]).ravel()
    idx = index[np.ix_(ry, rx)].ravel().astype(np.int64)
    dbz = N0Q_DBZ_OFFSET + N0Q_DBZ_STEP * idx
    echo = ((idx > 0) & (dbz >= dbz_threshold)).astype(float)
    n = np.bincount(cell, minlength=out.size)
    s = np.bincount(cell, weights=echo, minlength=out.size)
    with np.errstate(invalid="ignore", divide="ignore"):
        frac = np.where(n > 0, s / n, np.nan).reshape(out.shape)
    return frac[::-1] if lat_desc else frac


def read_nexrad_echo_fraction(
    path: Path | str,
    latitude: np.ndarray,
    longitude: np.ndarray,
    dbz_threshold: float = ECHO_DBZ,
) -> np.ndarray:
    """Read one n0q composite and aggregate it to the fraction of echo per grid cell.

    Only the window covering the grid is read.  Returns an array shaped
    (latitude, longitude) in the grid's own latitude order; NaN outside the mosaic.
    """
    pad = 0.1
    with rasterio.open(path) as src:
        win = (
            from_bounds(
                float(longitude.min()) - pad,
                float(latitude.min()) - pad,
                float(longitude.max()) + pad,
                float(latitude.max()) + pad,
                transform=src.transform,
            )
            .round_offsets()
            .round_lengths()
        )
        win = win.intersection(Window(0, 0, src.width, src.height))
        index = src.read(1, window=win)
        transform = src.window_transform(win)
    return _echo_fraction_on_grid(index, transform, latitude, longitude, dbz_threshold)


# --------------------------------------------------------------------------- ERA5
def download_era5_tcwv(
    bounds: tuple[float, float, float, float],
    datetimes: Iterable[datetime],
    output_dir: Path | str,
    margin_deg: float = 5.0,
) -> list[Path]:
    """Download hourly ERA5 total column water vapour for the days of `datetimes`.

    Needs the ``cdsapi`` package and Copernicus Climate Data Store credentials
    (``~/.cdsapirc``).  One file per month, ``era5_tcwv_YYYYMM.nc``, with all 24
    hours of each needed day over `bounds` (west, south, east, north) plus
    `margin_deg`.
    """
    try:
        import cdsapi  # noqa: PLC0415  (optional dependency)
    except ImportError as e:
        msg = "download_era5_tcwv needs the 'cdsapi' package: pip install cdsapi"
        raise ImportError(msg) from e
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    days = sorted({d.date() for dt in datetimes for d in _bracketing_epochs(dt)})
    w, s, e, n = bounds
    files = []
    for y, m in sorted({(d.year, d.month) for d in days}):
        out = output_dir / f"era5_tcwv_{y}{m:02d}.nc"
        files.append(out)
        if out.exists():
            continue
        request = {
            "product_type": ["reanalysis"],
            "variable": ["total_column_water_vapour"],
            "year": [str(y)],
            "month": [f"{m:02d}"],
            "day": [f"{d.day:02d}" for d in days if (d.year, d.month) == (y, m)],
            "time": [f"{h:02d}:00" for h in range(24)],
            "area": [n + margin_deg, w - margin_deg, s - margin_deg, e + margin_deg],
            "data_format": "netcdf",
            "download_format": "unarchived",
        }
        cdsapi.Client(quiet=True).retrieve(
            "reanalysis-era5-single-levels", request, str(out)
        )
    return files


def download_guide_inputs(
    datetimes: list[datetime],
    aoi_bounds: tuple[float, float, float, float],
    output_dir: Path = Path("tropo_guide"),
    era5: bool = True,
    nexrad: bool = True,
    margin_deg: float = 5.0,
) -> None:
    """Download what `crop_tropo(time_interpolation="guided")` needs.

    Hourly ERA5 total column water vapour goes to ``output_dir/era5`` (needs
    ``cdsapi`` and Climate Data Store credentials); NEXRAD n0q composites to
    ``output_dir/nexrad`` (public; contiguous US only).  Pass these two
    directories as `era5_file` and `nexrad_dir`.

    Parameters
    ----------
    datetimes : list[datetime]
        Acquisition times to be guided.
    aoi_bounds : tuple[float, float, float, float]
        AOI as (west, south, east, north) in degrees.
    output_dir : Path
        Parent directory for the two inputs.
    era5 : bool
        Download ERA5.
    nexrad : bool
        Download NEXRAD.
    margin_deg : float
        Margin in degrees added around the AOI for ERA5; use at least the
        `motion_margin_deg` given to `crop_tropo`.

    """
    output_dir = Path(output_dir)
    if era5:
        files = download_era5_tcwv(
            aoi_bounds, datetimes, output_dir / "era5", margin_deg
        )
        logger.info(f"ERA5: {len(files)} monthly file(s) in {output_dir / 'era5'}")
    if nexrad:
        have = download_nexrad_n0q(datetimes, output_dir / "nexrad")
        logger.info(f"NEXRAD: {len(have)} composite(s) in {output_dir / 'nexrad'}")


def load_era5_tcwv(source: Path | str | Sequence[Path | str]) -> xr.DataArray:
    """Hourly ERA5 TCWV from one file, a list of files, or a directory of them.

    Returns a DataArray on (time, latitude, longitude), sorted in time.
    """
    if isinstance(source, (str, Path)) and Path(source).is_dir():
        files = sorted(Path(source).glob("*.nc"))
    elif isinstance(source, (str, Path)):
        files = [Path(source)]
    else:
        files = [Path(f) for f in source]
    if not files:
        msg = f"No ERA5 files found in {source}"
        raise FileNotFoundError(msg)
    parts = []
    for f in files:
        ds = xr.open_dataset(f)
        if "valid_time" in ds.dims:
            ds = ds.rename({"valid_time": "time"})
        da = ds["tcwv"]
        extra = [d for d in da.dims if d not in ("time", "latitude", "longitude")]
        parts.append(da.isel(dict.fromkeys(extra, 0), drop=True).load())
    out = xr.concat(parts, dim="time").sortby("time")
    return out.isel(time=np.unique(out.time.values, return_index=True)[1])


def _era5_on_grid(
    tcwv: xr.DataArray, t: pd.Timestamp, latitude: np.ndarray, longitude: np.ndarray
) -> np.ndarray:
    """ERA5 TCWV at time `t` (linear between hours), bilinear on the grid."""
    t64 = np.datetime64(pd.Timestamp(t).tz_localize(None).to_datetime64(), "ns")
    times = tcwv.time.values
    if t64 < times.min() or t64 > times.max():
        msg = f"ERA5 does not cover {pd.Timestamp(t)} ({times.min()} to {times.max()})"
        raise ValueError(msg)
    # linear in time between the two bracketing hours (explicit: xarray's datetime
    # interpolation is version-dependent)
    j = int(np.searchsorted(times, t64, side="right"))
    j = min(max(j, 1), times.size - 1)
    ta, tb = times[j - 1], times[j]
    frac = float((t64 - ta) / (tb - ta)) if tb > ta else 0.0
    at_t = (1.0 - frac) * tcwv.isel(time=j - 1, drop=True) + frac * tcwv.isel(
        time=j, drop=True
    )
    # ERA5 is stored north to south; interpolation needs ascending coordinates
    at_t = at_t.sortby("latitude").sortby("longitude")
    out = at_t.interp(latitude=latitude, longitude=longitude).values.astype(float)
    if np.isnan(out).all():
        msg = (
            "ERA5 does not cover the TROPO grid"
            f" (lat {latitude.min():.2f} to {latitude.max():.2f},"
            f" lon {longitude.min():.2f} to {longitude.max():.2f})"
        )
        raise ValueError(msg)
    if np.isnan(out).any():
        logger.warning(
            "ERA5 covers only part of the TROPO grid; the ERA5 term is zero elsewhere"
        )
    return out


# --------------------------------------------------------------------------- guide
def _nan_smooth(field: np.ndarray, sigma_px: float) -> np.ndarray:
    m = np.isfinite(field)
    num = ndimage.gaussian_filter(np.where(m, field, 0.0), sigma_px, mode="nearest")
    den = ndimage.gaussian_filter(m.astype(float), sigma_px, mode="nearest")
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(den > 0.25, num / den, 0.0)


def interp_in_time_guided(
    ds0: xr.Dataset,
    ds1: xr.Dataset,
    t0: pd.Timestamp,
    t1: pd.Timestamp,
    t: pd.Timestamp,
    *,
    era5_tcwv: xr.DataArray,
    radar: dict[datetime, np.ndarray] | None = None,
    weights: dict[str, float] | None = None,
) -> xr.Dataset:
    """Time interpolation of the total delay cube, guided by ERA5 and NEXRAD.

    Drop-in alternative to the linear interpolation used by `crop_tropo`.

    Parameters
    ----------
    ds0, ds1 : xr.Dataset
        Cropped TROPO products with `hydrostatic_delay` and `wet_delay` on
        (time, height, latitude, longitude).
    t0, t1, t : pd.Timestamp
        Times of the two products and the requested time.
    era5_tcwv : xr.DataArray
        Hourly ERA5 total column water vapour on (time, latitude, longitude)
        covering `t0` to `t1` and the grid (see `load_era5_tcwv`).
    radar : dict[datetime, np.ndarray], optional
        Echo fraction (see `read_nexrad_echo_fraction`) on the grid of `ds0` for
        every time in `guide_radar_times(t0, t1, t)`.  Without it, or if any time
        is missing, only the ERA5 term is used, with its own weight.
    weights : dict, optional
        Override of `DEFAULT_GUIDE_WEIGHTS` (keys "era5", "cover", "motion").

    Returns
    -------
    xr.Dataset
        Copy of `ds0` with a `total_delay` variable on (height, latitude, longitude).

    """
    w = float((t - t0) / (t1 - t0))
    lat, lon = ds0.latitude.values, ds0.longitude.values
    hyd0 = ds0.hydrostatic_delay.squeeze("time", drop=True)
    hyd1 = ds1.hydrostatic_delay.squeeze("time", drop=True)
    wet0_da = ds0.wet_delay.squeeze("time", drop=True)
    wet0 = wet0_da.values.astype(float)
    wet1 = ds1.wet_delay.squeeze("time", drop=True).values.astype(float)
    wet_lin = (1.0 - w) * wet0 + w * wet1

    need = guide_radar_times(
        pd.Timestamp(t0).to_pydatetime(),
        pd.Timestamp(t1).to_pydatetime(),
        pd.Timestamp(t).to_pydatetime(),
    )
    use_radar = radar is not None and all(radar.get(x) is not None for x in need)
    wts = dict(DEFAULT_GUIDE_WEIGHTS)
    if not use_radar:
        wts = {"era5": DEFAULT_ERA5_ONLY_WEIGHT, "cover": 0.0, "motion": 0.0}
    if weights:
        wts.update(weights)

    # ERA5: fractional departure of the column from its own straight line
    e0, e1, et = (_era5_on_grid(era5_tcwv, x, lat, lon) for x in (t0, t1, t))
    with np.errstate(invalid="ignore", divide="ignore"):
        rho = et / ((1.0 - w) * e0 + w * e1) - 1.0
    rho = np.nan_to_num(rho, nan=0.0, posinf=0.0, neginf=0.0)
    wet = wet_lin * (1.0 + wts["era5"] * rho)

    if use_radar and (wts["cover"] or wts["motion"]):
        pixel_deg = float(np.abs(np.diff(lat)).mean())
        sig_cover = COVER_SMOOTH_KM / (pixel_deg * KM_PER_DEG)
        sig_track = TRACK_SMOOTH_KM / (pixel_deg * KM_PER_DEG)
        frames = {x: np.nan_to_num(np.asarray(radar[x], float)) for x in need}
        t0p, t1p = pd.Timestamp(t0).to_pydatetime(), pd.Timestamp(t1).to_pydatetime()
        ts = pd.Timestamp(t).round("5min").to_pydatetime()
        # echo-cover change against its own straight line, applied along the wet
        # delay's own height profile (equal to the fitted value at sea level)
        s0 = _nan_smooth(frames[t0p], sig_cover)
        s1 = _nan_smooth(frames[t1p], sig_cover)
        st = _nan_smooth(frames[ts], sig_cover)
        dcover = st - ((1.0 - w) * s0 + w * s1)
        k0 = int(np.abs(wet0_da.height.values - 0.0).argmin())
        profile = wet_lin / np.maximum(wet_lin[k0], 1e-6)
        wet = wet + wts["cover"] * dcover[None] * profile
        # echo motion, tracked hour by hour and summed over the interval
        hourly = guide_radar_times(t0p, t1p, t0p)
        reduce = max(round(TRACK_PIXEL_DEG / pixel_deg), 1)
        flow = np.zeros((2, lat.size, lon.size))
        for a, b in zip(hourly[:-1], hourly[1:]):
            flow += estimate_flow(
                _nan_smooth(frames[a], sig_track),
                _nan_smooth(frames[b], sig_track),
                reduce=reduce,
            )
        moved = (1.0 - w) * warp(wet0, -flow, w) + w * warp(wet1, flow, 1.0 - w)
        wet = wet + wts["motion"] * (moved - wet_lin)

    out = ds0.copy(deep=True)
    total = (1.0 - w) * hyd0 + w * hyd1
    out["total_delay"] = total + xr.DataArray(
        wet.astype(total.dtype), coords=wet0_da.coords, dims=wet0_da.dims
    )
    out["total_delay"].attrs = {
        "long_name": "Zenith total delay",
        "units": "meters",
        "time_interpolation": "guided",
        "guide_sources": "ERA5 + NEXRAD" if use_radar else "ERA5",
        "guide_weight_era5": float(wts["era5"]),
        "guide_weight_cover": float(wts["cover"]),
        "guide_weight_motion": float(wts["motion"]),
    }
    return out
