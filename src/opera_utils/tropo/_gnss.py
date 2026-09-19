"""GNSS-guided refinement of the OPERA TROPO correction.

Where a dense GNSS network exists, the zenith total delay (ZTD) it measures every
few minutes can correct what the weather model gets wrong.  The scheme here keeps
the model correction and adds a kriged field of (GNSS - model) residuals:

1. for every station, the model ZTD is read from the cropped TROPO cube at the
   station's own position and height;
2. the residual GNSS - model is differenced between a date and the reference
   date, using only stations present on both, so that constant per-station biases
   (antenna, height datum) cancel;
3. the differenced residuals are interpolated to the output grid with zero-mean
   kriging, whose range and nugget are chosen by leave-one-out at the stations;
4. the field is added to the model zenith delay before projection to line of sight.

Far from any station the kriged field decays to zero, so the result falls back to
the model correction.  In a test over the Los Angeles basin (192 stations, 7 km
spacing) this halved the residual of a Sentinel-1 interferogram relative to the
model alone; with stations more than about 15 km apart the benefit disappeared and
the result equalled the model.  It cannot be used without a reference date.

GNSS input is a table with columns ``id, lat, lon, height, datetime, ztd`` (meters),
one row per station and acquisition time.  `download_ngl_ztd` builds one from the
5-minute products of the Nevada Geodetic Laboratory (Blewitt et al., 2018).

Only ``numpy``, ``pandas``, ``scipy`` and ``xarray`` are needed for the correction.
"""

from __future__ import annotations

import gzip
import io
import logging
import re
import urllib.request
import zipfile
from collections.abc import Sequence
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
from scipy.interpolate import RegularGridInterpolator
from scipy.spatial.distance import cdist

logger = logging.getLogger(__name__)

__all__ = [
    "download_ngl_ztd",
    "gnss_residual_field",
    "krige",
    "load_gnss_ztd",
    "model_ztd_at_stations",
    "read_sinex_tropo",
    "station_residual_difference",
    "tune_kriging",
]

NGL_BASE = "https://geodesy.unr.edu/gps_timeseries"
NGL_HOLDINGS = "https://geodesy.unr.edu/NGLStationPages/DataHoldings.txt"
KRIGING_RANGES_KM = (8.0, 15.0, 25.0, 40.0, 70.0, 120.0)
KRIGING_NUGGETS = (0.02, 0.05, 0.1, 0.2, 0.4)
MIN_STATIONS = 5
GNSS_COLUMNS = ("id", "lat", "lon", "height", "datetime", "ztd")


# --------------------------------------------------------------------------- #
# kriging
# --------------------------------------------------------------------------- #
def _to_km(lat: np.ndarray, lon: np.ndarray, lat0: float) -> np.ndarray:
    """Project to a local plane in km (equirectangular about `lat0`)."""
    return np.column_stack(
        [
            np.asarray(lat, dtype=float) * 111.19,
            np.asarray(lon, dtype=float) * 111.19 * np.cos(np.radians(lat0)),
        ]
    )


def _covariance(a: np.ndarray, b: np.ndarray, range_km: float) -> np.ndarray:
    return np.exp(-cdist(a, b) / range_km)


def krige(
    xy_km: np.ndarray,
    values: np.ndarray,
    xy_target_km: np.ndarray,
    range_km: float,
    nugget: float,
    chunk: int = 50_000,
) -> np.ndarray:
    """Zero-mean (simple) kriging with an exponential covariance.

    The prediction decays to zero away from the data, which is what makes the
    GNSS-guided correction fall back to the model where there are no stations.

    Parameters
    ----------
    xy_km : np.ndarray
        Station coordinates, shape (n, 2), in km.
    values : np.ndarray
        Station values, shape (n,).  Remove their mean first if it is not
        meaningful.
    xy_target_km : np.ndarray
        Target coordinates, shape (m, 2), in km.
    range_km : float
        e-folding distance of the covariance.
    nugget : float
        Uncorrelated variance as a fraction of the signal variance.
    chunk : int
        Targets processed at once.

    Returns
    -------
    np.ndarray
        Predictions, shape (m,).

    """
    cov = _covariance(xy_km, xy_km, range_km) + nugget * np.eye(len(values))
    weights_rhs = np.linalg.solve(cov, np.asarray(values, dtype=float))
    out = np.empty(len(xy_target_km), dtype=float)
    for i in range(0, len(xy_target_km), chunk):
        out[i : i + chunk] = (
            _covariance(xy_target_km[i : i + chunk], xy_km, range_km) @ weights_rhs
        )
    return out


def _leave_one_out(
    xy_km: np.ndarray, values: np.ndarray, range_km: float, nugget: float
) -> np.ndarray:
    """Leave-one-out kriging predictions at the stations, in closed form."""
    cov = _covariance(xy_km, xy_km, range_km) + nugget * np.eye(len(values))
    inv = np.linalg.inv(cov)
    return values - (inv @ values) / np.diag(inv)


def tune_kriging(
    xy_km: np.ndarray,
    values: np.ndarray,
    ranges_km: Sequence[float] = KRIGING_RANGES_KM,
    nuggets: Sequence[float] = KRIGING_NUGGETS,
) -> tuple[float, float, float]:
    """Choose range and nugget by leave-one-out.

    Returns
    -------
    tuple[float, float, float]
        (range_km, nugget, leave-one-out RMS error).

    """
    best = (np.inf, ranges_km[0], nuggets[0])
    for rng in ranges_km:
        for nug in nuggets:
            err = values - _leave_one_out(xy_km, values, rng, nug)
            score = float(np.sqrt(np.mean(err**2)))
            if score < best[0]:
                best = (score, rng, nug)
    return best[1], best[2], best[0]


# --------------------------------------------------------------------------- #
# the correction
# --------------------------------------------------------------------------- #
def load_gnss_ztd(filename: str | Path) -> pd.DataFrame:
    """Read a GNSS ZTD table (CSV).

    Columns: `id, lat, lon, height, datetime, ztd`, with `ztd` in meters.
    """
    df = pd.read_csv(filename, parse_dates=["datetime"])
    missing = set(GNSS_COLUMNS) - set(df.columns)
    if missing:
        msg = f"{filename} is missing columns {sorted(missing)}"
        raise ValueError(msg)
    df["datetime"] = pd.to_datetime(df["datetime"], utc=True).dt.tz_localize(None)
    return df


def model_ztd_at_stations(
    total_delay: xr.DataArray, lat: np.ndarray, lon: np.ndarray, height: np.ndarray
) -> np.ndarray:
    """Model zenith total delay at station positions and heights.

    `total_delay` has dims (height, latitude, longitude), as written by
    `crop_tropo`.  Stations outside the cube, or above it, give NaN; stations
    below its lowest level are evaluated at that level.
    """
    lats = total_delay.latitude.values
    order = np.argsort(lats)
    rgi = RegularGridInterpolator(
        (total_delay.height.values, lats[order], total_delay.longitude.values),
        total_delay.values[:, order, :],
        bounds_error=False,
        fill_value=np.nan,
    )
    h = np.maximum(np.asarray(height, dtype=float), float(total_delay.height.values[0]))
    return rgi(np.column_stack([h, np.asarray(lat, float), np.asarray(lon, float)]))


def _rows_at(
    gnss: pd.DataFrame, when: datetime, tolerance: pd.Timedelta
) -> pd.DataFrame:
    dt = (gnss["datetime"] - pd.Timestamp(when)).abs()
    sel = gnss[dt <= tolerance].assign(_dt=dt[dt <= tolerance])
    return sel.sort_values("_dt").drop_duplicates("id").set_index("id")


def station_residual_difference(
    gnss: pd.DataFrame,
    total_delay: xr.DataArray,
    when: datetime,
    total_delay_ref: xr.DataArray,
    when_ref: datetime,
    tolerance_minutes: float = 10.0,
    outlier_sigma: float = 4.0,
) -> pd.DataFrame:
    """(GNSS - model) at `when` minus (GNSS - model) at `when_ref`, per station.

    Only stations with GNSS on both dates and inside both model cubes are kept;
    constant per-station biases cancel in the difference.  Values further than
    `outlier_sigma` robust standard deviations from the median are dropped.

    Returns
    -------
    pd.DataFrame
        Index `id`; columns `lat, lon, height, residual` (meters of zenith delay).

    """
    tol = pd.Timedelta(minutes=tolerance_minutes)
    a, b = _rows_at(gnss, when, tol), _rows_at(gnss, when_ref, tol)
    common = a.index.intersection(b.index)
    a, b = a.loc[common], b.loc[common]
    res = (
        a.ztd.values
        - model_ztd_at_stations(total_delay, a.lat, a.lon, a.height)
        - (
            b.ztd.values
            - model_ztd_at_stations(total_delay_ref, b.lat, b.lon, b.height)
        )
    )
    out = pd.DataFrame(
        {
            "lat": a.lat.values,
            "lon": a.lon.values,
            "height": a.height.values,
            "residual": res,
        },
        index=common,
    )
    out = out[np.isfinite(out.residual)]
    if len(out) >= MIN_STATIONS:
        med = out.residual.median()
        mad = 1.4826 * (out.residual - med).abs().median()
        if mad > 0:
            out = out[(out.residual - med).abs() <= outlier_sigma * mad]
    return out


def gnss_residual_field(
    residuals: pd.DataFrame, lat: np.ndarray, lon: np.ndarray
) -> tuple[np.ndarray, dict]:
    """Krige station residuals to target points.

    The median residual is removed first: a constant is invisible to InSAR and
    must not leak into areas without stations.

    Parameters
    ----------
    residuals : pd.DataFrame
        Output of `station_residual_difference`.
    lat, lon : np.ndarray
        Target coordinates in degrees, any (equal) shape.

    Returns
    -------
    field : np.ndarray
        Residual zenith delay at the targets [m], same shape as `lat`.  Zeros
        if fewer than 5 stations are available.
    info : dict
        `n_stations`, `range_km`, `nugget`, `loo_rms`, `residual_rms`.

    """
    lat, lon = np.asarray(lat, dtype=float), np.asarray(lon, dtype=float)
    info = {
        "n_stations": len(residuals),
        "range_km": np.nan,
        "nugget": np.nan,
        "loo_rms": np.nan,
        "residual_rms": np.nan,
    }
    if len(residuals) < MIN_STATIONS:
        logger.warning(
            f"Only {len(residuals)} GNSS stations: using the model correction alone"
        )
        return np.zeros(lat.shape, dtype=float), info
    lat0 = float(residuals.lat.mean())
    xy = _to_km(residuals.lat.values, residuals.lon.values, lat0)
    vals = residuals.residual.values - np.median(residuals.residual.values)
    rng, nug, loo = tune_kriging(xy, vals)
    info.update(
        range_km=rng,
        nugget=nug,
        loo_rms=loo,
        residual_rms=float(np.sqrt(np.mean(vals**2))),
    )
    ok = np.isfinite(lat) & np.isfinite(lon)
    field = np.zeros(lat.shape, dtype=float)
    field[ok] = krige(xy, vals, _to_km(lat[ok], lon[ok], lat0), rng, nug)
    return field, info


# --------------------------------------------------------------------------- #
# Nevada Geodetic Laboratory 5-minute tropospheric products
# --------------------------------------------------------------------------- #
def read_sinex_tropo(raw: bytes) -> pd.DataFrame:
    """Parse a (gzipped) SINEX_TROPO file into `datetime, ztd, ztd_sigma` [m]."""
    if raw[:2] == b"\x1f\x8b":
        raw = gzip.decompress(raw)
    text = raw.decode("ascii", errors="replace")
    block = re.search(r"\+TROP/SOLUTION(.*?)-TROP/SOLUTION", text, re.DOTALL)
    rows = []
    if block:
        for line in block.group(1).splitlines():
            parts = line.split()
            if line.startswith("*") or len(parts) < 4:
                continue
            m = re.fullmatch(r"(\d{2,4}):(\d{3}):(\d+)", parts[1])
            if not m:
                continue
            yy, doy, sec = (int(g) for g in m.groups())
            year = yy if yy > 100 else (2000 + yy if yy < 80 else 1900 + yy)
            try:
                ztd, sig = float(parts[2]) / 1000.0, float(parts[3]) / 1000.0
            except ValueError:
                continue
            when = pd.Timestamp(year=year, month=1, day=1) + pd.Timedelta(
                days=doy - 1, seconds=sec
            )
            rows.append((when, ztd, sig))
    return pd.DataFrame(rows, columns=["datetime", "ztd", "ztd_sigma"])


def _fetch(url: str, timeout: float = 120.0) -> bytes | None:
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "opera-utils"})
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return r.read()
    except Exception:
        return None


def _ngl_stations(
    bounds: tuple[float, float, float, float], start: datetime, end: datetime
) -> pd.DataFrame:
    raw = _fetch(NGL_HOLDINGS)
    if raw is None:
        msg = f"Could not download {NGL_HOLDINGS}"
        raise ConnectionError(msg)
    names = ["id", "lat", "lon", "height", "x", "y", "z", "begin", "end"]
    df = pd.read_csv(
        io.BytesIO(raw), sep=r"\s+", usecols=range(9), names=names, skiprows=1
    )
    df["lon"] = ((df.lon + 180.0) % 360.0) - 180.0
    west, south, east, north = bounds
    keep = (
        df.lat.between(south, north)
        & df.lon.between(west, east)
        & (df.begin <= f"{start:%Y-%m-%d}")
        & (df.end >= f"{end:%Y-%m-%d}")
    )
    return df.loc[keep, ["id", "lat", "lon", "height"]].reset_index(drop=True)


def _ngl_year_zip(station: str, year: int, cache_dir: Path) -> Path | None:
    dest = cache_dir / f"{station}.{year}.trop.zip"
    if dest.exists():
        return dest if dest.stat().st_size > 1000 else None
    for sub in ("IGS20/trop", "trop"):
        raw = _fetch(f"{NGL_BASE}/{sub}/{station}/{station}.{year}.trop.zip")
        if raw is not None and len(raw) > 1000:
            dest.write_bytes(raw)
            return dest
    dest.write_bytes(b"")  # remember the miss
    return None


def download_ngl_ztd(
    aoi_bounds: tuple[float, float, float, float],
    datetimes: list[datetime],
    output_file: Path = Path("gnss_ztd.csv"),
    cache_dir: Path = Path("ngl_tropo_cache"),
    window_minutes: float = 8.0,
    margin_deg: float = 0.25,
    num_workers: int = 8,
) -> Path:
    """Build a GNSS ZTD table from Nevada Geodetic Laboratory 5-minute products.

    Parameters
    ----------
    aoi_bounds : tuple[float, float, float, float]
        (west, south, east, north) in degrees.
    datetimes : list[datetime]
        Acquisition times (UTC).  ZTD is averaged over +/- `window_minutes`.
    output_file : Path
        CSV to write, with columns `id, lat, lon, height, datetime, ztd`.
    cache_dir : Path
        Where station-year archives (about 2.4 MB each) are kept.
    window_minutes : float
        Half-width of the averaging window.
    margin_deg : float
        Stations this far outside the AOI are included.
    num_workers : int
        Parallel downloads.

    Returns
    -------
    Path
        `output_file`.

    """
    cache_dir.mkdir(parents=True, exist_ok=True)
    times = sorted(
        (
            pd.Timestamp(d).tz_localize(None)
            if pd.Timestamp(d).tzinfo is None
            else pd.Timestamp(d).tz_convert("UTC").tz_localize(None)
        )
        for d in datetimes
    )
    west, south, east, north = aoi_bounds
    bounds = (
        west - margin_deg,
        south - margin_deg,
        east + margin_deg,
        north + margin_deg,
    )
    stations = _ngl_stations(bounds, times[0], times[-1])
    years = sorted({t.year for t in times})
    jobs = [(s, y) for s in stations.id for y in years]
    logger.info(f"{len(stations)} stations, {len(jobs)} station-year archives")
    with ThreadPoolExecutor(max_workers=num_workers) as ex:
        zips = dict(
            zip(jobs, ex.map(lambda j: _ngl_year_zip(j[0], j[1], cache_dir), jobs))
        )

    win = pd.Timedelta(minutes=window_minutes)
    rows = []
    for st in stations.itertuples():
        for t in times:
            zp = zips.get((st.id, t.year))
            if zp is None:
                continue
            try:
                with zipfile.ZipFile(zp) as zf:
                    raw = zf.read(f"{st.id}.{t.year}.{t.dayofyear:03d}.trop.gz")
            except (KeyError, zipfile.BadZipFile):
                continue
            sol = read_sinex_tropo(raw)
            near = sol[(sol.datetime >= t - win) & (sol.datetime <= t + win)]
            if len(near) >= 2:
                rows.append((st.id, st.lat, st.lon, st.height, t, near.ztd.mean()))
    table = pd.DataFrame(rows, columns=list(GNSS_COLUMNS))
    table.to_csv(output_file, index=False)
    logger.info(f"{table.id.nunique()} stations, {len(table)} rows -> {output_file}")
    return output_file
