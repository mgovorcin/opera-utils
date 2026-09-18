"""Crop OPERA TROPO products for an area of interest and interpolate them in time."""

from __future__ import annotations

import logging
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path
from typing import Literal

import pandas as pd
import rasterio as rio
import xarray as xr
from rasterio.warp import transform_bounds
from tqdm import tqdm

from ._helpers import (
    MissingTropoError,
    _bracket,
    _build_tropo_index,
    _create_total_delay,
    _interp_in_time,
    _open_crop,
)
from ._motion import DEFAULT_MOTION_WEIGHT, interp_in_time_motion

DEFAULT_MOTION_MARGIN_DEG = 5.0

logger = logging.getLogger(__name__)


def _process_one_datetime(
    dt: datetime,
    tropo_idx_series: pd.Series,
    lat_bounds: tuple[float, float],
    lon_bounds: tuple[float, float],
    height_max: float,
    output_dir: Path,
    skip_time_interpolation: bool,
    debug: bool = False,
    time_interpolation: Literal["linear", "motion"] = "linear",
    motion_weight: float = DEFAULT_MOTION_WEIGHT,
    motion_margin_deg: float = DEFAULT_MOTION_MARGIN_DEG,
) -> tuple[datetime, str]:
    """Worker: process one datetime and write output file.

    Parameters
    ----------
    dt : datetime
        Datetime to process.
    tropo_idx_series : pd.Series
        Series of TROPO product URLs/paths.
    lat_bounds : tuple[float, float]
        Latitude bounds as (north, south) in degrees.
    lon_bounds : tuple[float, float]
        Longitude bounds as (west, east) in degrees.
    height_max : float
        Maximum height in meters to include in cropping.
        Higher values with smaller atmospheric delay are ignored.
    output_dir : Path
        Directory to save cropped TROPO products.
    skip_time_interpolation : bool
        Skip time interpolation and use nearest file.
    debug : bool
        Debug mode. If True, write debug info during processing.
        Default: False.
    time_interpolation : {"linear", "motion"}
        How to interpolate between the two bracketing products.
    motion_weight : float
        Weight of the motion-aware prediction when `time_interpolation="motion"`.
    motion_margin_deg : float
        Extra margin, in degrees, read around the bounds to track the motion.

    Returns
    -------
    tuple[datetime, str]
        (datetime, status)
        status is "ok" | "skipped" | "missing" | "error:<msg>"

    """
    dt_pandas = pd.to_datetime(dt).tz_localize(None)
    time_str = dt_pandas.strftime("%Y%m%dT%H%M%S")
    output_file = output_dir / f"tropo_cropped_{time_str}.nc"

    if output_file.exists():
        return (dt, "skipped")

    if not skip_time_interpolation:
        try:
            early_url, late_url = _bracket(tropo_idx_series, dt_pandas)
        except MissingTropoError:
            return (dt, "missing")

        read_lat, read_lon = lat_bounds, lon_bounds
        if time_interpolation == "motion":
            # Weather moves ~300 km in 6 h: track it on a wider area, trim below
            read_lat = (
                min(lat_bounds[0] + motion_margin_deg, 90.0),
                max(lat_bounds[1] - motion_margin_deg, -90.0),
            )
            read_lon = (
                lon_bounds[0] - motion_margin_deg,
                lon_bounds[1] + motion_margin_deg,
            )

        if debug:
            tqdm.write(f"Cropping {early_url}")
        ds0 = _open_crop(early_url, read_lat, read_lon, height_max)
        if debug:
            tqdm.write(f"Cropping {late_url}")
        ds1 = _open_crop(late_url, read_lat, read_lon, height_max)

        t0 = ds0.time.to_pandas().item()
        t1 = ds1.time.to_pandas().item()
        if time_interpolation == "motion":
            td_interp = interp_in_time_motion(
                ds0, ds1, t0, t1, dt_pandas, motion_weight=motion_weight
            ).sel(
                latitude=slice(lat_bounds[0], lat_bounds[1]),
                longitude=slice(lon_bounds[0], lon_bounds[1]),
            )
        else:
            td_interp = _interp_in_time(ds0, ds1, t0, t1, dt_pandas)
    else:
        idx = tropo_idx_series.index.get_indexer([dt_pandas], method="nearest")[0]
        closest_url = tropo_idx_series.values[idx]
        ds = _open_crop(closest_url, lat_bounds, lon_bounds, height_max)
        da_total_delay = _create_total_delay(ds)
        td_interp = ds.copy()
        td_interp["total_delay"] = da_total_delay

    # Keep only total_delay and coord variables for output
    output_ds = xr.Dataset(
        {
            "total_delay": td_interp.total_delay,
            "latitude": td_interp.latitude,
            "longitude": td_interp.longitude,
            "height": td_interp.height,
        }
    )
    output_ds.to_netcdf(output_file, engine="h5netcdf")
    return (dt, "ok")


def crop_tropo(
    tropo_urls_file: Path,
    datetimes: list[datetime],
    aoi_bounds: tuple[float, float, float, float] | None = None,
    file_bounds: Path | str | None = None,
    output_dir: Path = Path("cropped_tropo"),
    skip_time_interpolation: bool = False,
    height_max: float = 10000.0,
    margin_deg: float = 0.3,
    num_workers: int = 2,
    time_interpolation: Literal["linear", "motion"] = "linear",
    motion_weight: float = DEFAULT_MOTION_WEIGHT,
    motion_margin_deg: float = DEFAULT_MOTION_MARGIN_DEG,
) -> None:
    """Crop OPERA TROPO products to AOI and interpolate to specific datetimes.

    Parameters
    ----------
    tropo_urls_file : Path
        File containing list of TROPO product URLs/paths (one per line).
    datetimes : list[datetime]
        List of datetime objects to get corrections for.
    aoi_bounds : tuple[float, float, float, float]
        AOI bounding box as (west, south, east, north) in degrees.
    file_bounds : Path | str | None
        Path to GeoTIFF file containing bounds to crop to (alternative to `aoi_bounds`).
    output_dir : Path
        Directory to save cropped TROPO products.
    skip_time_interpolation : bool
        Skip time interpolation and use nearest file.
    height_max : float
        Maximum height in meters to include in cropping.
    margin_deg : float
        Additional margin in degrees around AOI bounds.
    num_workers : int
        Processes to use. Default: 2
    time_interpolation : {"linear", "motion"}
        How to interpolate between the two products bracketing each datetime.
        "linear" (default) blends them pixel by pixel.  "motion" estimates how
        the wet delay field moved between them and interpolates along that
        motion, so that weather fronts are displaced instead of faded out; the
        hydrostatic delay is still interpolated linearly.
        See `opera_utils.tropo.interp_in_time_motion`.
    motion_weight : float
        Only for "motion": weight in [0, 1] of the motion-aware prediction
        against the linear one.  0 is linear interpolation.
    motion_margin_deg : float
        Only for "motion": additional margin in degrees, on top of `margin_deg`,
        read around the AOI to track the motion.  The output is trimmed back to
        the AOI plus `margin_deg`.  Margins below about 3 degrees leave too little
        context and give no benefit over linear interpolation.

    """
    if time_interpolation not in ("linear", "motion"):
        msg = (
            f"time_interpolation must be 'linear' or 'motion', got {time_interpolation}"
        )
        raise ValueError(msg)
    if time_interpolation == "motion" and skip_time_interpolation:
        msg = "time_interpolation='motion' cannot be used with skip_time_interpolation"
        raise ValueError(msg)
    if not 0.0 <= motion_weight <= 1.0:
        msg = f"motion_weight must be in [0, 1], got {motion_weight}"
        raise ValueError(msg)

    output_dir.mkdir(exist_ok=True, parents=True)

    if file_bounds is not None:
        if aoi_bounds is not None:
            msg = "Cannot specify both aoi_bounds and file_bounds"
            raise ValueError(msg)

        with rio.open(file_bounds) as src:
            aoi_bounds = src.bounds
            if src.crs != "epsg:4326":
                aoi_bounds = transform_bounds(src.crs, "epsg:4326", *aoi_bounds)

    assert aoi_bounds is not None
    west, south, east, north = aoi_bounds
    # For lat, use north -> south ordering for xarray slicing
    lat_bounds = (north + margin_deg, south - margin_deg)
    lon_bounds = (west - margin_deg, east + margin_deg)

    tropo_urls = Path(tropo_urls_file).read_text(encoding="utf-8").splitlines()
    tropo_idx_series = _build_tropo_index(tropo_urls)

    if num_workers < 1:
        msg = "num_workers must be >= 1"
        raise ValueError(msg)

    logger.info(f"Processing {len(datetimes)} datetime(s) with {num_workers} worker(s)")

    # Submit all tasks up front; tqdm tracks completions.
    futures = []
    status_counts: dict[str, int] = {"ok": 0, "skipped": 0, "missing": 0, "error": 0}

    if num_workers > 1:
        with ProcessPoolExecutor(max_workers=num_workers) as ex:
            for dt in datetimes:
                futures.append(
                    ex.submit(
                        _process_one_datetime,
                        dt,
                        tropo_idx_series,
                        lat_bounds,
                        lon_bounds,
                        height_max,
                        output_dir,
                        skip_time_interpolation,
                        False,
                        time_interpolation,
                        motion_weight,
                        motion_margin_deg,
                    )
                )

            for fut in tqdm(
                as_completed(futures),
                total=len(futures),
                desc="TROPO crop+interp",
                unit="scene",
            ):
                dt, status = fut.result()
                if status.startswith("error:"):
                    status_counts["error"] += 1
                    logger.error(f"{dt}: {status}")
                else:
                    status_counts[status] = status_counts.get(status, 0) + 1

    else:
        for dt in datetimes:
            _process_one_datetime(
                dt,
                tropo_idx_series,
                lat_bounds,
                lon_bounds,
                height_max,
                output_dir,
                skip_time_interpolation,
                time_interpolation=time_interpolation,
                motion_weight=motion_weight,
                motion_margin_deg=motion_margin_deg,
            )

    logger.info(
        "Done. "
        f"ok={status_counts['ok']}, "
        f"skipped={status_counts['skipped']}, "
        f"missing={status_counts['missing']}, "
        f"errors={status_counts['error']}"
    )
