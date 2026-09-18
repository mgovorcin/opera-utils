"""Motion-aware temporal interpolation of the OPERA TROPO wet delay.

TROPO products exist every 6 hours.  Linear interpolation in time between two
products does not move a weather front: it fades the moist air out at its old
position and in at its new one.  The functions here estimate how the wet delay
field moved between the two products (dense optical flow), slide both fields
along that motion to the requested time, and blend them.

The hydrostatic delay is smooth in time and space and is always interpolated
linearly.

Only ``numpy`` and ``scipy`` are needed.

Notes
-----
* Motion is estimated on the wet delay at one fixed height level (default
  1500 m).  At a fixed height the field carries no terrain imprint, and tests
  showed one motion field for the whole column to be as good as one per level.
* The prediction is blended with plain linear interpolation through
  ``motion_weight``.  Scored against 5-minute GNSS over the contiguous US, the
  pure motion-aware prediction (weight 1) overshoots by about a factor of two
  because sharp features in the model are rarely in exactly the right place; a
  weight near 0.5 reduced the error on every date tested.  Weight 0 is linear
  interpolation.
* A front travels roughly 300 km in 6 hours, so the fields must extend well
  beyond the area of interest.  `crop_tropo` reads a wider margin when motion
  interpolation is requested and trims the result.

"""

from __future__ import annotations

import numpy as np
import pandas as pd
import xarray as xr
from scipy import ndimage

__all__ = [
    "estimate_flow",
    "interp_in_time_motion",
    "interpolate_motion",
    "warp",
]

DEFAULT_FLOW_HEIGHT = 1500.0
DEFAULT_MOTION_WEIGHT = 0.5
# Track motion on pixels of about this size, with a Lucas-Kanade window of
# DEFAULT_FLOW_RADIUS such pixels.  Scored against GNSS, skill rose with the window
# size and levelled off near these values: finer scales of a 6-hourly global model
# field do not move coherently.
DEFAULT_FLOW_PIXEL_DEG = 0.56
DEFAULT_FLOW_RADIUS = 7


def _normalize(a: np.ndarray, b: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Scale both images to [0, 1] with a common, outlier-robust range."""
    both = np.concatenate([a.ravel(), b.ravel()])
    fill = float(np.nanmean(both))
    lo, hi = np.nanpercentile(both, [0.5, 99.5])
    scale = max(float(hi - lo), 1e-12)

    def f(z: np.ndarray) -> np.ndarray:
        return np.clip((np.nan_to_num(z, nan=fill) - lo) / scale, 0.0, 1.0)

    return f(a), f(b)


def _resize(img: np.ndarray, shape: tuple[int, int]) -> np.ndarray:
    """Bilinear resize of a 2D array to `shape`, aligning the image corners."""
    rows = np.linspace(0, img.shape[0] - 1, shape[0])
    cols = np.linspace(0, img.shape[1] - 1, shape[1])
    grid = np.meshgrid(rows, cols, indexing="ij")
    return ndimage.map_coordinates(img, grid, order=1, mode="nearest")


def _pyramid(img: np.ndarray, min_size: int) -> list[np.ndarray]:
    """Gaussian pyramid, finest level first, halving until `min_size`."""
    levels = [img]
    while min(levels[-1].shape) // 2 >= min_size:
        smooth = ndimage.gaussian_filter(levels[-1], sigma=1.0, mode="nearest")
        levels.append(smooth[::2, ::2])
    return levels


def _ilk_level(
    ref: np.ndarray,
    mov: np.ndarray,
    flow: np.ndarray,
    radius: int,
    num_warp: int,
) -> np.ndarray:
    """Refine `flow` at one pyramid level with iterative Lucas-Kanade."""
    grid = np.array(
        np.meshgrid(np.arange(ref.shape[0]), np.arange(ref.shape[1]), indexing="ij"),
        dtype=float,
    )
    sigma = (2 * radius + 1) / 4.0

    def window(z: np.ndarray) -> np.ndarray:
        return ndimage.gaussian_filter(z, sigma=sigma, mode="nearest")

    for _ in range(num_warp):
        moved = ndimage.map_coordinates(mov, grid + flow, order=1, mode="nearest")
        gy, gx = np.gradient(moved)
        # Linearize about the current flow: solve for the total flow, not an update
        err = gy * flow[0] + gx * flow[1] + ref - moved
        a11, a12, a22 = window(gy * gy), window(gy * gx), window(gx * gx)
        b1, b2 = window(gy * err), window(gx * err)
        # Tikhonov term: where the image has no texture the flow goes to zero,
        # which makes the result fall back to linear interpolation there.
        lam = 1e-3 * float(np.mean(a11 + a22)) + 1e-12
        a11, a22 = a11 + lam, a22 + lam
        det = a11 * a22 - a12 * a12
        flow = np.stack([(a22 * b1 - a12 * b2) / det, (a11 * b2 - a12 * b1) / det])
        flow = ndimage.gaussian_filter(
            flow, sigma=(0, sigma / 2, sigma / 2), mode="nearest"
        )
    return flow


def estimate_flow(
    ref: np.ndarray,
    mov: np.ndarray,
    *,
    reduce: int = 1,
    radius: int = DEFAULT_FLOW_RADIUS,
    num_warp: int = 6,
    min_size: int = 8,
) -> np.ndarray:
    """Dense displacement field such that ``ref[y, x] ~ mov[y + v, x + u]``.

    Coarse-to-fine iterative Lucas-Kanade optical flow.

    Parameters
    ----------
    ref, mov : np.ndarray
        2D fields on the same grid.  NaNs are filled with the mean.
    reduce : int
        Estimate the flow on a grid this many times coarser, then resample it
        to the input grid.  Use it to bring the pixel size to the scale at
        which the field moves coherently.
    radius : int
        Half-width, in (reduced) pixels, of the Lucas-Kanade window.
    num_warp : int
        Warping iterations per pyramid level.
    min_size : int
        Smallest image side allowed at the coarsest pyramid level.

    Returns
    -------
    np.ndarray
        Array of shape (2, ny, nx): displacement in pixels of the input grid,
        along rows (``v``) and along columns (``u``).

    """
    if ref.shape != mov.shape or ref.ndim != 2:
        msg = f"ref and mov must be 2D with equal shapes, got {ref.shape}, {mov.shape}"
        raise ValueError(msg)
    full_shape = ref.shape
    a, b = _normalize(np.asarray(ref, dtype=float), np.asarray(mov, dtype=float))
    if reduce > 1:
        small = (max(full_shape[0] // reduce, 2), max(full_shape[1] // reduce, 2))
        a = _resize(ndimage.gaussian_filter(a, reduce / 2.0, mode="nearest"), small)
        b = _resize(ndimage.gaussian_filter(b, reduce / 2.0, mode="nearest"), small)

    pyr_a, pyr_b = _pyramid(a, min_size), _pyramid(b, min_size)
    flow = np.zeros((2, *pyr_a[-1].shape))
    for lev_a, lev_b in zip(reversed(pyr_a), reversed(pyr_b)):
        if flow.shape[1:] != lev_a.shape:
            ratio = (lev_a.shape[0] / flow.shape[1], lev_a.shape[1] / flow.shape[2])
            flow = np.stack(
                [
                    _resize(flow[0], lev_a.shape) * ratio[0],
                    _resize(flow[1], lev_a.shape) * ratio[1],
                ]
            )
        flow = _ilk_level(lev_a, lev_b, flow, radius, num_warp)

    if flow.shape[1:] != full_shape:
        ratio = (full_shape[0] / flow.shape[1], full_shape[1] / flow.shape[2])
        flow = np.stack(
            [
                _resize(flow[0], full_shape) * ratio[0],
                _resize(flow[1], full_shape) * ratio[1],
            ]
        )
    return flow


def warp(img: np.ndarray, flow: np.ndarray, fraction: float) -> np.ndarray:
    """Sample `img` at ``x + fraction * flow(x)``.

    `img` may be 2D or 3D with the two trailing axes matching `flow`; the same
    displacement is then applied to every leading slice.
    """
    ny, nx = flow.shape[1:]
    grid = np.array(
        np.meshgrid(np.arange(ny), np.arange(nx), indexing="ij"), dtype=float
    )
    coords = grid + fraction * flow
    if img.ndim == 2:
        return ndimage.map_coordinates(img, coords, order=1, mode="nearest")
    out = np.empty(img.shape, dtype=float)
    for k in range(img.shape[0]):
        out[k] = ndimage.map_coordinates(img[k], coords, order=1, mode="nearest")
    return out


def interpolate_motion(
    field0: np.ndarray,
    field1: np.ndarray,
    weight: float,
    *,
    track0: np.ndarray | None = None,
    track1: np.ndarray | None = None,
    motion_weight: float = DEFAULT_MOTION_WEIGHT,
    reduce: int = 1,
) -> np.ndarray:
    """Interpolate between two fields along their estimated motion.

    Parameters
    ----------
    field0, field1 : np.ndarray
        Fields at the earlier and later time, 2D or 3D (levels first).
    weight : float
        Position of the requested time, 0 at `field0` and 1 at `field1`.
    track0, track1 : np.ndarray, optional
        2D fields used to estimate the motion.  Default: the inputs themselves,
        which must then be 2D.
    motion_weight : float
        0 gives linear interpolation, 1 the pure motion-aware result.
    reduce : int
        Passed to `estimate_flow`.

    Returns
    -------
    np.ndarray
        The interpolated field, same shape as the inputs.

    """
    if not 0.0 <= motion_weight <= 1.0:
        msg = f"motion_weight must be in [0, 1], got {motion_weight}"
        raise ValueError(msg)
    linear = (1.0 - weight) * field0 + weight * field1
    if motion_weight == 0.0 or weight <= 0.0 or weight >= 1.0:
        return linear
    if track0 is None or track1 is None:
        if field0.ndim != 2:
            msg = "track0 and track1 are required for 3D fields"
            raise ValueError(msg)
        track0, track1 = field0, field1

    flow01 = estimate_flow(track0, track1, reduce=reduce)
    flow10 = estimate_flow(track1, track0, reduce=reduce)
    # A pixel at the requested time came from x + w * F10 in field0 and will be
    # at x + (1 - w) * F01 in field1.
    moved = (1.0 - weight) * warp(field0, flow10, weight) + weight * warp(
        field1, flow01, 1.0 - weight
    )
    return (1.0 - motion_weight) * linear + motion_weight * moved


def interp_in_time_motion(
    ds0: xr.Dataset,
    ds1: xr.Dataset,
    t0: pd.Timestamp,
    t1: pd.Timestamp,
    t: pd.Timestamp,
    *,
    motion_weight: float = DEFAULT_MOTION_WEIGHT,
    flow_height: float = DEFAULT_FLOW_HEIGHT,
) -> xr.Dataset:
    """Time interpolation of the total delay cube, motion-aware for the wet part.

    Drop-in alternative to the linear interpolation used by `crop_tropo`.  The
    hydrostatic delay is interpolated linearly; the wet delay with
    `interpolate_motion`, using the motion of the wet delay at the height level
    nearest to `flow_height`.

    Parameters
    ----------
    ds0, ds1 : xr.Dataset
        Cropped TROPO products with `hydrostatic_delay` and `wet_delay` on
        (time, height, latitude, longitude).
    t0, t1, t : pd.Timestamp
        Times of the two products and the requested time.
    motion_weight : float
        0 gives linear interpolation, 1 the pure motion-aware result.
    flow_height : float
        Height, in meters, of the level on which motion is estimated.

    Returns
    -------
    xr.Dataset
        Copy of `ds0` with a `total_delay` variable on (height, latitude, longitude).

    """
    w = float((t - t0) / (t1 - t0))
    hyd0 = ds0.hydrostatic_delay.squeeze("time", drop=True)
    hyd1 = ds1.hydrostatic_delay.squeeze("time", drop=True)
    wet0 = ds0.wet_delay.squeeze("time", drop=True)
    wet1 = ds1.wet_delay.squeeze("time", drop=True)

    k = int(np.abs(wet0.height.values - flow_height).argmin())
    pixel_deg = float(np.abs(np.diff(wet0.latitude.values)).mean())
    reduce = max(round(DEFAULT_FLOW_PIXEL_DEG / pixel_deg), 1)
    wet = interpolate_motion(
        wet0.values.astype(float),
        wet1.values.astype(float),
        w,
        track0=wet0.values[k],
        track1=wet1.values[k],
        motion_weight=motion_weight,
        reduce=reduce,
    )
    out = ds0.copy(deep=True)
    total = (1.0 - w) * hyd0 + w * hyd1
    out["total_delay"] = total + xr.DataArray(
        wet.astype(total.dtype), coords=wet0.coords, dims=wet0.dims
    )
    # Arithmetic may carry over the description of the hydrostatic input
    out["total_delay"].attrs = {
        "long_name": "Zenith total delay",
        "units": "meters",
        "time_interpolation": "motion",
        "motion_weight": float(motion_weight),
        "flow_height": float(wet0.height.values[k]),
    }
    return out
