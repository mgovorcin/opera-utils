"""Tests for opera_utils.tropo._motion and its use in crop_tropo."""

from __future__ import annotations

from datetime import datetime

import numpy as np
import pandas as pd
import pytest
import xarray as xr
from scipy import ndimage

from opera_utils.tropo import _crop
from opera_utils.tropo._helpers import _interp_in_time
from opera_utils.tropo._motion import (
    estimate_flow,
    interp_in_time_motion,
    interpolate_motion,
    warp,
)

INTERIOR = slice(24, -24)


def _smooth_field(shape=(120, 160), sigma=6.0, seed=0):
    rng = np.random.default_rng(seed)
    big = ndimage.gaussian_filter(
        rng.standard_normal((shape[0] + 80, shape[1] + 80)), sigma
    )
    return big / big.std()


def _moving_pair(shift_rows, shift_cols, shape=(120, 160)):
    """Fields at times 0, 1/2 and 1 of a pattern translating by (rows, cols)."""
    big = _smooth_field(shape)

    def crop(z):
        return z[40:-40, 40:-40]

    def at(f):
        return crop(
            ndimage.shift(
                big, (f * shift_rows, f * shift_cols), order=3, mode="nearest"
            )
        )

    return at(0.0), at(0.5), at(1.0)


def _rmse(a, b):
    return float(np.sqrt(np.mean((a - b)[INTERIOR, INTERIOR] ** 2)))


class TestEstimateFlow:
    @pytest.mark.parametrize(("dr", "dc"), [(3.0, 6.0), (-4.0, 2.0), (0.0, -10.0)])
    def test_recovers_uniform_shift(self, dr, dc):
        f0, _, f1 = _moving_pair(dr, dc)
        flow = estimate_flow(f0, f1, radius=5)
        # ref[y, x] ~ mov[y + v, x + u]: the pattern moved by (+dr, +dc)
        assert np.median(flow[0][INTERIOR, INTERIOR]) == pytest.approx(dr, abs=0.3)
        assert np.median(flow[1][INTERIOR, INTERIOR]) == pytest.approx(dc, abs=0.3)

    def test_large_shift_needs_reduced_grid(self):
        f0, _, f1 = _moving_pair(0.0, 28.0, shape=(160, 240))
        flow = estimate_flow(f0, f1, reduce=4, radius=5)
        assert flow.shape == (2, *f0.shape)
        assert np.median(flow[1][INTERIOR, INTERIOR]) == pytest.approx(28.0, abs=1.5)

    def test_no_motion_gives_zero_flow(self):
        f0 = _moving_pair(0, 0)[0]
        flow = estimate_flow(f0, f0.copy())
        assert np.abs(flow).max() < 1e-6

    def test_flat_field_gives_zero_flow(self):
        flat = np.full((64, 64), 3.0)
        assert np.abs(estimate_flow(flat, flat + 1.0)).max() < 1e-6

    def test_nans_are_tolerated(self):
        f0, _, f1 = _moving_pair(2.0, 4.0)
        f0 = f0.copy()
        f0[:5, :5] = np.nan
        assert np.isfinite(estimate_flow(f0, f1)).all()

    def test_shape_mismatch_raises(self):
        with pytest.raises(ValueError, match="equal shapes"):
            estimate_flow(np.zeros((10, 10)), np.zeros((10, 12)))


class TestInterpolateMotion:
    def test_endpoints_are_returned_unchanged(self):
        f0, _, f1 = _moving_pair(3, 6)
        np.testing.assert_array_equal(interpolate_motion(f0, f1, 0.0), f0)
        np.testing.assert_array_equal(interpolate_motion(f0, f1, 1.0), f1)

    def test_zero_motion_weight_is_linear(self):
        f0, _, f1 = _moving_pair(3, 6)
        out = interpolate_motion(f0, f1, 0.25, motion_weight=0.0)
        np.testing.assert_allclose(out, 0.75 * f0 + 0.25 * f1)

    def test_beats_linear_on_a_translating_field(self):
        f0, fmid, f1 = _moving_pair(4.0, 8.0)
        linear = _rmse(0.5 * (f0 + f1), fmid)
        moved = _rmse(interpolate_motion(f0, f1, 0.5, motion_weight=1.0), fmid)
        assert moved < linear / 5

    def test_motion_weight_blends_the_two(self):
        f0, _, f1 = _moving_pair(4.0, 8.0)
        lin = interpolate_motion(f0, f1, 0.5, motion_weight=0.0)
        full = interpolate_motion(f0, f1, 0.5, motion_weight=1.0)
        half = interpolate_motion(f0, f1, 0.5, motion_weight=0.5)
        np.testing.assert_allclose(half, 0.5 * (lin + full))

    @pytest.mark.parametrize("bad", [-0.1, 1.5])
    def test_motion_weight_out_of_range_raises(self, bad):
        f0, _, f1 = _moving_pair(1, 1)
        with pytest.raises(ValueError, match="motion_weight"):
            interpolate_motion(f0, f1, 0.5, motion_weight=bad)

    def test_3d_requires_tracking_fields(self):
        f0, _, f1 = _moving_pair(1, 1)
        with pytest.raises(ValueError, match="track0 and track1"):
            interpolate_motion(np.stack([f0, f0]), np.stack([f1, f1]), 0.5)

    def test_3d_applies_one_motion_to_every_level(self):
        f0, fmid, f1 = _moving_pair(3.0, 6.0)
        scale = np.array([1.0, 0.5, 0.1])[:, None, None]
        out = interpolate_motion(
            scale * f0, scale * f1, 0.5, track0=f0, track1=f1, motion_weight=1.0
        )
        assert out.shape == (3, *f0.shape)
        for k in range(3):
            assert _rmse(out[k], scale[k, 0, 0] * fmid) < 0.05 * scale[k, 0, 0]

    def test_warp_zero_fraction_is_identity(self):
        f0 = _moving_pair(1, 1)[0]
        flow = np.ones((2, *f0.shape))
        np.testing.assert_allclose(warp(f0, flow, 0.0), f0)


def _tropo_ds(wet2d, time, lat, lon, heights=(0.0, 1500.0, 3000.0), hyd=2.3):
    scale = np.exp(-np.asarray(heights) / 2000.0)[:, None, None]
    wet = (scale * wet2d)[None].astype("float32")
    hydro = (hyd * np.exp(-np.asarray(heights) / 8000.0))[
        None, :, None, None
    ] * np.ones_like(wet)
    coords = {
        "time": [time],
        "height": list(heights),
        "latitude": lat,
        "longitude": lon,
    }
    dims = ("time", "height", "latitude", "longitude")
    return xr.Dataset(
        {
            "wet_delay": (dims, wet),
            "hydrostatic_delay": (dims, hydro.astype("float32")),
        },
        coords=coords,
    )


@pytest.fixture
def tropo_pair():
    f0, fmid, f1 = _moving_pair(3.0, 8.0, shape=(140, 200))
    lat = np.linspace(40.0, 40.0 - 0.07 * 139, 140)  # north to south, like the products
    lon = np.linspace(-100.0, -100.0 + 0.07 * 199, 200)
    t0, t1 = pd.Timestamp("2022-12-23T18:00"), pd.Timestamp("2022-12-24T00:00")

    def to_delay(z):
        return 0.10 + 0.03 * z  # metres, always positive

    return (
        _tropo_ds(to_delay(f0), t0, lat, lon),
        _tropo_ds(to_delay(f1), t1, lat, lon),
        to_delay(fmid),
        t0,
        t1,
    )


class TestInterpInTimeMotion:
    def test_output_matches_linear_layout(self, tropo_pair):
        ds0, ds1, _, t0, t1 = tropo_pair
        t = t0 + (t1 - t0) / 2
        lin = _interp_in_time(ds0, ds1, t0, t1, t)
        mot = interp_in_time_motion(ds0, ds1, t0, t1, t)
        assert (
            mot.total_delay.dims
            == lin.total_delay.dims
            == ("height", "latitude", "longitude")
        )
        assert mot.total_delay.shape == lin.total_delay.shape
        assert mot.total_delay.attrs["time_interpolation"] == "motion"
        assert mot.total_delay.attrs["flow_height"] == 1500.0

    def test_zero_weight_equals_linear(self, tropo_pair):
        ds0, ds1, _, t0, t1 = tropo_pair
        t = t0 + (t1 - t0) / 3
        lin = _interp_in_time(ds0, ds1, t0, t1, t)
        mot = interp_in_time_motion(ds0, ds1, t0, t1, t, motion_weight=0.0)
        np.testing.assert_allclose(mot.total_delay, lin.total_delay, rtol=1e-5)

    def test_closer_to_truth_than_linear(self, tropo_pair):
        ds0, ds1, wet_mid, t0, t1 = tropo_pair
        t = t0 + (t1 - t0) / 2
        hyd = float(ds0.hydrostatic_delay.isel(time=0, height=0).mean())
        truth = wet_mid + hyd
        lin = _interp_in_time(ds0, ds1, t0, t1, t).total_delay.isel(height=0).values
        mot = interp_in_time_motion(ds0, ds1, t0, t1, t, motion_weight=1.0)
        assert (
            _rmse(mot.total_delay.isel(height=0).values, truth) < _rmse(lin, truth) / 3
        )

    def test_hydrostatic_part_is_linear(self, tropo_pair):
        ds0, ds1, _, t0, t1 = tropo_pair
        ds1 = ds1.copy(deep=True)
        ds1["hydrostatic_delay"] = ds1.hydrostatic_delay + 0.012
        ds0["wet_delay"] = ds0.wet_delay * 0
        ds1["wet_delay"] = ds1.wet_delay * 0
        t = t0 + (t1 - t0) / 4
        out = interp_in_time_motion(ds0, ds1, t0, t1, t, motion_weight=1.0)
        expected = ds0.hydrostatic_delay.isel(time=0) + 0.25 * 0.012
        np.testing.assert_allclose(out.total_delay, expected, rtol=1e-5)


class TestCropTropoMotion:
    def _run(self, monkeypatch, tmp_path, tropo_pair, **kwargs):
        ds0, ds1, _, t0, t1 = tropo_pair
        calls = []

        def fake_open_crop(url, lat_bounds, lon_bounds, h_max):
            calls.append((lat_bounds, lon_bounds))
            ds = ds0 if "T180000" in url else ds1
            return ds.sel(
                latitude=slice(lat_bounds[0], lat_bounds[1]),
                longitude=slice(lon_bounds[0], lon_bounds[1]),
            )

        monkeypatch.setattr(_crop, "_open_crop", fake_open_crop)
        urls = pd.Series(
            ["tropo_20221223T180000Z.nc", "tropo_20221224T000000Z.nc"],
            index=pd.DatetimeIndex([t0, t1]),
        )
        status = _crop._process_one_datetime(
            datetime(2022, 12, 23, 21, 0, 0),
            urls,
            (36.0, 34.0),
            (-95.0, -92.0),
            10000.0,
            tmp_path,
            False,
            **kwargs,
        )
        out = xr.open_dataset(
            tmp_path / "tropo_cropped_20221223T210000.nc", engine="h5netcdf"
        )
        return status, calls, out

    def test_motion_reads_a_wider_area_and_trims_the_output(
        self, monkeypatch, tmp_path, tropo_pair
    ):
        status, calls, out = self._run(
            monkeypatch,
            tmp_path,
            tropo_pair,
            time_interpolation="motion",
            motion_margin_deg=3.0,
        )
        assert status[1] == "ok"
        assert calls[0] == ((39.0, 31.0), (-98.0, -89.0))
        assert float(out.latitude.max()) <= 36.0
        assert float(out.latitude.min()) >= 34.0
        assert float(out.longitude.min()) >= -95.0
        assert float(out.longitude.max()) <= -92.0
        assert np.isfinite(out.total_delay).all()

    def test_linear_is_unchanged_by_default(self, monkeypatch, tmp_path, tropo_pair):
        _, calls, out = self._run(monkeypatch, tmp_path, tropo_pair)
        assert calls[0] == ((36.0, 34.0), (-95.0, -92.0))
        ds0, ds1, _, t0, t1 = tropo_pair
        expected = _interp_in_time(
            ds0, ds1, t0, t1, t0 + (t1 - t0) / 2
        ).total_delay.sel(latitude=slice(36.0, 34.0), longitude=slice(-95.0, -92.0))
        np.testing.assert_allclose(out.total_delay, expected, rtol=1e-6)

    def test_motion_differs_from_linear(self, monkeypatch, tmp_path, tropo_pair):
        motion_dir, linear_dir = tmp_path / "motion", tmp_path / "linear"
        motion_dir.mkdir()
        linear_dir.mkdir()
        _, _, mot = self._run(
            monkeypatch, motion_dir, tropo_pair, time_interpolation="motion"
        )
        _, _, lin = self._run(monkeypatch, linear_dir, tropo_pair)
        assert float(np.abs(mot.total_delay - lin.total_delay).max()) > 1e-3

    @pytest.mark.parametrize(
        ("kwargs", "match"),
        [
            ({"time_interpolation": "cubic"}, "time_interpolation"),
            (
                {"time_interpolation": "motion", "skip_time_interpolation": True},
                "skip_time",
            ),
            ({"time_interpolation": "motion", "motion_weight": 2.0}, "motion_weight"),
        ],
    )
    def test_crop_tropo_rejects_bad_options(self, tmp_path, kwargs, match):
        urls = tmp_path / "urls.txt"
        urls.write_text("tropo_20221223T180000Z.nc\n")
        with pytest.raises(ValueError, match=match):
            _crop.crop_tropo(
                urls, [], aoi_bounds=(-95, 34, -92, 36), output_dir=tmp_path, **kwargs
            )
