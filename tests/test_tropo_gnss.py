"""Tests for opera_utils.tropo._gnss and its use in apply_tropo."""

from __future__ import annotations

import gzip
from datetime import datetime

import numpy as np
import pandas as pd
import pytest
import xarray as xr

from opera_utils.tropo import _apply
from opera_utils.tropo._gnss import (
    gnss_residual_field,
    krige,
    load_gnss_ztd,
    model_ztd_at_stations,
    read_sinex_tropo,
    station_residual_difference,
    tune_kriging,
)

T_REF = datetime(2020, 10, 2, 13, 52, 36)
T_SEC = datetime(2020, 10, 26, 13, 52, 36)


def _stations(n=60, seed=1):
    rng = np.random.default_rng(seed)
    return pd.DataFrame(
        {
            "id": [f"S{i:03d}" for i in range(n)],
            "lat": rng.uniform(33.2, 34.6, n),
            "lon": rng.uniform(-118.8, -116.6, n),
            "height": rng.uniform(0.0, 1500.0, n),
        }
    )


def _smooth_error(lat, lon):
    """A model error field that varies over ~50 km, in meters of zenith delay."""
    return 0.02 * np.sin((lat - 33.0) * 4.0) * np.cos((lon + 118.0) * 3.0)


def _cube(scale=1.0):
    h = np.array([0.0, 500.0, 1000.0, 2000.0, 4000.0])
    lat = np.linspace(35.0, 33.0, 41)
    lon = np.linspace(-119.0, -116.0, 61)
    profile = 2.4 * np.exp(-h / 7000.0) * scale
    data = profile[:, None, None] * np.ones((1, lat.size, lon.size))
    return xr.DataArray(
        data,
        dims=("height", "latitude", "longitude"),
        coords={"height": h, "latitude": lat, "longitude": lon},
        name="total_delay",
    )


def _gnss_table(st, err_sec, bias=None):
    """GNSS = model + per-station bias (+ model error on the second date)."""
    bias = np.zeros(len(st)) if bias is None else bias
    ref = model_ztd_at_stations(_cube(1.0), st.lat, st.lon, st.height) + bias
    sec = model_ztd_at_stations(_cube(1.01), st.lat, st.lon, st.height) + bias + err_sec
    rows = [
        pd.DataFrame({**st.to_dict("list"), "datetime": when, "ztd": z})
        for when, z in ((T_REF, ref), (T_SEC, sec))
    ]
    return pd.concat(rows, ignore_index=True)


class TestKriging:
    def test_exact_at_stations_without_nugget(self):
        rng = np.random.default_rng(0)
        xy = rng.uniform(0, 100, (30, 2))
        val = rng.standard_normal(30)
        np.testing.assert_allclose(krige(xy, val, xy, 20.0, 1e-10), val, atol=1e-6)

    def test_decays_to_zero_far_from_stations(self):
        xy = np.array([[0.0, 0.0], [5.0, 0.0], [0.0, 5.0]])
        far = np.array([[2000.0, 2000.0]])
        assert abs(krige(xy, np.array([3.0, 2.0, 4.0]), far, 20.0, 0.05)[0]) < 1e-6

    def test_nugget_smooths(self):
        xy = np.array([[0.0, 0.0], [1.0, 0.0]])
        val = np.array([1.0, -1.0])
        rough = krige(xy, val, xy, 50.0, 1e-8)
        smooth = krige(xy, val, xy, 50.0, 1.0)
        assert np.abs(smooth).max() < np.abs(rough).max()

    def test_tuning_recovers_a_smooth_field(self):
        st = _stations(120)
        xy = np.column_stack([st.lat * 111.0, st.lon * 92.0])
        val = _smooth_error(st.lat.values, st.lon.values)
        rng_km, _nugget, loo = tune_kriging(xy, val)
        assert rng_km >= 25.0
        assert loo < 0.3 * val.std()


class TestModelAtStations:
    def test_follows_the_vertical_profile(self):
        z = model_ztd_at_stations(
            _cube(), [34.0, 34.0], [-117.5, -117.5], [0.0, 1000.0]
        )
        np.testing.assert_allclose(
            z, 2.4 * np.exp(-np.array([0.0, 1000.0]) / 7000.0), rtol=2e-3
        )

    def test_outside_the_cube_is_nan_and_below_is_clamped(self):
        z = model_ztd_at_stations(_cube(), [40.0, 34.0], [-117.5, -117.5], [0.0, -50.0])
        assert np.isnan(z[0])
        assert z[1] == pytest.approx(2.4, rel=1e-6)

    def test_latitude_order_does_not_matter(self):
        cube = _cube()
        flipped = cube.isel(latitude=slice(None, None, -1))
        a = model_ztd_at_stations(cube, [33.7], [-117.2], [300.0])
        b = model_ztd_at_stations(flipped, [33.7], [-117.2], [300.0])
        np.testing.assert_allclose(a, b)


class TestStationResiduals:
    def test_constant_station_biases_cancel(self):
        st = _stations()
        err = _smooth_error(st.lat.values, st.lon.values)
        bias = np.random.default_rng(3).normal(0, 0.03, len(st))
        res = station_residual_difference(
            _gnss_table(st, err, bias), _cube(1.01), T_SEC, _cube(1.0), T_REF
        )
        np.testing.assert_allclose(
            res.residual.values, err[[int(i[1:]) for i in res.index]], atol=1e-9
        )

    def test_only_stations_on_both_dates_are_used(self):
        st = _stations(20)
        table = _gnss_table(st, np.zeros(20))
        table = table[~((table.id == "S003") & (table.datetime == pd.Timestamp(T_REF)))]
        res = station_residual_difference(table, _cube(1.01), T_SEC, _cube(1.0), T_REF)
        assert "S003" not in res.index
        assert len(res) == 19

    def test_outliers_are_dropped(self):
        st = _stations(40)
        err = _smooth_error(st.lat.values, st.lon.values)
        err[7] += 1.0
        res = station_residual_difference(
            _gnss_table(st, err), _cube(1.01), T_SEC, _cube(1.0), T_REF
        )
        assert "S007" not in res.index

    def test_time_tolerance(self):
        st = _stations(10)
        table = _gnss_table(st, np.zeros(10))
        table.loc[table.datetime == pd.Timestamp(T_SEC), "datetime"] += pd.Timedelta(
            minutes=30
        )
        res = station_residual_difference(table, _cube(1.01), T_SEC, _cube(1.0), T_REF)
        assert res.empty


class TestResidualField:
    def test_recovers_a_smooth_model_error(self):
        st = _stations(150)
        err = _smooth_error(st.lat.values, st.lon.values)
        res = station_residual_difference(
            _gnss_table(st, err), _cube(1.01), T_SEC, _cube(1.0), T_REF
        )
        lat, lon = np.meshgrid(
            np.linspace(33.4, 34.4, 25), np.linspace(-118.5, -116.9, 30), indexing="ij"
        )
        field, info = gnss_residual_field(res, lat, lon)
        truth = _smooth_error(lat, lon) - np.median(res.residual)
        assert np.sqrt(np.mean((field - truth) ** 2)) < 0.25 * truth.std()
        assert info["n_stations"] == len(res)
        assert info["loo_rms"] < info["residual_rms"]

    def test_too_few_stations_gives_zero(self):
        st = _stations(3)
        res = station_residual_difference(
            _gnss_table(st, np.full(3, 0.05)), _cube(1.01), T_SEC, _cube(1.0), T_REF
        )
        field, info = gnss_residual_field(res, np.array([34.0]), np.array([-117.5]))
        assert field[0] == 0.0
        assert info["n_stations"] == 3

    def test_a_constant_offset_is_not_spread(self):
        st = _stations(40)
        res = station_residual_difference(
            _gnss_table(st, np.full(40, 0.08)), _cube(1.01), T_SEC, _cube(1.0), T_REF
        )
        field, _ = gnss_residual_field(
            res, np.array([34.0, 45.0]), np.array([-117.5, -100.0])
        )
        np.testing.assert_allclose(field, 0.0, atol=1e-9)


class TestSinexTropo:
    TEXT = """%=TRO 2.00 NGL 20:300:00000
+TROP/SOLUTION
*SITE ____EPOCH___ TROTOT STDDEV  TGNTOT STDDEV  TGETOT STDDEV
 OKCB 20:300:49800 2412.3    1.2   -0.31   0.20    0.11   0.20
 OKCB 20:300:50100 2413.1    1.3   -0.30   0.20    0.12   0.20
 OKCB 20:300:50400 bad       1.3   -0.30   0.20    0.12   0.20
-TROP/SOLUTION
%=ENDTRO
"""

    def test_parses_plain_and_gzipped(self):
        for raw in (self.TEXT.encode(), gzip.compress(self.TEXT.encode())):
            df = read_sinex_tropo(raw)
            assert len(df) == 2
            assert df.ztd.iloc[0] == pytest.approx(2.4123)
            assert df.ztd_sigma.iloc[1] == pytest.approx(0.0013)
            assert df.datetime.iloc[0] == pd.Timestamp("2020-10-26 13:50:00")

    def test_empty_file(self):
        assert read_sinex_tropo(b"nothing here").empty


class TestTableAndApply:
    def test_load_checks_columns(self, tmp_path):
        good = _gnss_table(_stations(5), np.zeros(5))
        good.to_csv(tmp_path / "g.csv", index=False)
        assert len(load_gnss_ztd(tmp_path / "g.csv")) == 10
        good.drop(columns="height").to_csv(tmp_path / "bad.csv", index=False)
        with pytest.raises(ValueError, match="missing columns"):
            load_gnss_ztd(tmp_path / "bad.csv")

    def test_apply_tropo_requires_a_reference_date(self, tmp_path):
        with pytest.raises(ValueError, match="subtract_first_date"):
            _apply.apply_tropo(
                [tmp_path / "tropo_cropped_20201002T135236.nc"],
                tmp_path / "dem.tif",
                tmp_path / "inc.tif",
                subtract_first_date=False,
                gnss_ztd_file=tmp_path / "g.csv",
            )

    def test_zenith_field_on_a_projected_dem(self):
        rioxarray = pytest.importorskip("rioxarray")  # noqa: F841
        st = _stations(120)
        err = _smooth_error(st.lat.values, st.lon.values)
        x = np.arange(380_000.0, 520_000.0, 2000.0)
        y = np.arange(3_820_000.0, 3_700_000.0, -2000.0)
        dem = xr.DataArray(
            np.full((y.size, x.size), 200.0), dims=("y", "x"), coords={"y": y, "x": x}
        ).rio.write_crs("epsg:32611")
        field, info = _apply._gnss_zenith_field(
            _gnss_table(st, err),
            _cube(1.01).to_dataset(),
            _cube(1.0).to_dataset(),
            pd.Timestamp(T_SEC),
            pd.Timestamp(T_REF),
            dem,
        )
        assert field.shape == dem.shape
        assert field.dtype == np.float32
        assert info["n_stations"] > 100
        assert 0.3 * err.std() < field.std() < 1.5 * err.std()
