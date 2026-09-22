"""Tests for opera_utils.tropo._guide (ERA5 + NEXRAD guided time interpolation)."""

from __future__ import annotations

import sys
from datetime import datetime

import numpy as np
import pandas as pd
import pytest
import xarray as xr
from affine import Affine

from opera_utils.tropo import _crop, _guide
from opera_utils.tropo._guide import (
    DEFAULT_ERA5_ONLY_WEIGHT,
    DEFAULT_GUIDE_WEIGHTS,
    _echo_fraction_on_grid,
    guide_radar_times,
    interp_in_time_guided,
    load_era5_tcwv,
    read_nexrad_echo_fraction,
)
from opera_utils.tropo._helpers import _interp_in_time

T0, T1 = pd.Timestamp("2021-08-29T18:00"), pd.Timestamp("2021-08-30T00:00")
T = pd.Timestamp("2021-08-29T21:00")
LAT = np.round(np.linspace(40.0, 40.0 - 0.07 * 79, 80), 6)  # north to south
LON = np.round(np.linspace(-100.0, -100.0 + 0.07 * 119, 120), 6)
HEIGHTS = (0.0, 1500.0, 3000.0)


def _tropo_ds(wet2d, time, heights=HEIGHTS, hyd=2.3):
    scale = np.exp(-np.asarray(heights) / 2000.0)[:, None, None]
    wet = (scale * wet2d)[None].astype("float32")
    hydro = (hyd * np.exp(-np.asarray(heights) / 8000.0))[
        None, :, None, None
    ] * np.ones_like(wet)
    dims = ("time", "height", "latitude", "longitude")
    return xr.Dataset(
        {
            "wet_delay": (dims, wet),
            "hydrostatic_delay": (dims, hydro.astype("float32")),
        },
        coords={
            "time": [time],
            "height": list(heights),
            "latitude": LAT,
            "longitude": LON,
        },
    )


def _blob(col, row=40.0, width=6.0):
    yy, xx = np.mgrid[: LAT.size, : LON.size]
    return np.exp(-(((yy - row) ** 2 + (xx - col) ** 2) / (2 * width**2)))


def _era5(values_at):
    """ERA5-like TCWV on a coarse grid, hourly from T0 to T1: value = f(time)."""
    times = pd.date_range(
        T0 - pd.Timedelta(hours=1), T1 + pd.Timedelta(hours=1), freq="1h"
    )
    lat = np.arange(42.0, 33.0, -0.25)
    lon = np.arange(-102.0, -89.0, 0.25)
    data = np.stack([np.full((lat.size, lon.size), values_at(t)) for t in times])
    return xr.DataArray(
        data,
        coords={"time": times, "latitude": lat, "longitude": lon},
        dims=("time", "latitude", "longitude"),
    )


def _straight(t):
    return 30.0 + 10.0 * (t - T0) / (T1 - T0)


@pytest.fixture
def pair():
    wet0 = 0.10 + 0.05 * _blob(30.0)
    wet1 = 0.10 + 0.05 * _blob(60.0)
    return _tropo_ds(wet0, T0), _tropo_ds(wet1, T1)


def _radar(value=0.0, *, missing=None, override=None):
    out = {}
    for x in guide_radar_times(
        T0.to_pydatetime(), T1.to_pydatetime(), T.to_pydatetime()
    ):
        out[x] = np.full((LAT.size, LON.size), value)
    if override:
        out.update(override)
    if missing is not None:
        out[missing] = None
    return out


class TestRadarTimes:
    def test_hourly_plus_the_acquisition(self):
        times = guide_radar_times(
            datetime(2020, 10, 26, 12),
            datetime(2020, 10, 26, 18),
            datetime(2020, 10, 26, 13, 52, 36),
        )
        assert times[0] == datetime(2020, 10, 26, 12)
        assert times[-1] == datetime(2020, 10, 26, 18)
        assert datetime(2020, 10, 26, 13, 55) in times
        assert len(times) == 8

    def test_on_the_hour_is_not_duplicated(self):
        times = guide_radar_times(
            T0.to_pydatetime(), T1.to_pydatetime(), T.to_pydatetime()
        )
        assert len(times) == 7
        assert len(set(times)) == 7


class TestInterpInTimeGuided:
    def test_era5_on_its_straight_line_is_linear(self, pair):
        ds0, ds1 = pair
        out = interp_in_time_guided(ds0, ds1, T0, T1, T, era5_tcwv=_era5(_straight))
        lin = _interp_in_time(ds0, ds1, T0, T1, T)
        np.testing.assert_allclose(out.total_delay, lin.total_delay, rtol=1e-6)
        assert out.total_delay.attrs["time_interpolation"] == "guided"
        assert out.total_delay.attrs["guide_sources"] == "ERA5"
        assert out.total_delay.dims == ("height", "latitude", "longitude")

    def test_era5_departure_scales_the_wet_delay(self, pair):
        ds0, ds1 = pair
        era5 = _era5(lambda t: _straight(t) * (1.10 if t == T else 1.0))
        out = interp_in_time_guided(ds0, ds1, T0, T1, T, era5_tcwv=era5)
        lin = _interp_in_time(ds0, ds1, T0, T1, T).total_delay
        hyd = 0.5 * (
            ds0.hydrostatic_delay.isel(time=0, drop=True)
            + ds1.hydrostatic_delay.isel(time=0, drop=True)
        )
        wet_lin = lin - hyd
        expected = hyd + wet_lin * (1 + DEFAULT_ERA5_ONLY_WEIGHT * 0.10)
        np.testing.assert_allclose(out.total_delay, expected, rtol=1e-5)

    def test_hydrostatic_part_is_linear(self, pair):
        ds0, ds1 = pair
        ds1 = ds1.copy(deep=True)
        ds1["hydrostatic_delay"] = ds1.hydrostatic_delay + 0.012
        ds0 = ds0.copy(deep=True)
        ds0["wet_delay"] = ds0.wet_delay * 0
        ds1["wet_delay"] = ds1.wet_delay * 0
        era5 = _era5(lambda t: _straight(t) * (1.3 if t == T else 1.0))
        out = interp_in_time_guided(
            ds0, ds1, T0, T1, T, era5_tcwv=era5, radar=_radar(0.2)
        )
        expected = ds0.hydrostatic_delay.isel(time=0) + 0.5 * 0.012
        np.testing.assert_allclose(out.total_delay, expected, rtol=1e-5)

    def test_quiet_radar_uses_the_combined_weights(self, pair):
        ds0, ds1 = pair
        era5 = _era5(lambda t: _straight(t) * (1.10 if t == T else 1.0))
        out = interp_in_time_guided(
            ds0, ds1, T0, T1, T, era5_tcwv=era5, radar=_radar(0.0)
        )
        lin = _interp_in_time(ds0, ds1, T0, T1, T).total_delay
        hyd = 0.5 * (
            ds0.hydrostatic_delay.isel(time=0, drop=True)
            + ds1.hydrostatic_delay.isel(time=0, drop=True)
        )
        expected = hyd + (lin - hyd) * (1 + DEFAULT_GUIDE_WEIGHTS["era5"] * 0.10)
        np.testing.assert_allclose(out.total_delay, expected, rtol=1e-5)
        assert out.total_delay.attrs["guide_sources"] == "ERA5 + NEXRAD"

    def test_missing_radar_frame_falls_back_to_era5(self, pair):
        ds0, ds1 = pair
        out = interp_in_time_guided(
            ds0,
            ds1,
            T0,
            T1,
            T,
            era5_tcwv=_era5(_straight),
            radar=_radar(0.3, missing=datetime(2021, 8, 29, 20)),
        )
        assert out.total_delay.attrs["guide_sources"] == "ERA5"
        assert out.total_delay.attrs["guide_weight_era5"] == DEFAULT_ERA5_ONLY_WEIGHT

    def test_echo_cover_change_adds_wet_delay_along_the_profile(self, pair):
        ds0, ds1 = pair
        radar = _radar(
            0.0, override={T.to_pydatetime(): np.full((LAT.size, LON.size), 0.5)}
        )
        out = interp_in_time_guided(
            ds0, ds1, T0, T1, T, era5_tcwv=_era5(_straight), radar=radar
        )
        lin = _interp_in_time(ds0, ds1, T0, T1, T).total_delay
        added = (out.total_delay - lin).values
        np.testing.assert_allclose(
            added[0], DEFAULT_GUIDE_WEIGHTS["cover"] * 0.5, rtol=1e-4
        )
        # higher levels follow the wet delay profile, exp(-h / 2 km) in the fixture
        np.testing.assert_allclose(added[1] / added[0], np.exp(-0.75), rtol=1e-4)

    def test_radar_motion_moves_the_field(self, pair):
        ds0, ds1 = pair
        hours = guide_radar_times(
            T0.to_pydatetime(), T1.to_pydatetime(), T0.to_pydatetime()
        )
        # radar echoes travel with the moist blob: 30 px over 6 h, 5 px per hour
        radar = {
            x: np.clip(_blob(30.0 + 5.0 * i, width=5.0) * 1.5, 0, 1)
            for i, x in enumerate(hours)
        }
        radar[T.to_pydatetime()] = radar[hours[3]]
        motion_only = {"era5": 0.0, "cover": 0.0, "motion": 1.0}
        out = interp_in_time_guided(
            ds0,
            ds1,
            T0,
            T1,
            T,
            era5_tcwv=_era5(_straight),
            radar=radar,
            weights=motion_only,
        )
        truth = 0.10 + 0.05 * _blob(45.0)
        wet = out.total_delay.isel(height=0).values - 2.3
        lin = (
            _interp_in_time(ds0, ds1, T0, T1, T).total_delay.isel(height=0).values - 2.3
        )
        err = lambda a: float(np.sqrt(np.mean((a - truth) ** 2)))  # noqa: E731
        assert err(wet) < 0.5 * err(lin)

    def test_era5_must_cover_the_grid(self, pair):
        ds0, ds1 = pair
        era5 = _era5(_straight).assign_coords(longitude=lambda d: d.longitude + 40.0)
        with pytest.raises(ValueError, match="does not cover the TROPO grid"):
            interp_in_time_guided(ds0, ds1, T0, T1, T, era5_tcwv=era5)

    def test_era5_north_to_south_or_south_to_north_gives_the_same(self, pair):
        ds0, ds1 = pair
        era5 = _era5(lambda t: _straight(t) * (1.10 if t == T else 1.0))
        a = interp_in_time_guided(ds0, ds1, T0, T1, T, era5_tcwv=era5)
        b = interp_in_time_guided(
            ds0, ds1, T0, T1, T, era5_tcwv=era5.sortby("latitude")
        )
        np.testing.assert_allclose(a.total_delay, b.total_delay)

    def test_era5_must_cover_the_interval(self, pair):
        ds0, ds1 = pair
        era5 = _era5(_straight).sel(time=slice(None, T))
        with pytest.raises(ValueError, match="ERA5 does not cover"):
            interp_in_time_guided(ds0, ds1, T0, T1, T, era5_tcwv=era5)


class TestNexradReading:
    def test_fraction_on_grid(self):
        # mosaic at 0.01 deg starting at (-100.035, 40.035): cells of 7 x 7 pixels
        idx = np.zeros((300, 400), dtype=np.uint8)
        dbz30 = int((30.0 + 32.5) / 0.5)
        idx[:70, :70] = dbz30  # the first 10 x 10 cells, fully echo
        idx[:7, 70:77] = int((10.0 + 32.5) / 0.5)  # weak echo: below 20 dBZ
        transform = Affine(0.01, 0, -100.035, 0, -0.01, 40.035)
        frac = _echo_fraction_on_grid(idx, transform, LAT, LON)
        assert frac.shape == (LAT.size, LON.size)
        np.testing.assert_allclose(frac[:10, :10], 1.0)
        assert frac[0, 10] == 0.0
        # the mosaic covers only 3 x 4 degrees: cells beyond are NaN
        assert np.isnan(frac[-1, -1])

    def test_ascending_latitude_grid_is_handled(self):
        idx = np.zeros((100, 100), dtype=np.uint8)
        idx[:20, :] = 150  # northern rows
        transform = Affine(0.01, 0, -100.035, 0, -0.01, 40.035)
        lat_up = LAT[::-1]
        frac = _echo_fraction_on_grid(idx, transform, lat_up, LON)
        assert frac[-1, 0] > 0.9
        assert frac[0, 0] != frac[-1, 0]

    @pytest.mark.parametrize("driver", ["GTiff", "PNG"])
    def test_read_from_file_with_georeferencing(self, tmp_path, driver):
        rasterio = pytest.importorskip("rasterio")
        idx = np.zeros((300, 400), dtype=np.uint8)
        idx[:70, :70] = 150
        transform = Affine(0.01, 0, -100.035, 0, -0.01, 40.035)
        path = tmp_path / ("n0q.tif" if driver == "GTiff" else "n0q.png")
        opts = {"WORLDFILE": "YES"} if driver == "PNG" else {}
        with rasterio.open(
            path,
            "w",
            driver=driver,
            width=400,
            height=300,
            count=1,
            dtype="uint8",
            crs="EPSG:4326",
            transform=transform,
            **opts,
        ) as dst:
            dst.write(idx, 1)
        frac = read_nexrad_echo_fraction(path, LAT, LON)
        np.testing.assert_allclose(frac[:10, :10], 1.0)
        assert frac[20, 20] == 0.0


class TestDownloadsAndLoading:
    def test_nexrad_download_requests_the_needed_times(self, tmp_path, monkeypatch):
        fetched = []

        def fake_retrieve(url, dst):
            fetched.append(url)
            with open(dst, "wb") as f:
                f.write(b"x")

        monkeypatch.setattr(_guide.urllib.request, "urlretrieve", fake_retrieve)
        paths = _guide.download_nexrad_n0q(
            [datetime(2020, 10, 26, 13, 52, 36)], tmp_path
        )
        assert len(paths) == 8
        assert any("n0q_202010261355.png" in u for u in fetched)
        assert any("n0q_202010261200.wld" in u for u in fetched)
        # a second call downloads nothing new
        fetched.clear()
        _guide.download_nexrad_n0q([datetime(2020, 10, 26, 13, 52, 36)], tmp_path)
        assert fetched == []

    def test_era5_download_needs_cdsapi(self, tmp_path, monkeypatch):
        monkeypatch.setitem(sys.modules, "cdsapi", None)
        with pytest.raises(ImportError, match="cdsapi"):
            _guide.download_era5_tcwv(
                (-100, 30, -90, 40), [datetime(2020, 1, 1, 3)], tmp_path
            )

    def test_load_era5_from_a_directory(self, tmp_path):
        era5 = _era5(_straight)
        for k, part in enumerate(
            (era5.isel(time=slice(0, 4)), era5.isel(time=slice(3, None)))
        ):
            ds = part.rename(time="valid_time").to_dataset(name="tcwv")
            ds = ds.assign_coords(number=0).expand_dims(number=[0])
            ds.to_netcdf(tmp_path / f"era5_tcwv_{k}.nc")
        out = load_era5_tcwv(tmp_path)
        assert out.dims == ("time", "latitude", "longitude")
        assert out.time.size == era5.time.size  # the overlapping hour is kept once
        assert (np.diff(out.time.values) > np.timedelta64(0)).all()

    def test_load_era5_empty_directory_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            load_era5_tcwv(tmp_path)


class TestCropTropoGuided:
    def _run(self, monkeypatch, tmp_path, pair, **kwargs):
        ds0, ds1 = pair
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
            ["tropo_20210829T180000Z.nc", "tropo_20210830T000000Z.nc"],
            index=pd.DatetimeIndex([T0, T1]),
        )
        status = _crop._process_one_datetime(
            T.to_pydatetime(),
            urls,
            (38.0, 36.0),
            (-96.0, -94.0),
            10000.0,
            tmp_path,
            False,
            **kwargs,
        )
        out = xr.open_dataset(
            tmp_path / "tropo_cropped_20210829T210000.nc", engine="h5netcdf"
        )
        return status, calls, out

    def test_guided_reads_a_wider_area_trims_and_differs_from_linear(
        self, monkeypatch, tmp_path, pair
    ):
        era5 = _era5(lambda t: _straight(t) * (1.2 if t == T else 1.0))
        f = tmp_path / "era5.nc"
        era5.to_dataset(name="tcwv").to_netcdf(f)
        (tmp_path / "g").mkdir()
        (tmp_path / "l").mkdir()
        status, calls, out = self._run(
            monkeypatch,
            tmp_path / "g",
            pair,
            time_interpolation="guided",
            era5_file=f,
            motion_margin_deg=1.0,
        )
        assert status[1] == "ok"
        assert calls[0] == ((39.0, 35.0), (-97.0, -93.0))
        assert float(out.latitude.max()) <= 38.0
        assert float(out.latitude.min()) >= 36.0
        assert out.total_delay.attrs["time_interpolation"] == "guided"
        _, _, lin = self._run(monkeypatch, tmp_path / "l", pair)
        assert float((out.total_delay - lin.total_delay).max()) > 1e-3

    def test_guided_uses_nexrad_when_present(self, monkeypatch, tmp_path, pair):
        rasterio = pytest.importorskip("rasterio")
        era5 = _era5(_straight)
        f = tmp_path / "era5.nc"
        era5.to_dataset(name="tcwv").to_netcdf(f)
        radar_dir = tmp_path / "radar"
        radar_dir.mkdir()
        transform = Affine(0.01, 0, -101.0, 0, -0.01, 41.0)
        for x in guide_radar_times(
            T0.to_pydatetime(), T1.to_pydatetime(), T.to_pydatetime()
        ):
            idx = np.zeros((800, 1000), dtype=np.uint8)
            with rasterio.open(
                _guide.nexrad_path(radar_dir, x),
                "w",
                driver="PNG",
                width=1000,
                height=800,
                count=1,
                dtype="uint8",
                crs="EPSG:4326",
                transform=transform,
                WORLDFILE="YES",
            ) as dst:
                dst.write(idx, 1)
        (tmp_path / "o").mkdir()
        _, _, out = self._run(
            monkeypatch,
            tmp_path / "o",
            pair,
            time_interpolation="guided",
            era5_file=f,
            nexrad_dir=radar_dir,
            motion_margin_deg=1.0,
        )
        assert out.total_delay.attrs["guide_sources"] == "ERA5 + NEXRAD"

    @pytest.mark.parametrize(
        ("kwargs", "match"),
        [
            ({"time_interpolation": "guided"}, "era5_file"),
            (
                {
                    "time_interpolation": "guided",
                    "era5_file": "x.nc",
                    "skip_time_interpolation": True,
                },
                "skip_time",
            ),
            ({"time_interpolation": "cubic"}, "time_interpolation"),
        ],
    )
    def test_crop_tropo_rejects_bad_options(self, tmp_path, kwargs, match):
        urls = tmp_path / "urls.txt"
        urls.write_text("tropo_20210829T180000Z.nc\n")
        with pytest.raises(ValueError, match=match):
            _crop.crop_tropo(
                urls, [], aoi_bounds=(-96, 36, -94, 38), output_dir=tmp_path, **kwargs
            )


def test_smoothing_of_a_uniform_field_is_uniform():
    f = np.full((30, 40), 0.3)
    np.testing.assert_allclose(_guide._nan_smooth(f, 4.0), 0.3)
