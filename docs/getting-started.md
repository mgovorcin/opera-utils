# Getting started

## Install

`opera-utils` is available on conda-forge and PyPI:

=== "conda / mamba"

    ```bash
    mamba install -c conda-forge opera-utils
    ```

=== "pip"

    ```bash
    pip install opera-utils
    ```

The base install is intentionally light. Format-specific functionality is opt-in through extras:

| Extra | Unlocks | Pulls in |
| --- | --- | --- |
| `opera-utils[geopandas]` | Frame/burst *geometries* as GeoDataFrames, `disp-s1-intersects` CLI | `geopandas`, `pyogrio` |
| `opera-utils[asf]` | Searching/downloading CSLC, RTC, and static-layer products from ASF | `asf_search` |
| `opera-utils[disp]` | Searching, downloading, and reformatting DISP-S1 products | `xarray`, `dask`, `rasterio`, `rioxarray`, `zarr`, `h5netcdf`, `fsspec`, `s3fs`, `botocore`, `tqdm` |
| `opera-utils[nisar]` | Searching and downloading NISAR GSLC products | `fsspec`, `s3fs`, `aiohttp`, `tqdm` |
| `opera-utils[tropo]` | Cropping/applying OPERA TROPO tropospheric delay corrections | everything in `disp`, plus `scipy` |
| `opera-utils[all]` | Everything above | — |

```bash
pip install "opera-utils[disp]"     # just DISP-S1 support
pip install "opera-utils[all]"      # everything
```

!!! tip "Which extra do I need?"
    If you're not sure yet, start with `opera-utils[all]` — the extras mostly add dependencies for format-specific I/O (xarray, rasterio, geopandas), not core logic, so installing them all has little downside outside of image size.

## Set up data access credentials

Downloading real products (as opposed to just looking up frame/burst metadata) requires a free [NASA Earthdata Login](https://urs.earthdata.nasa.gov/users/new). The simplest setup is a `~/.netrc` entry:

```text title="~/.netrc"
machine urs.earthdata.nasa.gov
    login <your_username>
    password <your_password>
```

`opera-utils` also accepts `EARTHDATA_USERNAME`/`EARTHDATA_PASSWORD` environment variables for most (but not all) code paths — see [Set up Earthdata and AWS credentials](how-to-guides.md#set-up-earthdata-and-aws-credentials) for the full picture, including direct S3 access.

## 5-minute quickstart

These first three examples only need the base install — no extras, no credentials.

**Parse a Sentinel-1 burst ID out of a product filename:**

```pycon
>>> import opera_utils
>>> opera_utils.get_burst_id(
...     "OPERA_L2_CSLC-S1_T087-185683-IW2_20230322T161649Z_20240504T185235Z_S1A_VV_v1.1.h5"
... )
't087_185683_iw2'
```

**Look up the bursts and bounding box that make up a DISP-S1 frame** (the first call downloads and caches a small database file from GitHub):

```pycon
>>> import opera_utils
>>> opera_utils.get_burst_ids_for_frame(11114)
['t042_088905_iw1', 't042_088905_iw2', ..., 't042_088913_iw3']
>>> opera_utils.get_frame_bbox(11114)
(32610, Bbox(left=546450.0, bottom=4204110.0, right=833790.0, top=4409070.0))
```

**Same lookup from the command line:**

```bash
opera-utils disp-s1-frame-bbox 11114
# {"epsg": 32610, "bbox": [546450.0, 4204110.0, 833790.0, 4409070.0]}
```

Run `opera-utils --help` to see every subcommand available in your current install (the list grows as you install extras).

## Next steps

- Walk through a full workflow in [Tutorials](tutorials.md) — from finding a frame to plotting a displacement time series or forming a NISAR interferogram.
- Jump straight to a specific task in [How-to guides](how-to-guides.md).
- Read [Background theory](background-theory.md) to understand bursts, frames, and why DISP-S1 stacks need to be "rebased" before analysis.
