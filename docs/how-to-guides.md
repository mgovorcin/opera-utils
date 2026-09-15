# How-to guides

Task-oriented recipes. Each one assumes you've already done [Getting started](getting-started.md) (package installed, credentials set up). Extras required beyond the base install are called out per-section.

## Setup

### Install the extras you need

```bash
pip install "opera-utils[disp]"              # one extra
pip install "opera-utils[disp,nisar,tropo]"  # several
pip install "opera-utils[all]"               # everything
```

See the [extras table](getting-started.md#install) for what each one unlocks.

### Set up Earthdata and AWS credentials

Most download/streaming functions resolve Earthdata Login credentials in this order:

1. Credentials passed directly as function arguments (where supported, e.g. `open_h5(url, earthdata_username=..., earthdata_password=...)`).
2. A `~/.netrc` entry for `urs.earthdata.nasa.gov` (or `uat.urs.earthdata.nasa.gov` for the UAT/test environment).
3. `EARTHDATA_USERNAME` / `EARTHDATA_PASSWORD` environment variables.

```text title="~/.netrc"
machine urs.earthdata.nasa.gov
    login <your_username>
    password <your_password>
```

!!! warning "`opera_utils.download` is stricter"
    `opera_utils.download`'s `download_cslcs`/`download_rtcs`/etc. (ASF-backed CSLC/RTC downloads) read `~/.netrc` directly and do **not** fall back to the environment variables. Set up `.netrc` if you plan to use these.

For direct S3 access (`url_type="s3"`, only usable from an AWS `us-west-2` EC2 instance for Earthdata-hosted buckets), credentials are resolved from `AWS_ACCESS_KEY_ID`/`AWS_SECRET_ACCESS_KEY`/`AWS_SESSION_TOKEN` env vars, or automatically requested as temporary credentials from ASF using your Earthdata Login:

```python
from opera_utils.credentials import AWSCredentials, ASFCredentialEndpoints
import os

creds = AWSCredentials.from_asf(ASFCredentialEndpoints.OPERA)
os.environ.update(creds.to_env())
```

or print them for shell-sourcing:

```bash
python -m opera_utils.credentials
```

Temporary AWS credentials are cached for 50 minutes internally, so repeated calls don't re-hit the ASF endpoint.

## Bursts and frames

*Base install only.*

### Parse a burst ID out of a filename

```python
import opera_utils

opera_utils.get_burst_id("OPERA_L2_CSLC-S1_T087-185683-IW2_20230322T161649Z_20240504T185235Z_S1A_VV_v1.1.h5")
# 't087_185683_iw2'

opera_utils.normalize_burst_id("T087-185683-IW3")
# 't087_185683_iw3'
```

`get_burst_id` raises `ValueError` on a filename it can't parse — it doesn't skip silently.

### Group or filter a file list by burst ID or date

```python
opera_utils.group_by_burst(files)             # {burst_id: [files, ...]}
opera_utils.filter_by_burst_id(files, "t087_185683_iw2")
opera_utils.sort_by_burst_id(files, opera_utils.OPERA_BURST_RE)

opera_utils.get_dates("S1A_..._20191231T000000_....nc")   # [datetime(2019, 12, 31)]
opera_utils.group_by_date(files)
opera_utils.sort_files_by_date(files)
```

### Look up the bursts, bounding box, or EPSG for a DISP-S1 frame

```python
opera_utils.get_burst_ids_for_frame(11114)      # list[str]
opera_utils.get_frame_bbox(11114)               # (epsg, Bbox)
opera_utils.get_frame_to_burst_mapping(11114)   # full dict: epsg, xmin/ymin/xmax/ymax, is_land, burst_ids
```

```bash
opera-utils disp-s1-frame-bbox 11114
opera-utils disp-s1-frame-bbox 11114 --latlon --bounds-only
```

### Find which frame(s) a burst belongs to

Most bursts belong to exactly one frame; bursts in the along-track overlap between two adjacent frames belong to two:

```python
opera_utils.get_frame_ids_for_burst("t001_000009_iw1")   # e.g. [1, 2]
```

### Find frames that intersect an area of interest

*Requires `opera-utils[geopandas]`.*

```python
from opera_utils._types import Bbox

opera_utils.get_intersecting_frames(Bbox(-114.1, 31.5, -113.1, 32.5))  # GeoJSON FeatureCollection
```

```bash
opera-utils disp-s1-intersects --bbox -114.1 31.5 -113.1 32.5 --ids-only
opera-utils disp-s1-intersects --point -114 31
```

### Get frame or burst geometries as a GeoDataFrame

*Requires `opera-utils[geopandas]`.*

```python
gdf_frames = opera_utils.get_frame_geodataframe(frame_ids=[11114, 11115])
gdf_bursts = opera_utils.get_burst_geodataframe(burst_ids=["t087_185683_iw2"])
```

### Find spatially-consistent burst/date subsets before building a stack

Real archives often have incomplete (burst_id, date) coverage. Mixing burst IDs with different date coverage in one stack introduces spatial discontinuities in the resulting time series product. `get_missing_data_options` finds the maximal *rectangular* subsets (same dates for every burst kept):

```python
options = opera_utils.get_missing_data_options(slc_files=file_list)
best = options[0]  # sorted by total burst count, descending
print(best.burst_ids, best.dates, best.total_num_bursts)
```

```bash
opera-utils disp-s1-missing-data-options my_cslc_urls.txt --max-options 3
```

This writes one filtered filename list per option (`option_0_bursts_..._dates_....txt`) that you can feed into a download or processing step.

## CSLC / RTC search and download (ASF)

*Requires `opera-utils[asf]`.*

### Search for CSLC or RTC products

```python
from opera_utils.download import search_cslcs, L2Product

results = search_cslcs(
    bounds=(-102.71, 31.35, -102.6, 31.45),   # lon/lat WSEN
    start="2022-01-01", end="2022-06-01",
    product=L2Product.CSLC,
)
```

`start`/`end` accept `datetime`s, ISO strings, or natural-language strings like `"3 weeks ago"` (an `asf_search` feature). Duplicate reprocessed versions of the same `(startTime, burst_id)` are automatically filtered down to the latest.

```bash
opera-utils search-l2 --disp-s1-frame-id 20697 --start 2022-01-01 --end 2022-06-01
```

### Download CSLC, RTC, or static-layer products

```python
from opera_utils.download import download_cslcs, download_rtcs, download_cslc_static_layers

paths = download_cslcs(
    burst_ids=["t087_185683_iw2", "t087_185682_iw2"],
    output_dir="cslcs/",
    start="2022-10-01", end="2023-03-29",
)
```

These require a `~/.netrc` entry for `urs.earthdata.nasa.gov` (see [credentials setup](#set-up-earthdata-and-aws-credentials)).

## DISP-S1 search, download, and processing

*Requires `opera-utils[disp]`.*

### Search CMR for DISP-S1 products

```python
from datetime import datetime
from opera_utils.disp import search

products = search(frame_id=11115, end_datetime=datetime(2024, 1, 1))
```

```bash
opera-utils disp-s1-search --frame-id 11115 --end-datetime 2024-01-01
```

Omitting `frame_id` searches the whole archive and emits a warning — always scope by frame when you can.

### Download and spatially subset DISP-S1 products

```python
from opera_utils.disp._download import run_download

paths = run_download(
    frame_id=20697,
    bbox=(-102.71, 31.35, -102.6, 31.45),
    end_datetime=datetime(2024, 1, 1),
    num_workers=4,
    output_dir="subsets-west-texas",
)
```

```bash
opera-utils disp-s1-download \
    --frame-id 20697 --bbox -102.71 31.35 -102.6 31.45 \
    --end-datetime 2024-01-01 --num-workers 4 --output-dir subsets-west-texas
```

### Open a DISP-S1 product or remote HDF5 file directly

```python
from opera_utils.disp import DispProduct, open_h5

product = DispProduct.from_filename("OPERA_L3_DISP-S1_IW_F11116_VV_20160705T140755Z_20160729T140756Z_v1.0_20250318T222753Z.nc")
product.epsg, product.shape, product.bounds

with open_h5(product.filename) as hf:   # works for local, https://, and s3:// paths
    value = hf["displacement"][4000, 4000]
```

### Reformat a stack into one analysis-ready NetCDF/Zarr file

This is the step that rebases the moving reference date, applies corrections, masks low-quality pixels, and spatially re-references every epoch — see [Background theory](background-theory.md#why-disp-s1-stacks-need-rebasing) for why it's needed.

```python
from pathlib import Path
from opera_utils.disp import reformat_stack
from opera_utils.disp._enums import ReferenceMethod

reformat_stack(
    input_files=sorted(Path("subsets-west-texas").glob("OPERA_L3_DISP-S1*.nc")),
    output_name="stack-west-texas.zarr",   # or .nc
    reference_method=ReferenceMethod.BORDER,   # NONE | POINT | MEDIAN | BORDER | HIGH_COHERENCE
)
```

```bash
opera-utils disp-s1-reformat \
    --input-files subsets-west-texas/OPERA_L3_DISP-S1*.nc \
    --output-name stack-west-texas.zarr \
    --reference-method BORDER
```

Then open it normally:

```python
import xarray as xr
ds = xr.open_dataset("stack-west-texas.zarr", consolidated=False)   # or engine="h5netcdf" for .nc
```

### Choose a spatial reference method

Pass one of these to `reformat_stack`'s `reference_method`:

| Method | Behavior |
| --- | --- |
| `POINT` | Subtract the value at one `(row, col)` or `(lon, lat)` pixel. |
| `MEDIAN` | Subtract the spatial median over the whole frame, per epoch. |
| `BORDER` | Subtract the median of an outer ring of pixels (`reference_border_pixels` wide). |
| `HIGH_COHERENCE` | Subtract the median over pixels above a coherence threshold (`reference_coherence_threshold`). |
| `NONE` | No spatial re-referencing. |

### Rebase a raw displacement array yourself

If you need the moving-reference-date correction as a standalone step (e.g. on a custom array), use the low-level function directly:

```python
from opera_utils.disp import create_rebased_displacement

da_rebased = create_rebased_displacement(
    ds["displacement"], reference_datetimes=stack.reference_dates,
)
```

or, for a plain numpy `(time, rows, cols)` array:

```python
from opera_utils.disp._rebase import rebase_timeseries

rebased = rebase_timeseries(raw_array, reference_dates)
```

### Build a single-reference stack as per-epoch GeoTIFFs instead

For tools that prefer per-date rasters over one NetCDF/Zarr cube (e.g. GDAL/QGIS workflows, or chaining into [`dolphin`](https://github.com/isce-framework/dolphin)):

```python
from pathlib import Path
from opera_utils.disp.rebase_reference import main

main(nc_files=list(Path("subsets-west-texas").glob("*.nc")), output_dir=Path("aligned"))
```

```bash
python -m opera_utils.disp.rebase_reference aligned/ subsets-west-texas/OPERA_L3_DISP-S1*.nc
```

This writes per-epoch `displacement_{ref}_{sec}.tif` / `short_wavelength_displacement_{ref}_{sec}.tif`, plus quality-layer rasters and their averages, and auto-selects a high-quality, spatially-central reference pixel if you don't supply one.

### Correct for tectonic plate motion

`plate_motion.py` computes rigid ITRF2014 plate velocity at each pixel from a line-of-sight ENU raster and projects it onto LOS, as a standalone correction raster (not applied automatically by `reformat_stack`):

```bash
uv run --script src/opera_utils/disp/plate_motion.py \
    --los-enu-path los_enu.tif --plate-name NorthAmerica --out plate_motion_in_los.tif
```

```python
from opera_utils.disp.plate_motion import EulerPole, ITRF2014_PMM

plate = ITRF2014_PMM["northamerica"]
euler = EulerPole(plate.omega_x, plate.omega_y, plate.omega_z)
ve, vn, vu = euler.velocity_enu_m_per_yr(lat_deg=32.0, lon_deg=-103.0)
```

### Export a reformatted stack to MintPy

```python
from pathlib import Path
from opera_utils.disp.mintpy import disp_nc_to_mintpy

disp_nc_to_mintpy(
    Path("stack-west-texas.zarr"),
    sample_disp_nc=Path("subsets-west-texas/OPERA_L3_DISP-S1_..._v1.0_....nc"),
    outdir=Path("mintpy"),
)
```

```bash
python -m opera_utils.disp.mintpy stack-west-texas.zarr \
    --sample-disp-nc subsets-west-texas/OPERA_L3_DISP-S1_....nc --outdir mintpy
```

Writes `timeseries.h5`, `avgSpatialCoh.h5`, `velocity.h5`, reliability masks, and (given a `los_enu_path` or `geometry_dir`) `geometryGeo.h5` — ready for standard MintPy tools.

### Derive an incidence-angle / slant-range raster from a LOS-ENU static layer

```bash
python -m opera_utils.disp.create_incidence_range --los-enu geometry/los_enu.tif
```

## NISAR GSLC search, download, and processing

*Requires `opera-utils[nisar]`.*

### Search for GSLC or GUNW products

```python
from opera_utils.nisar import search, search_gunw

gslcs = search(bbox=(40.62, 13.56, 40.72, 13.64))
gslcs = search(relative_orbit_number=172, track_frame_number=8, orbit_direction="A")

pairs = search_gunw(bbox=(40.62, 13.56, 40.72, 13.64))   # pre-formed interferogram pairs
```

```bash
opera-utils nisar-gslc-search --track-frame 076_A_022
opera-utils nisar-gunw-search --relative-orbit-number 151 --track-frame-number 11 --orbit-direction A
```

### Download and subset GSLC products

```python
from opera_utils.nisar import run_download

paths = run_download(
    bbox=(40.62, 13.56, 40.72, 13.64),
    polarizations=["HH"],
    output_dir="gslc_subsets",
    num_workers=4,
)
```

```bash
opera-utils nisar-gslc-download \
    --output-dir gslc_subsets --bbox -118.5 34.0 -118.0 34.5 \
    --polarizations HH --polarizations VV --num-workers 4
```

For a full-file download with no subsetting, use `download_gslcs` instead (needs `~/.netrc`, not the env-var fallback).

### Read a GSLC product

```python
from opera_utils.nisar import GslcProduct

product = GslcProduct.from_filename("gslc_subsets/NISAR_L2_PR_GSLC_004_076_A_022_...h5")
shape = product.get_shape(frequency="A", polarization="HH")
row, col = product.lonlat_to_rowcol(lon=40.65, lat=13.60)
subset = product.read_subset(rows=slice(row, row + 500), cols=slice(col, col + 500))
```

### Mask invalid or water pixels

```python
from opera_utils.nisar import get_gslc_mask, get_gunw_mask

mask = get_gslc_mask("gslc_subsets/NISAR_..._001.h5", output_path="gslc_valid_mask.tif")  # True = valid
mask = get_gunw_mask("NISAR_GUNW_....h5", layer="unwrappedInterferogram")
```

### Compute incidence angle / LOS geometry

*Additionally requires `rioxarray`, `rasterio`, `scipy`, `xarray`.*

```python
from opera_utils.nisar import prepare_incidence_angle

inc_path, los_path = prepare_incidence_angle(
    gslc_path="gslc_subsets/NISAR_..._001.h5",
    dem_path="dem.tif",
    ref_tif="ifg_multilooked.tif",
    n_workers=8,
)
```

### Look up NISAR frame info from a bounding box, frame number, or product file

*Additionally requires `opera-utils[geopandas]`; needs a `.gpkg` frame database from `fetch_nisar_frame_to_bounds_file`-adjacent releases.*

```python
from opera_utils.nisar import nisar_frame_info

nisar_frame_info(Path("nisar_frames.gpkg"), bbox=(40.6, 13.5, 40.8, 13.7))
nisar_frame_info(Path("nisar_frames.gpkg"), frame_number=8, plot=True)
```

```bash
opera-utils nisar-frame-info nisar_frames.gpkg --frame-number 8
```

## Tropospheric delay correction

*Requires `opera-utils[tropo]`.*

The workflow: crop raw OPERA TROPO products to your AOI and interpolate to your interferogram dates, project onto the ground/LOS geometry, then match and subtract from interferograms by date. See [Background theory](background-theory.md#tropospheric-delay-correction) for the LOS sign convention.

```python
from datetime import datetime
from pathlib import Path
from opera_utils.tropo import crop_tropo, apply_tropo, match_and_apply_tropo

crop_tropo(
    tropo_urls_file=Path("tropo_urls.txt"),
    datetimes=[datetime(2020, 1, 1, 3, 0, 0)],
    aoi_bounds=(-120.5, 34.0, -119.5, 35.0),
    output_dir=Path("cropped_tropo"),
)

apply_tropo(
    cropped_tropo_list=sorted(Path("cropped_tropo").glob("tropo_cropped_*.nc")),
    dem_path=Path("dem_warped_utm.tif"),
    incidence_angle_path=Path("incidence_angle.tif"),
    output_dir=Path("tropo_corrections"),
)

match_and_apply_tropo(
    ifg_files=sorted(Path("interferograms").glob("*.iono_corrected.tif")),
    tropo_dir=Path("tropo_corrections"),
    reference_point_file="reference_point.txt",
)
```

```bash
opera-utils tropo-crop --tropo-urls-file tropo_urls.txt --datetimes 2020-01-01T03:00:00 --aoi-bounds -120.5 34.0 -119.5 35.0
opera-utils tropo-apply --cropped-tropo-list cropped_tropo/*.nc --dem-path dem_warped_utm.tif --incidence-angle-path incidence_angle.tif
```

## Raster and HDF5 utilities

*Some functions below require `rasterio`/GDAL (part of `opera-utils[disp]`, or install GDAL separately).*

### Stitch/merge rasters across bursts or subswaths

Requires `gdal_merge.py` on `PATH` in addition to the `osgeo` Python bindings.

```python
from opera_utils import stitching

stitching.merge_images(
    file_list=["los_east_0.h5", "los_east_1.h5"],
    outfile="los_east.tif",
    strides={"x": 6, "y": 3},
    out_nodata=0,
)
```

### Build a CSLC static-layer geometry stack for a frame

```python
from opera_utils.geometry import create_geometry_files

files = create_geometry_files(frame_id=11114, output_dir="geometry_out/", max_download_jobs=4)
```

### Read raster metadata or arrays without opening a full GIS stack

```python
from opera_utils._io import load_gdal, get_raster_bounds, get_raster_crs

arr = load_gdal("stitched.tif", band=1, subsample_factor=2)
bounds = get_raster_bounds("stitched.tif")
```

### Explore an HDF5 product interactively in Jupyter

*Requires `ipywidgets` and `matplotlib`.*

```python
from opera_utils.h5explorer import HDF5Explorer

h = HDF5Explorer("OPERA_L2_CSLC-S1_....h5")
h.data.VV.shape   # tab-completable
```

## Command-line reference

`opera-utils --help` lists every subcommand available in your current install. Always-available (base install):

| Command | Wraps |
| --- | --- |
| `search-l2` | `opera_utils.download.search_l2` (requires `asf_search`) |
| `disp-s1-frame-bbox` | frame → EPSG/bbox lookup |
| `disp-s1-intersects` | frame ↔ AOI intersection (requires `geopandas`) |
| `disp-s1-missing-data-options` | spatially-consistent burst/date subset selection |

Added by extras: `disp-s1-download` / `disp-s1-search` / `disp-s1-reformat` (`[disp]`); `tropo-crop` / `tropo-apply` (`[tropo]`); `nisar-gslc-download` / `nisar-gslc-search` / `nisar-gunw-search` / `nisar-frame-info` (`[nisar]`).

Every subcommand's flags are generated from the underlying Python function's parameters (`snake_case` → `--kebab-case`), so `<command> --help` always reflects the true, current signature.
