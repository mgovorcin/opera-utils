# Tutorials

These tutorials walk through complete, realistic workflows end to end. If you just need a single function or command, see the [How-to guides](how-to-guides.md) instead.

## Tutorial 1 — From a Sentinel-1 filename to a DISP-S1 frame

*Requires only the base install — no extras, no credentials.*

Say you have a directory of OPERA CSLC-S1 files and want to know which DISP-S1 frame(s) they belong to, and what area that covers.

**1. Parse the burst ID out of each filename:**

```python
import opera_utils

files = [
    "OPERA_L2_CSLC-S1_T087-185683-IW2_20230322T161649Z_20240504T185235Z_S1A_VV_v1.1.h5",
    "OPERA_L2_CSLC-S1_T087-185684-IW2_20230322T161649Z_20240504T185235Z_S1A_VV_v1.1.h5",
]
burst_ids = [opera_utils.get_burst_id(f) for f in files]
# ['t087_185683_iw2', 't087_185684_iw2']
```

**2. Group files by burst, or by acquisition date** — useful before building a stack:

```python
by_burst = opera_utils.group_by_burst(files)
by_date = opera_utils.group_by_date(files)
```

**3. Find which frame(s) a burst belongs to** (most bursts belong to one frame; bursts in the overlap between two adjacent frames belong to two):

```python
>>> opera_utils.get_frame_ids_for_burst("t087_185683_iw2")
[11115]
```

**4. Go the other direction — get every burst ID and the bounding box for a frame:**

```python
>>> burst_ids = opera_utils.get_burst_ids_for_frame(11115)
>>> epsg, bbox = opera_utils.get_frame_bbox(11115)
```

**5. Find every frame that intersects an area of interest** (requires `opera-utils[geopandas]`):

```python
from opera_utils._types import Bbox

aoi = Bbox(left=-114.1, bottom=31.5, right=-113.1, top=32.5)
matches = opera_utils.get_intersecting_frames(aoi)  # GeoJSON FeatureCollection
```

or from the shell:

```bash
opera-utils disp-s1-intersects --bbox -114.1 31.5 -113.1 32.5 --ids-only
```

**6. Before processing a real stack, check for spatially-inconsistent coverage.** Real archives are often missing some (burst, date) pairs; mixing burst IDs with different date coverage introduces spatial discontinuities in a DISP-S1-style product:

```python
options = opera_utils.get_missing_data_options(slc_files=files)
best = options[0]  # sorted by total burst count, descending
print(best.num_burst_ids, "bursts x", best.num_dates, "dates")
```

See [Background theory](background-theory.md#bursts-and-frames) for why bursts, frames, and this consistency check exist in the first place.

## Tutorial 2 — Build and plot a DISP-S1 displacement time series

*Requires `pip install "opera-utils[disp]"` and an Earthdata Login (see [Getting started](getting-started.md#set-up-data-access-credentials)).*

This walks through the same pipeline as [`scripts/create-mintpy.sh`](https://github.com/opera-adt/opera-utils/blob/main/scripts/create-mintpy.sh) and the west-Texas example in the project README, using the Python API directly.

**1. Search CMR for DISP-S1 products over a frame:**

```python
from datetime import datetime
from opera_utils.disp import search

products = search(
    frame_id=20697,
    start_datetime=datetime(2021, 1, 1),
    end_datetime=datetime(2022, 1, 1),
)
print(len(products), "products found")
```

**2. Download and spatially subset them to your area of interest:**

```python
from opera_utils.disp._download import run_download

paths = run_download(
    frame_id=20697,
    bbox=(-102.71, 31.35, -102.6, 31.45),
    start_datetime=datetime(2021, 1, 1),
    end_datetime=datetime(2022, 1, 1),
    num_workers=4,
    output_dir="subsets-west-texas",
)
```

or, equivalently, from the command line:

```bash
opera-utils disp-s1-download \
    --frame-id 20697 \
    --bbox -102.71 31.35 -102.6 31.45 \
    --start-datetime 2021-01-01 \
    --end-datetime 2022-01-01 \
    --num-workers 4 \
    --output-dir subsets-west-texas
```

**3. Reformat the stack into one analysis-ready file.** This is the step that does the heavy lifting: it rebases the displacement time series onto a single reference date (DISP-S1's reference date otherwise jumps forward periodically — see [Background theory](background-theory.md#why-disp-s1-stacks-need-rebasing)), applies solid-earth/ionospheric corrections, masks low-quality pixels, and spatially re-references every epoch:

```python
from opera_utils.disp import reformat_stack
from opera_utils.disp._enums import ReferenceMethod
from pathlib import Path

reformat_stack(
    input_files=sorted(Path("subsets-west-texas").glob("OPERA_L3_DISP-S1*.nc")),
    output_name="stack-west-texas.zarr",
    reference_method=ReferenceMethod.BORDER,
)
```

CLI equivalent:

```bash
opera-utils disp-s1-reformat \
    --input-files subsets-west-texas/OPERA_L3_DISP-S1*.nc \
    --output-name stack-west-texas.zarr \
    --reference-method BORDER
```

**4. Open it with xarray and plot:**

```python
import xarray as xr
import matplotlib.pyplot as plt

ds = xr.open_dataset("stack-west-texas.zarr", consolidated=False)
ds.displacement.isel(time=-1).plot()
plt.show()

# a single pixel's time series, in cm
(ds.displacement.isel(y=500, x=500) * 100).plot()
```

**5. (Optional) Export to MintPy** for standard InSAR time-series tools (velocity fitting, plotting, further corrections):

```python
from pathlib import Path
from opera_utils.disp.mintpy import disp_nc_to_mintpy

disp_nc_to_mintpy(
    Path("stack-west-texas.zarr"),
    sample_disp_nc=sorted(Path("subsets-west-texas").glob("*.nc"))[0],
    outdir=Path("mintpy"),
)
```

For more detail on each step — rebasing math, reference methods, plate-motion correction, tropospheric correction, and the alternative per-epoch GeoTIFF workflow — see the [DISP-S1 section of the How-to guides](how-to-guides.md#disp-s1-search-download-and-processing).

## Tutorial 3 — Form a NISAR GSLC interferogram

*Requires `pip install "opera-utils[nisar]"` and an Earthdata Login.*

This follows [`docs/examples/nisar_gslc_interferogram.py`](https://github.com/opera-adt/opera-utils/blob/main/docs/examples/nisar_gslc_interferogram.py).

**1. Search for two acquisitions over the same track/frame:**

```python
from opera_utils.nisar import search

products = search(bbox=(40.62, 13.56, 40.72, 13.64))
for p in products:
    print(p.filename, p.track_frame_id, p.start_datetime)
```

**2. Download and subset both to your area of interest:**

```python
from opera_utils.nisar import run_download

paths = run_download(
    bbox=(40.62, 13.56, 40.72, 13.64),
    polarizations=["HH"],
    output_dir="gslc_subsets",
    num_workers=1,
)
```

**3. Read the complex SLC data for each acquisition and form the interferogram:**

```python
import h5py
import numpy as np

with h5py.File(paths[0]) as f0, h5py.File(paths[1]) as f1:
    slc1 = f0["/science/LSAR/GSLC/grids/frequencyA/HH"][:]
    slc2 = f1["/science/LSAR/GSLC/grids/frequencyA/HH"][:]

ifg = slc1 * np.conj(slc2)
```

**4. Multilook, then derive amplitude, phase, and coherence** — see the full example script for the multilooking helper and pixel-spacing-aware look-count calculation, and plotting.

**5. Mask invalid pixels** before further analysis:

```python
from opera_utils.nisar import get_gslc_mask

mask = get_gslc_mask(paths[0])   # True = valid
```

See [How-to guides](how-to-guides.md#nisar-gslc-search-download-and-processing) for masking GUNW products, computing incidence-angle/LOS geometry, and searching pre-formed GUNW interferogram pairs directly (skipping GSLC-level interferogram formation).

## Notebooks

The repository also ships runnable notebooks under [`docs/notebooks/`](https://github.com/opera-adt/opera-utils/tree/main/docs/notebooks) and `docs/`:

- [`notebooks/tutorial-disp-s1-stack.ipynb`](https://github.com/opera-adt/opera-utils/blob/main/docs/notebooks/tutorial-disp-s1-stack.ipynb) — the DISP-S1 workflow above, interactively.
- `accessing-disp-s1-with-xarray.ipynb`, `demo-fetch-disp-s1.ipynb` — earlier exploratory notebooks.

!!! warning
    `accessing-disp-s1-with-xarray.ipynb` and `demo-fetch-disp-s1.ipynb` predate the current `opera_utils.disp` API (they reference a `stack_to_dataarray` helper and a `scripts/fetch_disp.py` script that no longer exist in this version). Use Tutorial 2 above, or `notebooks/tutorial-disp-s1-stack.ipynb`, instead — they exercise the current, tested API.
