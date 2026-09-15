# Background theory

This page explains the concepts behind `opera-utils`'s API — what a burst ID and frame actually are, why DISP-S1 stacks need "rebasing" before you can treat them as a time series, and how the plate-motion and tropospheric corrections fit in. For task recipes, see the [How-to guides](how-to-guides.md); for exact function signatures, see the [Code reference](reference/summary.md).

## Bursts and frames

### Bursts

Sentinel-1's TOPS acquisition mode splits each pass into three sub-swaths (IW1, IW2, IW3), and each sub-swath into a sequence of **bursts** along track. OPERA (via the [`burst_db`](https://github.com/opera-adt/burst_db) project) assigns every `(track, along-track position, sub-swath)` combination a stable **burst ID** of the form

```
t<track>_<burst_number>_iw<subswath>
```

e.g. `t087_185683_iw2` — 3-digit relative orbit/track number, a 6-digit ESA-assigned burst index, and the sub-swath. The same ID also appears in the hyphenated, uppercase form used in official OPERA product filenames (`T087-185683-IW2`); `opera_utils.normalize_burst_id` collapses either form to the canonical one, and every other burst-ID function in the package normalizes internally, so you rarely need to think about which form you have.

Burst IDs are geographically stable across repeat passes — the same physical patch of ground is always `t087_185683_iw2`, acquisition after acquisition — which is what makes it possible to group multi-temporal stacks by burst.

### Frames

A **frame** is a fixed, pre-defined grouping of roughly two dozen burst IDs spanning a track's IW1–IW3 sub-swaths over a contiguous along-track extent. Frame boundaries are defined once — not per-acquisition — in the same `burst_db` reference database, identified by a single integer `frame_id`. OPERA generates one DISP-S1 product stack per frame, from a consistent set of that frame's constituent bursts across time.

Two facts follow directly from this design:

- **Adjacent frames overlap** by a few burst IDs at their shared edge, so a burst can belong to more than one frame. `opera_utils.get_frame_ids_for_burst` reflects this by returning a *list* (usually length 1, occasionally 2).
- **Frames carry their own CRS.** Each frame has a UTM EPSG code and a projected bounding box (`opera_utils.get_frame_bbox`) — DISP-S1 products are natively gridded in that per-frame UTM projection, not lat/lon.

`opera-utils` ships two complementary lookup tables (auto-downloaded and cached from `burst_db` GitHub releases via [`pooch`](https://www.fatiando.org/pooch/)):

- **frame → bursts**: burst ID list, EPSG, and bounding box per frame (`opera_utils.get_frame_to_burst_mapping`, `get_burst_ids_for_frame`, `get_frame_bbox`).
- **burst → frames**: the (usually singleton) list of frames each burst belongs to (`opera_utils.get_burst_to_frame_mapping`, `get_frame_ids_for_burst`).

Actual polygon footprints (not just bounding boxes) are available separately as GeoJSON/GeoDataFrames, letting `opera_utils.get_intersecting_frames` do a true spatial — not just bbox — intersection against an area of interest.

### The missing-data problem

Real acquisition archives are rarely complete for every `(burst_id, date)` pair — a burst might be missing for one pass due to a processing gap or an orbit change. If you build a stack using different sets of dates for different burst IDs within a frame, the resulting product has **spatial discontinuities**: some parts of the frame have data on dates that other parts don't. `opera_utils.get_missing_data_options` searches for the maximal *rectangular* subsets of a burst/date incidence matrix — same date set used for every burst ID kept — so you can choose a spatially-consistent stack before processing, trading off total burst count against coverage completeness.

## The DISP-S1 product

DISP-S1 is OPERA's Level-3 surface displacement product: one file per frame per date-pair, produced from a **time-series InSAR** algorithm (small-baseline-style) run incrementally as new acquisitions arrive.

### Why DISP-S1 stacks need rebasing

Because DISP-S1 is generated incrementally, the algorithm doesn't keep every historical acquisition as the reference forever — periodically, the reference date jumps forward to a more recent acquisition (each such segment is a "ministack"). Practically, that means a raw sequence of DISP-S1 products is **not** already referenced to a single, fixed date: epoch *N* might be displacement relative to `2021-01-01`, while epoch *N+1* is relative to `2021-06-15`.

`opera_utils.disp`'s rebasing step turns this into a single continuous time series by summing the "crossover" displacement each time the reference date changes:

$$
d_i^{\text{rebased}} = d_i^{\text{raw}} + \sum_{\substack{j \le i \\ \text{ref}(j) \ne \text{ref}(j-1)}} d_{j-1}^{\text{rebased}}
$$

Intuitively: every time the reference jumps, add back in the displacement that had already accumulated up to the old reference, so the running total stays anchored to the very first epoch. This is implemented as a plain-numpy loop in `opera_utils.disp._rebase.rebase_timeseries` (operating on a `(time, rows, cols)` array) and as an xarray/dask-friendly, chunked wrapper in `opera_utils.disp.create_rebased_displacement`. `reformat_stack` (the main pipeline entry point) applies this automatically; `rebase_reference.py` implements the same accumulation for the alternative per-epoch-GeoTIFF workflow.

A `NaNPolicy` controls how missing data is handled when it falls on a reference-crossover epoch: `propagate` (default — a NaN there poisons every subsequent epoch, since the true cumulative offset is unknown) or `omit` (treat it as zero displacement).

### Spatial referencing

Displacement is inherently *relative* — every pixel's value only means something in comparison to a reference. After rebasing (which fixes the *temporal* reference), you also need to choose a *spatial* reference: which pixel's value should be treated as "no motion." `opera_utils.disp._reference.get_reference_values` implements four strategies (`ReferenceMethod`):

- **`POINT`** — a single, user-chosen pixel (by row/col or lon/lat) assumed to be stable.
- **`MEDIAN`** — the spatial median across the whole frame, per epoch (a "most pixels are stable" assumption).
- **`BORDER`** — the median of an outer ring of pixels, useful when deformation is expected to be concentrated away from the frame edges.
- **`HIGH_COHERENCE`** — the median over pixels above a temporal-coherence threshold, i.e. only "trustworthy" pixels vote.

### Corrections applied during reformatting

Beyond rebasing and referencing, `reformat_stack` optionally applies the solid-earth-tide and ionospheric-delay correction layers that ship inside each DISP-S1 product, and masks pixels using one or more quality layers (e.g. `recommended_mask`, `phase_similarity`) combined with a boolean reduction (`logical_or`/`logical_and` across layers) via `combine_quality_masks`.

Two corrections are **not** applied automatically and are left as separate, explicit steps:

- **Tectonic plate motion** — see below.
- **Tropospheric delay** — see below.

## Tectonic plate motion correction

InSAR line-of-sight displacement includes the *rigid* motion of the tectonic plate the ground sits on, which is usually not the signal of interest (it's a near-constant secular rate, well described by existing plate models, as opposed to the local deformation signal being studied). `opera_utils.disp.plate_motion` computes this rigid contribution using the ITRF2014 plate motion model: each plate's motion is described by an **Euler pole** — an angular velocity vector $(\omega_x, \omega_y, \omega_z)$ such that a point at Earth-centered position $\vec{r}$ moves with velocity

$$
\vec{v} = \vec{\omega} \times \vec{r}
$$

`EulerPole.velocity_enu_m_per_yr` evaluates this at a given lat/lon and converts to local East-North-Up components, which are then projected onto the sensor's line-of-sight direction to produce a correction raster in the same units as the displacement product. This is deliberately a separate, standalone script/module rather than a step inside `reformat_stack`, since it depends on choosing a plate (not always obvious near plate boundaries) and a LOS geometry raster.

## Tropospheric delay correction

Radar signals travel more slowly through denser/wetter atmosphere, adding a spurious phase delay correlated with topography and weather that varies from acquisition to acquisition. `opera_utils.tropo` corrects for this using OPERA TROPO products — pre-computed zenith hydrostatic + wet delay cubes on a `(time, height, lat, lon)` grid (derived from numerical weather models).

The correction pipeline mirrors the physical picture:

1. **Crop + interpolate in time** (`crop_tropo`): TROPO products are only available at fixed synoptic times (every few hours), so each requested acquisition datetime is bracketed by the nearest available products (within ±6 hours) and the delay is linearly interpolated in time.
2. **Project to the ground and to line-of-sight** (`apply_tropo`): the 3D zenith-delay cube is interpolated onto the DEM surface (so "zenith delay at this ground point's actual elevation," not just at a fixed height), then converted from zenith to line-of-sight using the incidence angle $\theta$:

$$
\text{delay}_{\text{LOS}} = -\,\frac{\text{delay}_{\text{zenith}}}{\cos\theta}
$$

   The sign convention matches DISP-S1: a *positive* correction means a *decrease* in path delay, which reads as apparent uplift toward the satellite.

3. **Match and subtract** (`match_and_apply_tropo`): each interferogram is matched to the tropospheric correction for its *secondary* date (by filename convention) and the correction is subtracted, optionally after re-referencing the correction itself to the same pixel used as the interferogram's spatial reference — otherwise the correction would reintroduce a spatially-constant offset that spatial referencing was meant to remove.

## Data access architecture

Across CSLC, DISP-S1, and NISAR GSLC, the package uses the same two building blocks for remote access:

- **Search via NASA CMR** (`opera_utils._cmr.fetch_cmr_pages`) — an unauthenticated, paginated UMM-JSON query against `cmr.earthdata.nasa.gov`, used by `opera_utils.disp.search`, `opera_utils.nisar.search`, and `opera_utils.nisar.search_gunw`. No credentials are needed just to search; they're only needed to actually fetch product bytes.
- **Cloud-optimized HDF5 reads** (`opera_utils._remote.open_h5`) — rather than downloading a whole file before reading it, `open_h5` wraps an `fsspec` byte stream (HTTPS or S3) in `h5py.File` with paged, cached reads, so downstream code can subset/stream a large remote HDF5 product much like a local one. This is what makes spatial subsetting during download (`run_download`'s `process_file`) practical: only the requested slice's bytes are actually transferred.

Both are backed by the shared `opera_utils.credentials` module (Earthdata Login resolution, and — for S3 — either static AWS credentials or short-lived ones obtained from ASF's `s3credentials` endpoint using your Earthdata Login).

## NISAR product identifiers

NISAR GSLC and GUNW filenames encode a **track/frame ID** in the form `RRR_D_TTT` (relative orbit, `A`/`D` orbit direction, frame number) — conceptually analogous to a DISP-S1 frame ID, but NISAR's own fixed frame grid rather than OPERA's `burst_db` one. GUNW filenames additionally encode a reference/secondary acquisition-date pair (`pair_id`, `"YYYYMMDD_YYYYMMDD"`), the same "which two acquisitions form this interferogram" idea as a DISP-S1 date pair.
