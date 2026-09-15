# opera-utils

`opera-utils` is a Python toolkit for working with data products from the [NASA JPL OPERA project](https://www.jpl.nasa.gov/go/opera) and related NISAR products, built around three things:

- **Filename and metadata parsing** for OPERA CSLC-S1, DISP-S1, and NISAR GSLC/GUNW products (burst IDs, frame IDs, acquisition dates, orbit/polarization info).
- **The DISP-S1 frame/burst database**, so you can go from a frame ID to its constituent Sentinel-1 bursts (and back), or find which frames cover an area of interest, without re-deriving OPERA's tiling scheme yourself.
- **Search, download, and processing helpers** for pulling CSLC/RTC/DISP-S1/NISAR-GSLC products from ASF and NASA CMR, subsetting them to an area of interest, and turning a stack of individual granules into an analysis-ready time series (rebased to a single reference date, spatially re-referenced, corrected for solid-earth/ionospheric/tropospheric/plate-motion effects, and optionally exported to [MintPy](https://github.com/insarlab/MintPy)).

Most functionality lives in the base package; heavier, format-specific pipelines (DISP-S1, NISAR GSLC, tropospheric correction) are opt-in via [extras](getting-started.md#install) so you only pull in `xarray`/`dask`/`rasterio`/etc. when you need them.

## Where to start

This documentation follows the [Diátaxis](https://diataxis.fr/) framework — pick the page that matches what you're trying to do:

| I want to... | Go to |
| --- | --- |
| Install the package and run my first commands | [Getting started](getting-started.md) |
| Follow a guided, end-to-end walkthrough | [Tutorials](tutorials.md) |
| Look up how to do one specific task | [How-to guides](how-to-guides.md) |
| Look up an exact function signature | [Code reference](reference/summary.md) |
| Understand *why* things work the way they do (bursts, frames, reference-date rebasing, plate motion) | [Background theory](background-theory.md) |

## Quick links

- Source & issue tracker: [github.com/opera-adt/opera-utils](https://github.com/opera-adt/opera-utils)
- PyPI: [pypi.org/project/opera-utils](https://pypi.org/project/opera-utils/)
- conda-forge: [github.com/conda-forge/opera-utils-feedstock](https://github.com/conda-forge/opera-utils-feedstock)
- Command-line entry point: `opera-utils --help`
