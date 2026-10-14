# CSMX-WRF

WRF scripts and evaluation code for cyclic forcing of the ML-based CSMX soil moisture dataset in convection-permitting simulations over eastern China.

## Overview

This repository contains the workflow used to drive the Weather Research and Forecasting (WRF) model with **CSMX**, a machine-learning-based soil moisture dataset, through *cyclic (restart) forcing*. Soil moisture is re-initialized from CSMX at regular intervals during long convection-permitting integrations, and the resulting simulations are evaluated against observations to assess the impact on land–atmosphere coupling over eastern China.

## Repository contents

| File | Description |
| --- | --- |
| `recycle.py` | Core driver that reads CSMX soil moisture, regrids it onto the WRF grid, overwrites `SMOIS`/`SH2O` in the restart file, updates `namelist.input`, and submits/monitors the next WRF segment on a job scheduler (SLURM/LSF). |
| `namelist.input` | Example WRF configuration (4 km convection-permitting domain, `sf_surface_physics = 4`, 4 soil layers) used for the simulations. |
| `SM.ipynb` | Soil moisture evaluation against in-situ observations and reference datasets (ERA5-Land, GLDAS-Noah, CLDAS). |
| `Prec.ipynb` | Precipitation evaluation against CMORPH (climatology, bias, and Hovmöller/zonal-mean analysis). |
| `SurfaceClimatology.ipynb` | Climatological evaluation of surface variables (T2, Q2, ET, SH, wind, precipitation) with Taylor diagrams and skill metrics. |
| `SurfaceExtremeClimatology.ipynb` | Extreme climatology analysis (TXx, TNn, RX1day) with exceedance probability, Mann–Whitney U tests, and year-block bootstrap. |
| `diagnosis.ipynb` | Sub-regional diagnosis of land–surface fluxes (SM, SH, LH, runoff) and their response to the cyclic forcing. |
| `WB.ipynb` | Water-balance (moisture recycling) analysis for the sub-regions. |

## Workflow

1. **Cyclic forcing** (`recycle.py`): at every forcing interval (default 48 h), the CSMX soil moisture for the current date is interpolated to the WRF grid (nearest-neighbor), quality-controlled against `smcmax`, and written into the top soil layer of the WRF restart file (`wrfrst_d01`). The time-control section of `namelist.input` is then rewritten for the next segment and the job is resubmitted until the full simulation period is covered.

2. **Evaluation**: the Jupyter notebooks compare a control run (CTRL) against the CSMX-forced run (FSM) using reference observations (GLEAM, CMFD, CMORPH, ERA5-Land, GLDAS, CLDAS) across five sub-regions of eastern China (HHP, LPT, QDW, YHP, SEH).

## Requirements

- Python 3 with `numpy`, `pandas`, `scipy`, `netCDF4`, `xarray`, `matplotlib`, `cartopy`, `cmaps`, `shapefile`, `statsmodels`, `joblib`, `tqdm`
- A compiled WRF installation and the WPS static files (`geo_em.d01.nc`, `soil_properties.nc`)
- A job scheduler (SLURM/LSF) configured for the target machine

## Usage

Run the cyclic-forcing workflow:

```bash
nohup python -u recycle.py > recycle.out 2>&1 &
```

Before running, set the paths and experiment parameters at the bottom of `recycle.py` (`DOMAIN_PATH`, `RUN_PATH`, `SM_PATH`, `Replace_Interval_Hours`, start/end dates). The scheduler commands (`sbatch`/`bjobs` or `scancel`) must be adapted to the local system.

## Citation

If you use this code, please cite the associated paper (IJOC submission).
Dong, Y., Zhang, G., Wang, S., Yuan, H., Guo, W., & Tang, J. (2026). Impacts of machine learning-based soil moisture data on convective-permitting regional climate simulations over eastern China. *International Journal of Climatology*. Manuscript submitted for publication.
