# Eddy tilt mechanism analysis

These notebooks test why the already-measured EAC eddies tilt. `TiltDis` and `TiltDir` are loaded as authoritative measurements and are never recomputed.

Run order:

1. `00_build_stratification_cache.ipynb`
2. `01_vertical_shear.ipynb`
3. `02_beta_stratification_burger.ipynb`
4. `03_pv_topographic_steering.ipynb`
5. `04_depth_dependent_propagation.ipynb`
6. `05_wind_ekman_sensitivity.ipynb` (optional; requires a wind-stress cache)
7. `06_joint_mechanism_comparison.ipynb`
8. `07_temporal_tilt_evolution.ipynb`
9. `08_stratification_polarity_tilt.ipynb`

The stratification cache is file-parallel and restartable. Version 4 calculates
N2 from surface-referenced potential density; the retired version 3 used
pressure-dependent in-situ density and must not be used. For every eddy-day,
version 4 calculates depth-mean N2 at each ocean grid column inside the same
elliptical maximum-tangential-velocity contour used by `compute_core_mean`,
then takes the unweighted core mean. Centre-column values and core coverage
diagnostics are retained for sensitivity and quality control. The cache also
contains fixed-depth N2 integrals and maxima, pycnocline-window mean and maximum
N2, pycnocline depth, and density-threshold mixed-layer depth. It can be run
from the notebook or from Katana with:

```bash
python build_stratification_cache.py --workers 5 --point-batch-size 128
```

Start with 4–8 workers and check job memory before increasing concurrency.
The point batch size controls vectorized core-column reads and is not a Dask grid
chunk; the calculation never constructs full-domain xgcm metrics.
Each completed ROMS file is written to a partition folder whose name includes a
hash of the scientific settings. Changing the algorithm, requested depths,
pycnocline window, or MLD threshold therefore cannot silently reuse stale
partitions. The resulting cache is
`n2_eddy_day_v4_potential_density_core.parquet`.

`mechanism_tools.py` contains only helpers unique to this workflow. The notebooks reuse `seacofs_tilt_tools.py`, `ml_subsurface_tools.py`, and `beta_effect_background_flow/*` for existing functionality.

Notebooks 02, 03, and 06 also divide AEs and CEs separately into low, medium,
and high eddy-level core-N2 tertiles. Their primary directional comparison uses
N2 residuals from polarity-specific quadratic latitude trends, preventing the
equatorward stratification gradient from defining the categories by itself.
The resulting boxplots, rose plots, alignment summaries, heatmaps, and joint
model forest plots are all eddy-equal where directional observations repeat in
time.

The notebooks assume they are launched with this directory as the working directory on Katana. Cache paths are declared near the top of the relevant notebooks.

Notebook 07 tests whether measured tilt grows within eddies in absolute and
normalised age, then separates trajectories into sustained planetary-PV,
sustained topographic-PV, and transitional/mixed exposure. It includes
eddy-level slope tests, clustered nonlinear age models, increment/sign tests,
individual trajectories, PV-regime transitions, and a component-wise test of
the predicted surface-minus-deep differential displacement. It uses the
existing vertical-profile dictionary and does not require a new cache.

Notebook 08 isolates the stratification–polarity–tilt question. It compares
AE and CE core N2 over 0–200 m and 0–500 m, confirms the polarity contrast in
`TiltDis`, and separates raw, covariate-adjusted, between-eddy, and within-eddy
stratification associations. It also tests fixed-depth integrals and maxima,
pycnocline-following N2, pycnocline depth, and mixed-layer depth. It reuses the
core-N2 cache built by notebook 00.

After changing or rebuilding the stratification cache, rerun notebook 00 first,
then notebooks 02, 03, 06, and 08. Notebooks 01, 04, 05, and 07 do not consume
this cache and do not need to be rerun for an N2-method change.
