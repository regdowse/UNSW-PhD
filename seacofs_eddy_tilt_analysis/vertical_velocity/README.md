# Vertical velocity snapshots

Open `01_snapshot_core_extrema.ipynb` on Katana. Set one or more `(Eddy, Day)` pairs in the controls cell. It reads the original SEACOFS `w` field (upward m/s) at `s_w` interfaces. This is distinct from the fitted eddy-table `w`, which is relative vorticity.

For each fitted profile depth, the notebook places the same maximum-tangential-velocity ellipse used by `seacofs_tilt_tools.core_grid_indices` at that depth's fitted centre. `FRACTIONS` scales its axes. As in the modular dataset pipeline, the native `w(s_w, eta_rho, xi_rho)` snapshot is transposed to `(x, y, z)` and flipped in `z`. The extra surface level is dropped from `w`, leaving 30 levels aligned directly with the 30 levels in `grid.z_r`. Each wet rho cell uses the closest local `z_r` depth; samples exceeding `MAX_MISMATCH_M` are excluded. The table records valid coverage and the reference depths and locations of the largest upward and downward velocities. A second table reduces these to the largest upward and downward values over the selected column.

The maps and depth curves are meant to reveal whether extrema are spatially coherent and stable across core fractions. Single-cell maxima are deliberately exploratory; they can be sensitive to noise, boundary cells, and the number of cells searched. The notebook reports coverage and also shows the core mean for comparison. It does not calculate a population correlation with tilt yet.

## Climatology

Run `02_vertical_velocity_climatology.ipynb` on Katana from this folder to compare
AE and CE vertical velocity in the upper 1,000 m. It selects a reproducible,
polarity-balanced subset by default; the controls permit all eligible eddies.
The primary outcome is the signed core mean at each fitted depth. It also
retains depth-specific RMS, maxima, minima and upward-cell fraction, plus one
sample-weighted mean, RMS, maximum and minimum for each eddy-day and core
fraction. Repeated days are averaged within each eddy before population
comparisons; bootstrap intervals resample eddies. Coverage and core-cell
thresholds are applied before analysis.

The notebook reads each selected model-day velocity volume once and reuses it
across fitted depths. Its one Parquet cache and settings JSON live under
`/srv/scratch/z5297792/SEACOFS_26yr_eddy_dataset_modular/vertical_velocity_climatology/`.
`REBUILD=True` overwrites the cache after settings change. The column mean is
an arithmetic mean over the sampled depth-cell values, not a volume flux.

Run `03_tilt_vertical_velocity.ipynb` after notebook 02. It uses the saved
per-depth cache and fixes the footprint at core fraction 1.5. It plots vertical
velocity profiles for a reproducible sample of individual eddy-days, then
compares each day's `TiltDis` and `TiltDir` with its upper-1,000 m signed core
mean and largest absolute vertical velocity. Upward and downward extrema are
retained separately. The notebook uses descriptive eddy-day scatterplots and
direction sectors; it does not treat repeated days as independent evidence.
No model NetCDF files are reread.
