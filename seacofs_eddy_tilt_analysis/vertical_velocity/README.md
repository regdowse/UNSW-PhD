# Vertical velocity snapshots

Open `01_snapshot_core_extrema.ipynb` on Katana. Set one or more `(Eddy, Day)` pairs in the controls cell. It reads the original SEACOFS `w` field (upward m/s) at `s_w` interfaces. This is distinct from the fitted eddy-table `w`, which is relative vorticity.

For each fitted profile depth, the notebook places the same maximum-tangential-velocity ellipse used by `seacofs_tilt_tools.core_grid_indices` at that depth's fitted centre. `FRACTIONS` scales its axes. Each wet rho cell uses the nearest local native `s_w` interface; samples exceeding `MAX_MISMATCH_M` are excluded. The table records valid coverage and the physical depths and locations of the largest upward and downward velocities. A second table reduces these to the largest upward and downward values over the selected column.

The maps and depth curves are meant to reveal whether extrema are spatially coherent and stable across core fractions. Single-cell maxima are deliberately exploratory; they can be sensitive to noise, boundary cells, and the number of cells searched. The notebook reports coverage and also shows the core mean for comparison. It does not calculate a population correlation with tilt yet.
