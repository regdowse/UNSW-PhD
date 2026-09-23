# Eddy and environmental controls of tilt

Run either notebook from this folder (or from inside `seacofs_eddy_tilt_analysis`) on Katana.

1. **eddy_tilt_controls.ipynb**: surface/depth-dependent relative vorticity, Omega, ellipse axis ratio, radius, rotation-speed scale, Rossby number; profile means/maxima/peak depths/endpoint changes; vertical extent, centreline shape and lifecycle; corrected core stratification; existing cached vertical velocity when available.
2. **environment_tilt_controls.ipynb**: log10 environmental PV-gradient magnitude and planetary/topographic components, their ratio, coherence, planetary beta, latitude, bathymetry/slope, existing bulk background shear, regional and seasonal sensitivity.

Both retain authoritative `TiltDis` and `TiltDir`. They show absolute and radius-normalised tilt, availability audits, eddy-day binned medians/IQRs with whole-eddy bootstrap intervals, equal-eddy correlation sensitivity, and adjusted associations with eddy-clustered errors. Major-axis orientation, propagation, along/cross-shear mechanism tests and depth-following PV extensions are excluded.

## Inputs and execution

Use your existing Katana environment with NumPy, pandas, SciPy, Matplotlib, statsmodels, netCDF4, IPython and a Parquet engine. There is no new ESP dependency or daily-model-file processing here. The established grid loader reads its reference grid.

Core inputs use `seacofs_tilt_tools.Paths`: processed Eddy-Day table, measured tilt and confirmed vertical profiles. Additional inputs:

- Corrected N2: `/srv/scratch/z5297792/SEACOFS_26yr_eddy_dataset/tilt_mechanisms/n2_eddy_day_v4_potential_density_core.parquet`. Version/method and per-metric coverage are validated. Core N2 is not background stratification.
- Vertical velocity: `/srv/scratch/z5297792/SEACOFS_26yr_eddy_dataset_modular/vertical_velocity_climatology/core_vertical_velocity_upper_1000m.parquet`. Reuse only; missing data skips this section explicitly. The existing builder is `../vertical_velocity/02_vertical_velocity_climatology.ipynb`; no new builder is added. Default fraction 1.5, at least 8 valid cells and 70% coverage.
- Environmental notebook: existing `tilt.DEFAULT_SURFACE_PV_CACHE`, validated against Gaussian core fraction 1, nonlinear averaging and peak-|vorticity| selection through 1000 m. It must already exist; stale settings cause a clear error instead of rebuilding or overwriting it. Cache-only environmental columns are joined to current tilt/property rows by Eddy-Day.
- Optional shear: existing `BackgroundConfig().background_table_path` in `background_flow_cache_all_eddies_v4`. This is the agreed shear comparison only. Missing cache skips the section.

A present but invalid cache fails rather than silently substituting another method. Optional missing sections are printed and recorded in output provenance. No claim is made that remote Katana caches exist based on local inspection.

Start with `N_BOOT=300`, increasing it for final figures. Each notebook runs independently. Tables and settings/source-file metadata save under `/srv/scratch/z5297792/SEACOFS_26yr_eddy_dataset_modular/eddy_env_tilt_controls/{eddy_properties,environment}/`. Each run creates a new UTC-stamped subfolder, avoiding stale result tables when a cache section is unavailable; figures display in the notebook. No results are committed as real-data findings.

## Definitions and limits

- Fitted `w` is relative vorticity; water vertical velocity uses explicit `vv_*` names. Signed rotation is converted to magnitude with AE/CE retained separately. Rossby number is `|zeta/f|`.
- `|Omega| Rc` is a speed scale in m/s, **not** a reconstructed maximum rotational speed. Axis ratio derives from eigenvalues of the positive-definite ellipse matrix.
- Depth targets 0, 100, 200, 500, 1000 m map to existing levels within 65 m; the resolution table shows actual levels/rejected targets. Profiles are capped at 1000 m. No new vertical interpolation occurs. Available/matched comparisons are paired by Eddy AND Day. Common-depth maxima also require every intervening cached level and every selected fitted property.
- Full-profile maxima are per eddy-day and their depth range varies. Thickness-weighted means require finite values throughout the retained sampled profile; isolated missing levels are not bridged. Maxima use available finite levels. Peak-depth ties select the shallowest.
- Profile shape is daily centreline geometry, not a replacement tilt calculation. Zero-chord profiles have undefined straightness/bending diagnostics. Lifecycle uses first/last observed Day, not the processed `Age` duration field.
- N2 extrema are core means of each water column's vertical maximum. Vertical velocity column means/RMS are weighted by sampled valid-cell count, not thickness/volume. The upper/deep velocity difference is not a measured stretching-budget term.
- Beta and latitude are not independent forcings. Surface-minus-0–500 m background velocity difference is a bulk shear proxy (m/s), not a derivative (s⁻¹).
- Binned curves are daily descriptions: long tracks can influence them more. IQRs show spread; bootstrap bars are pointwise uncertainty for the median, not simultaneous bands. Separate per-eddy-median correlations give an equal-eddy, between-eddy sensitivity. Adjusted models use inverse within-eddy row counts and clustered errors; they do not remove dependence between nearby distinct eddies.
- Models are separate per polarity/focal predictor, on that predictor's complete cases. Different coefficients may describe different samples. The focal predictor is not repeated as a control, and beta models omit latitude. Singular fits and insufficient samples are reported, not interpreted. Q-values are exploratory BH false-discovery adjustments within each model table.
- Normalising tilt by radius introduces mathematical coupling with radius-related predictors. Profile geometry and extent also share measurements with tilt. Interpret these as associations, not independent causal tests.

## Local validation

`MPLBACKEND=Agg python -m unittest discover -s seacofs_eddy_tilt_analysis/eddy_env_tilt_controls -p 'test_*.py' -v`

Synthetic checks cover signed extrema, ellipse geometry, depth-weighted means, matched samples, duplicated keys, cache coverage, bin boundaries, observed lifecycle and clustered model recovery. Full scientific execution requires Katana data.
