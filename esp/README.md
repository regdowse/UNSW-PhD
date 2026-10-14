# ESP manuscript revision notebooks

Open `01_noise_and_sampling.ipynb` on Katana, **restart the kernel**, and run
in order. `ESP_ROOT` defaults to `/home/z5297792/ESP_zonodo`. This directory
is never edited. Dependencies: NumPy, pandas, SciPy, Matplotlib and Jupyter.
No SEACOFS production data are needed.

## Controls and outputs

- `SAVE_OUTPUTS=False` (default): display tables/figures without creating a result
  directory or writing CSV/figure/provenance files. Jupyter can still save its
  own notebook outputs; this switch does not disable notebook autosaving.
- `SAVE_OUTPUTS=True`: write a fresh timestamped run below
  `results/01_noise_and_sampling/` (ignored by Git).
- `QUICK_RUN=True`: 20 repetitions; `False`: 100. Zero-noise controls run once.
- `SOLO_CENTRE_GUESS_M=10000`, `SOLO_SEARCH_HALF_WIDTH_M=60000`: explicit
  approximate along-track position and fixed search interval for the local
  SOLO variant. They are not derived from synthetic truth and stay fixed while
  varying the inner fitting window.

Coordinates/radii are metres, velocity m/s, rotation/vorticity s^-1. The noise
reference is the peak absolute *normalised* tangential speed
`abs(Omega)*Rc/sqrt(2*e)`. Noise SD applies separately to u and v. Spatial
covariance is `sigma² exp(-distance²/(2 L²))` within each component; u/v error
fields are independent. Scales are stress-test choices, not instrument estimates.

## Experiments and local variants

1. Fixed synthetic sampling with independent versus correlated errors at
   0/2/5/10/20% noise.
2. Inner-window sensitivity at 5% noise with a fixed outer observation domain.

Original inner estimators are loaded from the exact external `functions.py`.
DOPPIO and LATTE retain those functions. `local_estimators.py` provides:

- `solo_local_initialisation`: find a supported local sign change nearest the
  supplied approximate centre; fit an odd cubic locally within a fixed search
  interval; refit inside the inner window with bounded crossing and scaled
  coordinates. Same SOLO parameter equations. This changes the initial fitting
  domain, constraints and numerical scaling, not only the initial guess. No
  true eddy parameters are passed to this estimator. Its performance is
  conditional on the supplied approximate location/search interval.
- `outer_fit_diagnostic`: same Gaussian objective, sign filter, peak-based seed,
  optimiser bounds and original post-fit 100 km radius limit, but explicit
  optimisation/radius/insufficient-data/nonfinite statuses instead of silently
  returning seeds. Attempted estimates survive radius-limit rejection.

All methods use that local diagnostic outer fit in the updated notebook.
A baseline check compares it to the original outer function. The original and
local SOLO variants receive identical random observations, paired by experiment,
noise type, level and repetition. Different methods have different geometry and
sample counts; this is not a controlled algorithm ranking. The original executed
notebook remains in `01_original_backend_results.ipynb`; its helper references
are historical and outputs should be read as a saved baseline, not regenerated
with changed code. The exact previous implementation is available in Git.

Figures show median and central 90% simulation spread among valid fits and
valid-fit fraction over all attempts. Shading is not a confidence interval.
No accuracy cutoff removes large finite errors. SOLO has no shape/orientation
score. Orientation is axial (180-degree symmetry). Independent errors at
multiple levels share draws; window experiments use a separate random stream.
Settings selected from these sensitivity experiments need independent validation.

`validation_helpers.py` still supports the original outer function for checks;
its seed-return heuristic applies only in that compatibility mode. Original
inner APIs do not expose all convergence flags. 'Valid' means passing numerical
and physical diagnostics, not certified convergence or accuracy. Warnings and
failure reasons are retained. The full backend hash, local helper hashes,
settings and package versions accompany saved runs.

## Verification

From this directory in a scientific Python environment:

```bash
python -m unittest discover -p 'test_*.py' -v
```

Checks cover geometry/units, curl/divergence, covariance, angle wrapping,
failed-fit denominators, local crossing handling, explicit outer statuses,
and identical noise for original/local SOLO. Local reference-backend execution
is separate from scientific validation with the installed Katana backend.

## Notebook 02: background flow and radial-profile departures

Open `02_background_and_model_departure.ipynb`. It reuses notebook 1's fixed
sampling, local SOLO option, original DOPPIO/LATTE inner estimators, and explicit
outer-fit diagnostics. No existing backend or notebook-1 code is changed.

Three one-factor-at-a-time, noise-free tests: signed uniform currents in two
directions; signed simple shear in two orientations; and a heavier-tailed radial
profile family approaching the Gaussian. Defaults comprise 96 deterministic
conditions, not a Monte Carlo ensemble. There is no background correction.

`SAVE_OUTPUTS=False` displays only; `True` saves three figures, condition tables,
a baseline table and provenance under `results/02_background_and_model_departure/`.
The current and shear scales, fixed windows and SOLO initial guess are editable
in the first code cell. Radius-limit failures preserve attempted fitted radii.

Parameter errors target the underlying eddy component. Vector RMSE on a separate
fixed evaluation grid is reported against both the eddy and the total flow. The
profile family preserves centre vorticity and peak normalised-speed radius, not
peak speed. Compare peak-contour equal-area radius rather than inventing a true
Gaussian Rc for a non-Gaussian vortex. This tests radial misspecification only;
non-elliptical asymmetry, radial flow and temporal evolution are not included.

`departure_helpers.py` contains generation, scoring and plotting. Its five tests
check the Gaussian limit, peak radius and centre gradient, nondivergence, shear
curl and orientation, zero-background controls, and exclusion of fitting locations
from the evaluation grid. Run them with `python -m unittest test_departure_helpers -v`.

### Stage-specific diagnostics (notebook 02)

The user-selected `RC_LIMIT_M=1_000_000` m is retained. No comparison of radius
limits is run. `inner_valid` records a completed finite, positive-definite,
unit-determinant inner fit; `valid` records full-fit acceptance. Neither flag
is an accuracy criterion or a complete optimiser-convergence certificate.
Centre/geometry scores remain available after an outer-fit rejection, whereas
refined vorticity, radius and velocity scores require the accepted outer fit.
Notebook 02 plots use the appropriate flag for each metric. Existing notebook 01
ensemble summaries continue to use full-fit validity, as before.


## Notebook 03: ADCP–satellite platform influence

Open `03_platform_weighting.ipynb`, restart the kernel and run in order. It reads
`ESP_ROOT/data/sadcp_data.mat` and `satellite_data.nc`; paths, the ADCP window,
37 m target depth, fitting radii, two satellite weights and `SAVE_OUTPUTS` are
at the top. Requires netCDF4 (or h5py for the supplied NetCDF4 file) in addition
to the usual scientific packages. No original backend files are edited.

Compare ADCP-only, satellite-only, original equal-per-observation combined, and
selected weighted combined fits. Defaults are satellite inner weight 0 and
outer weight 1: ADCP determines inner geometry, then both platforms have equal
mean-error terms in the outer fit. The outer fit refines rotation amplitude as
well as radius, so final vorticity can still respond to satellite influence.
Sweep each weight separately over 0, 0.1, 0.3, 1 and 3, holding the other fixed.

Each weighted stage minimises ADCP mean squared residual plus lambda times
satellite mean squared residual. Inner residuals are linear-core velocity
vectors; outer residuals are normalised tangential Gaussian-profile errors.
Counts are stage-specific and outer counts follow the original sign filter.
Equal point weights correspond to lambda = N_satellite / N_ADCP, not lambda=1.
The weights specify influence, not observation-error variances. All available
points inside the selection domains are retained without thinning, but this
does not correct within-platform clustering or account for covariance.

The three original fits use their own original adaptive selections. Weighted
fits freeze the original combined core/outer observation sets, isolating weight
sensitivity from moving selection boundaries. The original outer sign filter
still depends on the fitted inner geometry; retained counts are reported.
`platform_weighting.py` supplies the local weighted affine solution of the
same five-parameter LATTE inner model, and an explicit-status weighted outer
fit. Checks compare both local formulations with the imported backend before
interpreting results. External original inner functions remain in use for the
three reference fits. Reference-data checks reproduce the released example.

This is one original 4.5-day ADCP-window comparison, not a time series. Preserve
native satellite -> daily 06 UTC -> ADCP mean-time interpolation. Actual native
satellite timestamps are recorded. The inner selection default is 35 km from
the original code, while the manuscript describes 30 km; reconcile this in the
revision. The outer selection remains 100 km; the configurable numerical Rc
limit remains 1,000 km. No radius-limit comparison is run.

Tables/figures show parameters, stage-specific sample support, common-set
per-platform residuals, and one-factor weight sweeps. Residuals overlap training
data and are descriptive, not independent validation. Surface geostrophic and
37 m total velocities remain physically different; temporal evolution, depth
shear and ageostrophic flow are not corrected. Uncertainties/effective sample
sizes are not estimated, and no weight is selected automatically as best.

`SAVE_OUTPUTS=True` writes a new timestamped directory under
`results/03_platform_weighting/`: parameter/residual/support/observation and
membership CSVs, three figures, and provenance including data/backend/local
hashes, settings, checks, fixed masks and evaluation IDs. False displays only.
The notebook is committed with cleared outputs for execution on Katana.

Nine focused tests check exact inner recovery, objective normalisation,
replication invariance, zero-weight exclusion, equal-point outer equivalence,
degenerate geometry, interpolation boundaries and explicit outer failures.
Run `python -m unittest test_platform_weighting -v` or the full suite above.


## Notebook 04: drifter–satellite complementary sampling

Open `04_drifter_satellite_weighting.ipynb`, restart the kernel and run all cells.
Defaults select 29 October 2023 at 06 UTC, the day illustrated at index 10 in
`LATTE_app2_multiplatform_data.ipynb`. Read `ESP_ROOT/data/drifter_data.nc` and
`satellite_data.nc`, using netCDF4. Paths, date, inclusive +/-24-hour window,
fitting radii, satellite weights and `SAVE_OUTPUTS` are at the top.

The notebook replays the original combined daily trajectory through that day
because selection depends on the previous day's centre/ellipse. Both inner
passes retain the previous day's Q, matching the original code. A satellite-only
trajectory follows that same rule independently. For a controlled inner
comparison, the drifter-only reference uses the selected 06 UTC window and the
fixed combined-core drifter subset, rather than the original standalone drifter
series's 12 UTC windows and unrestricted spatial selection. At the default day,
18 drifter observations are available and 9 enter that fixed inner subset.

Compare drifter-only inner, satellite-only full, original combined and weighted
combined fits. Defaults exclude satellite from the inner fit and give equal
platform mean-error terms outside. Sweep each satellite weight over 0, 0.1, 0.3
and 1 with the other fixed. Weighted selections remain fixed; the outer sign
filter remains geometry-dependent. The outer-zero scenario is diagnostic only:
a finite radius does not establish outer support. Unsupported fitted speed peaks
are marked in the radius plots. No automatic selection of an optimal weight.

The sampling map and radial distributions use the original combined ellipse
for both platforms. The fixed-core comparison and full Gaussian comparison
have separate residual tables; both overlap fitting data and are descriptive.
Sample counts, active objective mass, time/depth support and number of distinct
drifters are reported per fit/stage. Depths are nominal (drifter 15 m, satellite
surface). Error variances, covariance and effective sample size are not inferred;
ageostrophic flow, vertical shear and temporal evolution remain limitations.

`drifter_weighting.py` supplies original-example loading/selection, comparison
and plots. It reuses the notebook-3 inner/outer solvers with an explicit
`reference_platform='Drifter'`; the default remains ADCP so notebook 3's behaviour
is unchanged. All original inner reference fits call the external backend;
weighted variants remain local. Formulation checks run against that backend.

`SAVE_OUTPUTS=True` writes fresh CSV tables, membership, three figures and
provenance under `results/04_drifter_satellite_weighting/`. False writes no result
files. Notebook outputs are cleared for the first Katana run. Eight focused tests
cover reference-label equivalence, fixed masks, outer/inner separation, missing
drifters, inner-only reporting, time boundaries, lagged-Q selection and radial
support. Local reference-data execution reproduced the released illustrated
centre and radius within 0.1 m; this does not certify scientific accuracy or a
changed Katana backend. Run `python -m unittest test_drifter_weighting -v`.


## Notebook 05: drifter window sensitivity

Open `05_drifter_window_sensitivity.ipynb`, restart and run all cells. It reports
**days with an inner-core result out of common comparison days**, parameter time
series, and matched-date differences from the current 48-hour window. Total
window durations default to 12, 24, 48, 72 and 96 hours (48 means +/-24 hours).
`SAVE_OUTPUTS=False` displays only. True writes three CSVs, one figure and
provenance under `results/05_drifter_window_sensitivity/`.

Compare drifter only, original equal-per-observation combined, and combined
platform-normalised weighting with fixed satellite multiplier 0.1. No weight
sweep or outer fit. Original drifter-only fitting uses all temporal-window
observations; combined fitting uses the original 35 km two-pass spatial selection
with lagged Q. Cross-configuration differences therefore include spatial-selection
effects. Within each configuration, only temporal window length changes.

Use the same daily 06 UTC dates for every case, trimming the largest half-window
at both ends of the available drifter archive. Defaults yield 96 dates from
22 October 2023 to 25 January 2024. This is an interior-period comparison, not
coverage of the trimmed archive endpoints. Each combined trajectory starts
at the same supplied approximate centre and identity Q, then uses its own
previous fit. It is not a replay of the earlier original full-period trajectory.

Stop a run at its first rejected inner fit, without restart or borrowing another
track's centre. Subsequent rows are explicitly not attempted, rather than
additional failed fits. There is no longest-track or percentage metric. Counts
require finite valid elliptical geometry, full design rank and expected rotation
sign; combined runs retain the original 100 km between-pass refinement limit.
No accuracy or new inter-day-jump cutoff is imposed. Counts do not certify eddy
identity or parameter accuracy; large finite excursions remain visible.

`drifter_windows.py` reuses the original backend and notebook-4 data preparation.
Five focused tests verify common dates/total-duration semantics, stopping without
restart, full count denominators, matched-date differences, missing references
and rank rejection. Run `python -m unittest test_drifter_windows -v`. Existing
backend/helper files and saved notebook outputs are not modified.
