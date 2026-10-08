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
