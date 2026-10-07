# ESP manuscript revision: notebook 01

Open `01_noise_and_sampling.ipynb` on Katana and run the cells in order.
Set `ESP_ROOT` to `/home/z5297792/ESP_zonodo` (the default). No production
SEACOFS data or tilt caches are needed. Requires NumPy, pandas, SciPy,
Matplotlib and a Jupyter Python kernel.

The notebook imports `solo`, `doppio`, `latte` and `out_core_param_fit` from
the exact external `functions.py`. It never edits that checkout. Shared
experiment/diagnostic code is in `validation_helpers.py`; no corrected
replacement ESP estimators are used in this version.

## Scope

1. Noise sensitivity: fixed synthetic Gaussian vortices and sampling, independent
   versus Gaussian spatially correlated velocity errors, 0/2/5/10/20% noise.
2. Inner-window sensitivity: fixed 5% noise and outer sampling, varying the
   inner window relative to the known synthetic Rc.

SOLO uses a circular vortex; DOPPIO/LATTE use the same elliptical vortex.
These are representative configurations, **not** exact reproductions of the
paper's tables, nor a matched-observation ranking of the three methods.
The noise-free baseline is a self-consistency check, not real-world validation.
No background flow, time evolution, platform fusion or depth/tilt analysis is
included. Default noise and correlation scales are stress-test choices, not
instrument-error estimates.

## Running and outputs

`QUICK_RUN=True` uses 20 realisations per condition for a first inspection;
`False` (default) uses 100. Noise-free controls are evaluated only once.
Errors at different levels share random realisations to make comparisons
paired. The window experiment uses a separate random stream. Geometry is
fixed across all repetitions. Window settings are exploratory sensitivity
tests; selected settings would need an independent validation ensemble.

All coordinates and radii passed to ESP are in **metres**, velocities in m/s,
and Omega/vorticity in s^-1. The SEACOFS workflow's kilometre convention is
not used here. The noise reference is the peak absolute **normalised**
tangential speed, `abs(Omega)*Rc/sqrt(2*e)`, not the maximum physical speed
over an elliptical contour. Each component has the stated noise SD.

Outputs go to a new timestamped `results/01_noise_and_sampling/...` directory:
raw attempts, summaries, failure counts, settings/backend provenance, and PNG/PDF
figures. Outputs are ignored by Git. Both stages are refitted for every attempt.
Outer input locations stay fixed; the backend's rotation-sign filtering can
still change the retained outer observations. Their count is saved.

The original backend may silently return initial values when its outer fit
fails or exceeds `Rc_max`. The wrapper passes the same peak-based initial
guess explicitly and flags an unchanged parameter pair as `outer_seed_return`.
This is a conservative heuristic, **not** a guaranteed optimiser-status API;
an unchanged genuine optimum can be flagged. The default 100 km `Rc_max`
is retained and recorded. Original inner optimiser success flags are not
exposed by these APIs; 'valid' means passes the listed numerical/physical
diagnostics, not certified convergence or scientific accuracy. Warnings are
retained in the raw results. Negative-determinant/non-positive-definite Q,
wrong rotation sign, insufficient observations and nonfinite results are flagged.
No accuracy cutoff is used to remove large but finite parameter errors.

Figures show medians and the 5th--95th percentiles **among valid fits**, alongside
valid-fit percentages using all attempts. Shading is simulation spread, not a
confidence interval. SOLO's imposed shape has no shape/orientation score.
Orientation is axial (180-degree symmetry); angular errors concern the chosen
noncircular synthetic truth and do not establish identifiability near circularity.

## Checks

From this directory, in a scientific Python environment:

```bash
python -m unittest -v test_validation_helpers
```

These checks cover geometry, units, curl/divergence, axial angles, noise
covariance, shared cross-transect observations, invalid fits, fallback handling
and summary denominators. They do not substitute for running the notebook
against the installed Katana backend. Inspect the baseline and failure tables
before drawing conclusions. Keep existing backend issues visible; proposed
estimator changes belong in separately named local functions and comparisons.
