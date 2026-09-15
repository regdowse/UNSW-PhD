# Shear-relative tilt climatology

Population geometry and conditional evolution for **AE/CE x planetary/topographic
eddy-days**. This extends `shear_lift_tilt` with angle-conditioned evolution,
separate tilt/shear rotation, persistence matrices, explicit continuous episodes,
and day-weighted estimates with whole-eddy bootstrap uncertainty. It contains no
case studies and does not require a full Ertel-PV calculation.

## Run on Katana

After `git pull origin main`, open this folder and run:

1. `00_prepare_inputs.ipynb`
2. `01_shear_tilt_climatology.ipynb`

The existing `beta_effect_background_flow/01_build_background_cache.ipynb`
must already have produced the **v4** `eddy_day_background.parquet`. Input
paths and all analysis settings are editable near the top of the notebooks.
Uses the same Python environment as the existing analysis (numpy, pandas,
matplotlib, Parquet support and dependencies of `seacofs_tilt_tools`).

Notebook 00 loads measured tilt and background flow, calculates the current
nonlinear ESP-Gaussian environmental PV gradient inside **FRAC=1**, and stores
an input table. This initial all-eddy PV calculation can be expensive. A complete
cache is reused only when source fingerprints, PV options and helper hash match.
An interrupted initial calculation must be rerun; it does not write partial
science outputs. No new model velocity or temperature extraction is performed.

Cache and results live outside the repository:

```text
/srv/scratch/z5297792/SEACOFS_26yr_eddy_dataset_modular/shear_tilt_climatology/
  climatology_inputs.parquet
  input_metadata.json
  results/<UTC run timestamp>/
    settings.json
    *.csv
    *.png
```

Notebook 01 exports sample/episode counts, orientation distributions and
scorecards, lagged conditional curves, transition probabilities, PV-relative
tilt changes, latitude-stratified geometry and sensitivity tables. Saved figures
are generated from the real Katana inputs; repository notebooks have no outputs.

## Definitions and primary choices

- **Tilt:** authoritative coherent lower-to-upper `TiltDir`/`TiltDis`.
- **Flow contrast:** centred 91-day moving seasonal climatology, 0-200 m minus
  200-500 m. The existing climatology is a day-weighted combination of monthly
  climatological fields, not a daily rolling average of the event year.
- Recover the lower layer as `(500*U_0_500 - 200*U_0_200)/300`; require the
  sampled centre water column to reach 500 m. `water_depth_m` is centre
  bathymetry, whereas `h` is the Gaussian-weighted core mean used for regimes.
- **Geometry:** positive parallel means downshear; negative means upshear;
  positive transverse means left of shear. Angles are counterclockwise from
  shear, not compass-bearing subtraction. No AE/CE sign is imposed.
- **Regimes:** trailing seven consecutive-day median of log(topographic /
  planetary gradient magnitude), at least five finite values. Planetary <=
  `-log(2)` and `h >= 3000 m`; topographic >= `log(2)`. Missing current ratios
  stay unknown. Mixed and unknown observations remain available for audits.
- **Quality:** tilt >= 5 km, velocity contrast > 0.005 m/s. Both endpoints and
  every intervening day must qualify for a temporal pair.
- **Horizons:** exact 1, 3 and 5 days, primary 3 days. Regime changes, missing
  days, invalid tilt/shear or shallow columns break episodes. Eddy identity
  includes polarity. A track may enter the same regime multiple times.
- **Uncertainty:** 500 whole-eddy bootstrap replicates; suppress confidence
  intervals/conditional curves below 20 eddies or 50 observations per group.
  Point estimates and counts remain in tables. CIs are pointwise, not a
  multiplicity-adjusted declaration of significance.
- **Weights:** qualifying eddy-days are primary. Equal total weight per eddy
  within each reported population/bin is a separate sensitivity, never an
  implicit median over disconnected episodes.

Trailing rather than centred regime smoothing is intentional for temporal
analysis: future PV cannot determine today's label. Centred smoothing is a
predefined sensitivity matching earlier workflows. All pairs are rebuilt when
thresholds or regime definitions change.

## Interpretation limits

Whole-column coherent tilt and the layer-mean velocity contrast have different
vertical supports. Their comparison is a geometric diagnostic, **not a matched
two-layer motion budget**. Future work requiring that budget needs matched
depth definitions and separately validated centre estimates; these notebooks do
not silently replace your measured tilt with raw centre displacements.

Tilt-vector increments use fixed starting shear axes. Angular diagnostics
separate rotation of tilt from rotation of the background reference direction.
Wrapped endpoint angles cannot detect whole turns between samples. A diagonal
transition matrix indicates persistence, not proof of wave phase locking.

Conditioning on a noisy initial vector can produce apparent relaxation through
regression to the mean; requiring both endpoints above the tilt threshold also
selects trajectories. Multi-lag, threshold, weighting and background-definition
checks are provided, but they do not eliminate this issue. No residual is labelled
an internal-PV force. A preferred tilt does not prove baroclinic instability.

The PV proxy retains the existing topographic/planetary interpretation and
reconstructed eddy vorticity contribution to the topographic term. It is not a
pure external field or full 3-D PV. The optional directional comparison excludes
coherence < 0.2 and zero gradient magnitude, without assuming that tilt changes
must align with the gradient.

Latitude-stratified plots and sample counts show potential compositional
differences. A large climatology does not remove eddy interactions, systematic
background bias or cross-eddy dependence. Whole-eddy bootstrap accounts for
within-track dependence only. Regional/year blocks and explicit interaction
filters would be further publication-level robustness checks.

## Local verification

```bash
python -m unittest discover -s seacofs_eddy_tilt_analysis/shear_tilt_climatology -p 'test_*.py' -v
```

Tests cover geographic sign conventions, rotating reference axes, wrapped
angles, weighting, missing records, invalid intervals, regime re-entry,
shallow columns and duplicate keys. Notebook execution is also checked with
synthetic inputs locally; physical conclusions require the real Katana run.
