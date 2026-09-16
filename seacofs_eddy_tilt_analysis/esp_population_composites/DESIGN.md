# ESP population composites — proposed analysis

Original design, recorded 16 September 2026. The initial four-group workflow is
now implemented; see README.md for delivered scope and remaining extensions. See
`NOTEBOOK_INVENTORY.md` for the local repository notebook survey.

## Scientific target

Estimate the mean reconstructed horizontal velocity structure through depth for
specified eddy populations, preserving displacement of deep centres relative to
one fixed shallow reference. Initially compare AE/CE and planetary/topographic
populations; extend to Rossby number, stratification, region, shear and age.
ESP supplies horizontal velocities at multiple depths, not vertical velocity or
temperature, density or N² fields. Existing scalar N² caches can define groups;
hydrographic composites require separate spatial model data.

Reconstruct each eddy-day before averaging. A mean of ESP parameters is generally
not the mean velocity field. Keep the mean fitted centreline alongside the mean
field: the centre of a superposed velocity field need not equal the mean of its
constituent fitted centres. A broad or weak composite may reflect variable tilt
directions rather than weak individual eddies.

## Existing building blocks and limitations

- `case_studies/eddy_cross_sections_3d/esp_composite_tools.py` provides single-eddy
  daily reconstruction and optional radius scaling. Its Day-only joins and
  grouping must not be reused on pooled eddies; population keys are `(Eddy, Day)`
  and `(Eddy, Day, Depth)`.
- `depth_resolved_pv_tilt/` already uses matched snapshots for depth comparisons.
- `tilt_mechanisms/` supplies validated v4 potential-density stratification
  summaries; these describe eddy cores, not independent background stratification.
- `shear_tilt_climatology/` supplies environmental shear and clustered analysis.
- `surface_pv_esp_gaussian/` and `surface_pv_footprint_sensitivity/` distinguish
  net gradient direction from strong but cancelling local topographic exposure.
- Surface processing smooths parameters, fills missing days and renumbers IDs.
  Its final columns do not preserve interpolation provenance. Do not infer an
  observed fit merely from finite parameters or join raw IDs without a mapping.
- Vertical fitting already checks polarity, centre jumps and local support.
  The final confirmation stage only requires a minimum row count (default three);
  it does not supply depth-resolved reconstruction-error certification.

## Selection: quality first, scientific grouping second

Create an audited table with one row per eddy-day and explicit exclusion reasons.
Require unique keys, finite required parameters, positive radius, positive-definite
ellipse matrix, valid reference centre and the exact requested fitted depths.
Check polarity and units against the reconstruction implementation. Report source
and cache versions; validate joins rather than silently dropping duplicate rows.

Use one fixed shallow reference depth across the primary sample. Do not let the
reference silently migrate deeper as shallow measurements disappear. Retain only
days with all selected analysis depths for each primary composite. Publish a
shallower companion analysis to show the selection effect of requiring deep
profiles. Available-depth composites are secondary and must show support with
depth. Never substitute zeros for missing profiles or interpolate vertical fits.

Reconstruct provenance where possible using upstream processing plus a stable ID
mapping; otherwise label it unknown. A genuinely refitted subsurface profile on
a gap-filled surface day is distinct from an interpolated vertical profile.
Establish that distinction before choosing an observed-only sensitivity sample.

Audit reconstruction skill against background-removed SEACOFS velocity on a
stratified sample spanning polarity, region, depth, radius, Rossby number and
candidate quality ranges. Use joint vector RMSE/normalised error and support;
do not rely solely on correlation. Calibrate any additional fit-quality threshold
from these diagnostics before examining the desired tilt result. Universal skill
filtering requires computing that metric for every selected eddy-day.

Keep small and zero tilts in structure/magnitude composites. Do not select for
onshore tilt, predicted PV alignment, straight centrelines or a visually clean
eddy. These are outcomes. Large real bends, unusual Rossby numbers and weak
stratification are not automatically bad fits. Record retention by population.
Directional statistics may use a separate declared small-tilt sensitivity cutoff.

## Coordinate frames and physical support

- Geographic east/north: primary for planetary/geographic tilt questions.
- Slope frame: positive cross-slope axis points towards shallower water,
  `-grad(h)`; rotate positions AND velocities by the same rigid transformation.
  Use one orientation for the entire eddy-day column, not a different rotation
  at each depth. Compare a declared shallow-core orientation with a column-based
  orientation as a sensitivity.
- Shear frame: later comparisons relative to a declared upper-minus-lower shear.
- Tilt-aligned frame: only for conditional shape, never evidence of preferred tilt.

A poorly defined slope direction disqualifies an observation from slope-aligned
analysis, not from every composite. Quantify slope strength, spatial directional
coherence and valid seabed support. Seamounts with cancelling gradients are a
separate population; do not assign them a confident onshore axis from noise.
Keep geographic versions of the same selected samples for interpretation.

Recenter each complete column once. Recentring each depth separately erases tilt.
Use the shallow reference Rc for radius-normalised horizontal coordinates at all
depths; depth-dependent scaling can change the apparent vertical structure.
Compare radius-scaled composites with kilometre-coordinate composites. Retain
dimensional velocities first; velocity normalisation answers a different question
and should be a separately labelled sensitivity. Maintain positive-down exact
physical depths initially.

ESP can mathematically extend over land or below local bathymetry. Distinguish
idealised parameterised fields from ocean-supported fields. For the latter,
transform and apply each day's wet/bathymetric mask and report spatial support;
never interpret extrapolated vectors as observed water-column velocities.

## Population definitions and weighting

Define groups at eddy-day level, so an eddy can contribute eligible days to more
than one regime. Retain mixed PV conditions. Use the existing 2:1 ratio as an
initial regime definition, with sensitivity to thresholds and surface versus
depth-averaged environmental PV. Keep the eddy-internal PV gradient diagnostic
separate. Topographic dominance is not proof of dynamic topographic control.

Use |Ro| for intensity bins while preserving signed Ro and AE/CE labels. Declare
whether Ro is shallow or column-averaged. Show actual bin boundaries and both
common physical boundaries and within-polarity quantiles when useful: they answer
different questions. Do not build all possible crossed subgroups initially.

For stratification, start with v4 core N² and its existing 80% valid-core-coverage
criterion; label the averaging depth. Treat environmental N as a future distinct
input. Compare raw and geographically/seasonally matched groups, auditing radius,
depth coverage and regime balance. Match only covariates relevant to the question;
do not automatically control away the defining mechanism.

Primary mean: average eligible daily fields within each eddy and population,
then average the resulting eddy means equally. Day-weighted means describe a
random eligible eddy-day and are a useful sensitivity. Report both unique eddies
and days, per-eddy contributions, and support at every depth/grid cell.

Bootstrap whole eddy IDs, keeping all their eligible days together. When comparing
groups, resample the global ID list jointly so shared eddies remain coupled.
Consider spatial/year blocks as a sensitivity to dependence between eddies.
Publish centreline uncertainty and distributions, not only smooth mean arrows.
Do not treat the number of eddy-days as the independent sample size.

## Proposed notebook sequence

1. `00_population_audit.ipynb`: join/provenance checks, fit/coverage diagnostics,
   selection attrition, exact-depth catalogue, frozen group definitions.
2. `01_polarity_regime_composites.ipynb`: AE/CE by planetary/topographic regime,
   appropriate reference frames, matched depth samples and support maps.
3. `02_rossby_stratification_composites.ipynb`: one contrast at a time, common
   support, raw/matched comparisons and explicit physical bin boundaries.
4. `03_uncertainty_sensitivity.ipynb`: eddy bootstrap, QC/depth/normalisation/
   weighting sensitivities, directional spread and representative member cases.

For scale, stream daily reconstructions and cache per-eddy/group sufficient
statistics and compact centrelines, rather than retaining every daily 3-D field.
Store selection manifests, configuration hashes, input versions and exclusion
counts. Freeze selection before reconstruction. Synthetic verification must
cover rotated vectors, retained tilt, colliding Day values across eddies,
unequal track lengths, missing depths, masks and shared-ID bootstrap contrasts.

First deliverable should be four well-audited composites and their uncertainty,
not a large gallery of potentially incomparable combinations. Thresholds beyond
existing conventions remain provisional until the Katana data audit is run.
