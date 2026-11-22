# Manuscript workflow review / agent handoff

Analysis decisions and historical progress; moved here on 2026-11-17.
Current revision work is tracked in [revision/checklist.md](../revision/checklist.md).
Use the latest topic reports for numerical claims; dated results below may be historical.

Review date: 2026-11-16 (Australia/Sydney). Local source: UNSW-PhD `main`,
HEAD `bb1d6c6` (figure tags), with user's uncommitted `manuscript/info.md` update.
This document records the investigation and decisions. Historical planning-only
statements below are superseded by the current implementation status here.

## Dataset refresh (2026-11-17)

User requested clearing and rerunning all seven manuscript notebooks after the
source dataset update. Current source: 133,038 eddy-days / 3,097 eddies;
106,865 valid tilts and 97,012 directional days. Earlier numerical summaries
below describe the previous dataset; use the current topic outputs/report.md
and CSVs for manuscript revision.

The rebuilt PV cache passes current input/settings/code/content checks and was
reused. Background monthly fields, grid, sigma depths, calculation kernels and
archive identities are unchanged; only the source eddy-day table differed.
Validated all existing product hashes, retained the monthly fields/partitions,
and refreshed the three background sampling tables for 133,038 rows each.
The previous manifest is retained in scratch, with the sampling refresh recorded
in the current manifest. No archive pass or PV calculation was repeated.

Example IDs changed with the dataset: Figure 2 CE 613/3370 -> 643/3370;
AE 76/1675 -> 81/1675. These are identified from matching surface-vorticity,
footprint counts (420/474), gradient direction/magnitude and coherence fingerprints,
not by assuming a surviving ID represents the same eddy. Figure 3
2103/7923 -> 2191/7923 is identified from the retained reference profile's centre
(the next nearest current surface centre is >43 km away). Configs now pin these
current IDs. Figure 3 recomputes to 3.593093 km / 66.566643 degrees, matching saved
tilt within numerical precision; current projected extent is 14.498783 km.
These updated values supersede the former 3.7/21.1 km labels.

All notebooks explicitly enable inline Matplotlib, clear previous outputs and
run in fresh myenv310 kernels. Executed notebooks retain figures and report/table
outputs as requested. Figures, CSV/Parquet results, reports and manifests are
regenerated in each topic's outputs directory. All 49 scientific tests pass.
No changes to analysis definitions or plotting styles; no commit or push.

## Tilt-method illustration status (2026-11-16)

`tilt_method/` implements Figure 3 with JSON, a short notebook and local delta
calculation/plotting helpers. Fixed Eddy 2103 / Day 7923 (window 7921–7925), no
random example search. No legacy or upstream runtime imports, PV/background
product dependency or scratch writes. This completes requested Figures 2–17;
Figure 1 and the extra post-composite figure remain excluded as agreed.

Production-equivalent settings: centred five calendar days, Gaussian sigma 1 day,
10 m increments, 1000 m cap, minimum 200 m fit-label range / five fit points,
inverse unweighted sample-variance depth weights with epsilon 1e-10, 20-degree
bearing offset. Interpolate only inside each day's capped observed range; retain
reference-day-supported intervals before averaging, variance and reconstruction.
Keep existing upper-interval labels and mixed numerical units (horizontal km,
depth m) in the weighted 3-D SVD. This is an illustration of the current estimator,
not a change to production tilt. The no-interpolation control/composite decisions
do not change this estimator's own interpolation step.

Fresh result: tilt 3.6888208011596264 km, 64.9503132314173 degrees; saved distance
differs by 4.44e-16 km, circular bearing difference zero. Maximum projected extent
21.070722511684558 km. All five days have 22 capped fitted levels; reference
support 0–858.918478 m, fit labels 0–840 m (85 points, five contributions each).
The sub-5-km methods example is independent of the population directional cutoff.

Original three-panel style retained and visually compared with the tagged saved
image. The ribbon is relative fit weight, not uncertainty. Extent is max-minus-min
reference position along the fitted axis, distinct from maximum pairwise 2-D or
surface-relative distance. Display positive axis is shallow-to-deep; arrow and
saved bearing are deep-to-shallow. Outputs include PNG/PDF, concise writing report,
summary, source/support and step tables, manifest with input/config/code hashes.
Eight new scientific checks bring the full suite to 49 passing tests.
Fresh-kernel execution passed. Full-example checks verify 134 source / 110 capped
profile points, 425 increments and 85 reconstruction points; temporal means,
sample variances, fit weights, endpoint distance, projected extent and manifest
hashes. A separate read-only check agrees with the current production functions
within 1e-12. The stored notebook has no saved output clutter. No source dataset
was regenerated; no commit or push performed.

## Background implementation status (2026-11-16)

User authorized Figure 11's agreed implementation. `background/` now contains
JSON settings, a short notebook, local climatology/layer/bootstrap kernels and
the tagged three-panel plotting style, including `seaborn-v0_8-whitegrid`.
No legacy runtime imports. Independent scratch root:
`/srv/scratch/z5297792/eddy_tilt_manuscript/background_flow/monthly_v1/`.
The full independent rebuild completed: all 307 partitions were freshly computed,
then reduced into monthly means/counts and 31/61/91-day all-source tables. No
legacy cache or partition was used as an input. Fresh-kernel notebook execution
passed and the regenerated figure was compared with the tagged saved output.

Archive preflight found 307 files / 9,189 unique consecutive daily records,
1994-01-02 through 2019-02-28. This is 26 calendar years, with incomplete first
and final years, not 26 full years. All source eddy-days are within the archive.
Current deep-planetary directional population: 6,103 days; all centre cells are
wet and >=500 m deep (minimum 1423.99 m). Footprint-mean h remains the 3000 m
population selector; centre-column bathymetry checks layer coverage separately.

Preserve the original static physical sigma-centre / midpoint-interface layer
approximation. Read all 30 model levels; surface is the top model level. Require
complete finite 0–200 and 0–500 m coverage, rather than renormalising missing
depths. Derive 200–500 m from those layer integrals, and retain geographic
east/north fields. Monthly sums/counts use all archive days, including days
without eddies. The circular 365-day seasonal calendar maps Feb 29 to Feb 28.
31/61/91-day windows use identical monthly fields and complete matched days;
91 is primary. Pooled day medians and whole-track bootstrap (5000 draws, seed
42 with exported statistic offsets), reusing seeds across windows.

Fresh partitions are atomic and restartable only with matching input/kernel
signature. A completed product is validated before use; changed inputs/code
require a fresh product directory. Raw archive provenance records paths, sizes,
mtime_ns and all ocean-time values, not a claimed content hash of ~0.8 TB.
Small source files, kernels and generated products have SHA-256 identities.
Eight new checks bring the suite to 41 passing tests. A fresh one-file partition
and real-grid 200/500 m thickness checks passed before the full run.

Full-data coverage validation retains all 6,103 deep-planetary directional days:
AE 4,347 / 489 eddies; CE 1,756 / 378. No layer/window exclusions. Every targeted
cell has all expected archive days for every month, component and layer. Fresh
91-day values agree with the historical background table after keyed comparison:
surface exactly; layer means differ by at most 8.33e-17 m/s (round-off). This is
validation after independent regeneration, not reuse of the historical cache.

Matched internal checks preserve the zonal-contrast population result. Median
surface-minus-200–500 m zonal velocity (m/s): AE 31/61/91-day =
0.04813/0.04744/0.04719; CE = 0.04803/0.04774/0.04763. Every lower 95% whole-track
CI bound is positive. Primary 91-day CIs: AE [0.04474,0.04955], CE [0.04430,0.05049].
The median of paired daily differences is distinct from a difference of medians:
31 minus 91 days gives AE -0.00039, CE -0.00054 m/s; individual absolute changes
can reach ~0.042 m/s. Do not imply individual days are unchanged or that window
agreement establishes causation. Internal checks remain in tables, not paper panels.

Outputs: Figure 11 PNG/PDF, layer median/CIs, window contrast summaries and paired
differences, monthly coverage, matched members, selection/exclusion audit, writing
report and manifest. Numerical/bootstrap settings are in JSON, with actual seeds
in the CSVs. The original whitegrid style, 1x3 layout, colours, labels, maps and
98th-percentile shared colour saturation are retained. A small automatic y-margin
keeps the original panel-a label from colliding with the new top tick; no estimates
or CI are clipped. Figure 3 is now implemented in the topic above.
No commit or push performed.

Final checks passed: 18,309 matched window rows, all 36 layer summaries, six
contrast summaries and 144 monthly-support checks; daily layer algebra, pooled
medians, identical Figure 10 target membership, complete CIs and manifest hashes.
All 41 tests pass. The final display adjustment was re-exported and visually
checked; it changes only automatic axis headroom, not numerical tables.

## Composite implementation status (2026-11-16)

`composites/` contains JSON settings, a short notebook and local calculation/plotting
helpers for Figure 17. Reads confirmed fitted profiles, saved tilt and validated
fresh surface PV. Uses shared unsmoothed magnitude categories and 3000 m target
selection: deep planetary / all-depth magnitude-defined mixed / shallow topographic.
Other-depth dominated days remain excluded, never relabelled mixed. Shared directional
QC includes saved tilt >=5 km. No legacy imports, external ESP backend, interpolation,
velocity reconstruction or scratch rebuild is needed.

Retains valid actual fitted centres within 0–1000 m inclusive, referenced to each
day's own valid fitted surface centre. Rotate individual offsets into that day's
total environmental-PV frame before averaging; along is with the signed gradient,
perpendicular is +90 degrees counterclockwise. Displacements are surface-to-depth,
distinct from saved fitted deep-to-shallow TiltDis/TiltDir. Each available day has
equal weight at each exact depth. Pointwise 95% CIs use whole-track resampling
(500 draws, seed 731), with the same draw across depths within each group.

Retained days/eddies: planetary AE 4,347/489, CE 1,756/378; mixed AE 14,601/931,
CE 9,411/874; topographic AE 5,297/335, CE 8,258/586. Current deepest contributing
fit is 858.918478 m in every group; depthwise support is exported rather than
assuming constant membership. At that depth: planetary AE 3,925/414, CE 1,589/342;
mixed AE 12,118/818, CE 8,157/802; topographic AE 3,100/257, CE 6,090/523.
Eddy counts overlap across regimes. Surface-reference means and CIs are zero by
construction; a composite's vector magnitude is not the mean member magnitude.

Original six-panel layout, colours, confidence shading, labels and symmetric
row-wise axes retained and compared with the tagged output. Column titles clarify
deep planetary / shallow topographic selection. Local outputs contain PNG/PDF,
depthwise means/CIs/spread/support, inventory, selection/exclusion audit, projected
member centres, concise writing report and provenance manifest. Seven new tests
bring the suite at that stage to 33 passing checks. Background status is above;
Figure 3's method illustration is now implemented above.
Fresh-kernel notebook execution passed. Full-data checks covered 936,326 fitted
profile rows / 43,670 selected days and all 132 group-depth estimates: zero surface
offsets, distance-preserving rotation, exact pooled means/support, identical target
membership to Figures 10/12 and current manifest hashes. No commit or push performed.

## Controls implementation status (2026-11-16)

`controls/` has JSON settings, short notebook and local analysis/plotting helpers.
Figures 13/14 use fixed existing levels 0, 181.2647932696809 and
515.416489345845 m, with the same days/properties at all levels. No ranking or
interpolation. Original |Omega| and AR plots and surface-Rc normalisation retained.
Actual structure sample: AE 43,889 days / 1,344 eddies, CE 45,495 / 1,491.

User approved Figure 15's five exact levels 0, 48.89763802286193,
105.61553821400828, 181.2647932696809, 327.44085148324933 m. Depth comparison
matches days after the directional AR/tilt QC at every level. Surface histograms
and AR-class panels have their own documented samples. Day-weighted means and
1000 whole-eddy bootstrap resamples (seed 42) replace old equal-eddy estimates;
original six-panel style retained. Depth-matched sample: AE 38,460 days / 1,309
eddies, CE 39,549 / 1,483. Preserve the inherited support thresholds: >=5 days
per eddy in the actual sample/class and >=20 eddies per displayed estimate.

Figure 16 joins current surface latitude with validated fresh surface-fit PV.
Plots remain unadjusted medians/IQR and fits to bin medians; no causal inference.
User approved raw y-limit 50 km (old 40 clips IQR reaching 43.7 km), while the
normalised limit stays 0.6. An additional Figure 13 shared-axis clipping issue
(IQR reaches 48.1 km) was also approved for correction: the shared upper bound
now expands to show all IQR shading while preserving layout and styling.

Outputs include figures, bin estimates/counts, matched members, QC/eligibility
audits, ellipse bootstrap estimates/coverage/paired depth contrasts, descriptive
environment fit coefficients and report.md. No new scratch products or legacy
runtime imports. Tests are in the shared ignored tests directory. Controls are
not a dominance/depth-selected population. Composite status is recorded above;
background status is recorded above.
Seven controls tests bring the full suite to 26 passing checks. Fresh-kernel
notebook execution passed. All four figures were compared with saved tagged
outputs. Additional checks verify exact depth-key matching, histogram percentages,
matched ellipse support and readable exports. The final axis correction does not
change numerical CSVs; raw structure/environment IQR bands are checked for clipping.

## Surface-PV implementation details (2026-11-16)

Implemented independent `pv/pv.py`, `pv/pv_kernels.py`, `pv/pv_plots.py`, JSON
settings and a short notebook. Required scientific/plotting kernels are local;
no legacy runtime imports. Explicit fresh build completed at
`/srv/scratch/z5297792/eddy_tilt_manuscript/surface_pv/surface_fit_v1/`:
`surface_pv.parquet` and a provenance manifest. Existing caches were untouched.

127,426 rows, unique Eddy-Day keys; all have finite environmental PV gradients.
No cached tilt outcomes, temporal smoothing or depth-population filtering. Surface
geometry and surface w drive PV. Auxiliary signed/absolute peak w, selected depth
and support count are retained separately for future work. Hash-validated loads
fail on changed inputs, settings, calculation code or product content; rebuilds
must be explicit. Formula implementation is the original environmental gradient
with Gaussian-envelope weights, not full 3-D Ertel PV.

Validation: five new scientific tests (ten including census) pass. All 44 shared
PV covariates agree with the historical surface-fit comparison cache on all
127,426 matched days (rtol=1e-12, atol=0). Auxiliary peak amplitude/depth match
the old peak-w cache. This is a comparison after independent regeneration, not
use of either old cache as an input. Fixed Figure 2 examples reproduce
1.146674e-11 and 7.482916e-14 m^-2 s^-1. Visually compared the regenerated PNG
with tagged notebook cell 7: original layout, bathymetry, ellipse/profile lines,
magenta local arrows, green mean arrows, legends and colourbar retained.

Outputs: `pv/outputs/` contains Figure 2 PNG/PDF, examples.csv, summary.csv,
sample_audit.csv, concise report.md and run manifest. No inference/CI is attached
to this method illustration. Tests are in the shared Git-ignored tests directory.
The notebook passed fresh-kernel loading, product validation and figure export.
Figures 7–10/12 and regime classification are now implemented below. Figure 11
is implemented in the background topic above.

Figure 2 arrow-frame correction (approved by user): rotate both local magenta
vectors and the green mean vector from geographic east/north into model x/y
before drawing. For grid angle a in radians, x=cos(a)E+sin(a)N and
y=-sin(a)E+cos(a)N. The reference grid has a=-20 degrees, so display arrows turn
+20 degrees relative to the old plot. Magnitudes, colours, mean-arrow display
length, layout and underlying geographic products are unchanged. Tests cover
cardinal directions, inverse rotation, preserved magnitudes, actual Quiver
components and mean-arrow endpoints. Thirteen tests pass.

Scientific input/function hashes now determine product validity separately from
plot/report code. The manifest was upgraded after verifying unchanged numerical
functions, kernels, inputs and product hash; original provenance is retained in
its update record. No PV recomputation or Parquet rewrite was performed.
Regenerated Figure 2 PNG/PDF and writing report; visually checked the corrected
directions. Example/summary/audit CSVs and the PV Parquet are byte-for-byte
unchanged. The original figure style and arrow reference scales are retained.

## PV-dominance implementation status (2026-11-16)

`pv_dominance/` contains JSON settings, short notebook and analysis/plotting Python.
Reads validated fresh PV covariates, joins current tilt by unique Eddy-Day keys,
and reuses geography/rose helpers from manuscript census. No legacy runtime
imports, old-cache inputs or temporal smoothing. Dominance remains magnitude-only;
targeted populations additionally use the shared 3000 m water-depth threshold.
Factor-two categories include their boundaries; invalid, negative and zero/zero
inputs are unclassified. Magnitude categories are independent of tilt validity.
Directional comparisons use tilt >=5 km and valid total PV/tilt bearings.

Full-data classification: AE planetary/mixed/topographic days 5,662/20,259/39,031;
CE 2,388/13,131/46,955. No unclassified days in current inputs. Directional
planetary sample: AE 4,352 days / 494 eddies; CE 1,757 / 379. Directional
topographic: AE 27,814 / 1,282; CE 33,945 / 1,470. Eddy counts overlap by regime.

Figure styles ported from tagged cells for Figures 7–10/12. User approved three
specific corrections during this step: Figure 7 geographic bearings rotate by
the actual grid angle (-20 degrees) for model-grid display; Figure 9 limits
expand to cover unsmoothed tails (raw log ratio -5.9894 to 6.8938); Figure 10
retains AE radial minimum 270 and expands CE 170 to 200 to accommodate peak 188.
Figure 10 selects deep planetary days (h>=3000 m); Figure 12 selects shallow
topographic days (h<=3000 m). Here h is Gaussian-footprint-mean water depth from
PV, not centre-column or profile depth. Other-depth dominated days keep their
dominance category and are excluded from these targets, never relabelled mixed.
Figures 7–9 and the all-depth inventory remain unchanged. The same JSON threshold
drives population selection and the plotted isobaths. The contour uses local grid
bathymetry, so its position need not exactly match footprint-mean selection.
Figure 12 keeps the original longitude split at 154.75 E; equality is off-shelf.
It should not be described as a new bathymetric shelf classification.

`outputs/` contains PNG/PDF figures, concise writing report, regime proportions,
directional summaries (all and longitude groups), gradient quartiles, figure
sample/exclusion counts, spatial mean-log-ratio bins/support and row-level audit.
These are descriptive results, not independent-day significance tests. Five new
tests cover boundary/invalid/zero inputs, scale/component symmetry, instantaneous
classification without depth filters, direction thresholds, circular differences
and key coverage. A depth-boundary/exclusion test brings the total to nineteen
passing tests, including earlier topics. Scratch inputs
and PV products are read-only in this step. No commit/push performed.
Fresh-kernel notebook execution passed. Compared all five figures visually with
saved tagged outputs. Verified category percentages/counts, spatial-bin support,
unclipped planetary roses and exact run-manifest code hashes. Figure 7 excludes
two AE source days without a region label; neither had valid tilt. Figure 11
is implemented in the background topic; controls/composites are also implemented above.
The all-depth directional counts above describe the pre-depth-filter population;
current targeted counts and exclusions are in report.md, figure_samples.csv and
depth_selection_counts.csv. direction_summary.csv explicitly labels all-depth
versus deep_planetary/shallow_topographic populations; audit retains both regime
and target_population, including exclusion reasons and h.
After depth selection and directional QC: Figure 10 AE 4,347 days / 489 eddies,
CE 1,756 / 378; Figure 12 AE 5,297 / 335, CE 8,258 / 586. Regenerated Figures
10/12 and checked their styles visually. Figures 7–9, all-depth regime counts,
gradient summaries and map-bin tables remain byte-for-byte unchanged.

## Census implementation and style details (2026-11-16)

Presentation rule (user correction): preserve the original tagged figures' style
for every topic: layouts, sizes, colours, fonts, legends, axes, histogram/box
construction and map arrangement. Agreed scientific changes may change plotted
values, not the presentation by default. Discuss any necessary style departure.
Correctness takes priority over visual matching: raise issues for discussion
rather than reproducing an error silently; implement approved corrections.
Use this existing document as the decision/progress record; no extra log is needed.

Census style restoration: `census/census_plots.py` contains the required original
plotting routines adapted locally from the tagged notebook cells and shared
legacy helper. It has no legacy runtime imports. Figures 4/5 restore original
dimensions, outlined boxes, mirrored histograms, polar legends and mean arrows.
Figure 6 restores model-grid maps, region colours/boundaries and inset roses.
Its -20 degree rotation is a display transformation into grid coordinates;
analytical bearings and tables remain geographic. The shared -33 degree region
definition, >=5 km direction selection and equal-day primary weighting remain.
JSON records figure dimensions, original radial scale and display rotation.
The manifest includes the plotting module hash as well as the analysis code.
Validation after restoration: visually compared all three PNGs with saved tagged
legacy notebook outputs; five scientific tests and fresh-kernel notebook execution
passed. Exported rose counts conserve each regional sample, original Figure 5
radial limits do not clip bars, and both module hashes match the run manifest.

User authorized the census topic. Implemented `census/census.py`, `census/census.json` and a
short `census/census.ipynb`; see README.md and implementation_plan.md. No legacy imports.
Figures 4–6, CSV summaries/audit and manuscript report are generated locally at
`census/outputs/` (ignored by Git). Run manifest records config and source/code
SHA-256 hashes. No derived scratch cache is required for this topic.

Validated against full source data: 127,426 days / 2,982 eddies, 101,656 valid
tilts, 91,880 directional days. AE/CE median tilts: 18.03/15.18 km. Five scientific
invariant tests pass; notebook passed fresh-kernel execution. Original notebook
kept without saved output clutter; executed validation copy is in /tmp.

Additional inconsistency resolved explicitly: old census regions use -33 degrees
latitude, old regional rose figure uses -33.5. New single JSON region definition
uses -33 for all census summaries/figures. Geographic bearings are retained for
analysis; Figure 6 uses the original model-grid display rotation. Counts and
differences are in report.md.
The new region grid matches the old -33 definition on every wet cell. It leaves
land cells unlabelled (old shelf expansion labels 46 dry cells). Two source days
are consequently unassigned, both without valid tilt; plotted samples unaffected.
Later topic implementations are recorded above. No commit or push was part of
the census implementation.

## Relocation completed 2026-11-16

User subsequently authorized moving the existing analysis into `legacy_analysis/`
as a reference snapshot. This relocation is complete: all 40 top-level entries
outside `manuscript/` and root `AGENTS.md` were moved, including hidden checkpoints
and local generated files. The root now contains those two entries plus
`legacy_analysis/`. Verified 481 relocated files/directories/symlinks against a
pre-move inventory (SHA-256 and size for files; link targets for symlinks).
Legacy file contents were not changed. Scratch products were not moved or rebuilt.
`info.md` source paths and root `AGENTS.md` navigation were updated.

All historical notebook/helper paths in this review now resolve under
`seacofs_eddy_tilt_analysis/legacy_analysis/`, unless explicitly identified as
manuscript files or sibling upstream dataset code. The independence and scratch
plans below are still proposals; no analysis implementation has begun. No Git
commit or push is part of this relocation task.

## User intent and authorization

- Build a plan first. No extraction, refactor, migration, archiving, cache rebuild,
  or scientific changes yet. User authorized notebook execution if needed and
  Markdown notes here. This review did not execute whole notebooks.
- `manuscript/` should be independent of the surrounding exploratory code.
  Copy/adapt necessary functions into it eventually; no imports from folders
  intended for archiving. Avoid building an elaborate package/framework.
- Topic-based notebooks; analysis and substantial plotting in Python modules;
  notebooks expose settings, brief explanations, calls and figures.
- Established eddy dataset remains the input. Recompute additional analysis
  products into new dedicated scratch locations. Preserve old caches for comparison.
- Independence should expose inconsistencies, not force new results to match old
  figures. Record causes of differences; scientific choices require discussion.
- `info.md` is a starting map, not an exhaustive dependency specification.

## Agreed decisions after review (supersede proposals below)

- PV gradients use the surface ESP parameters consistently: surface ellipse
  geometry and surface vorticity amplitude, retaining the Gaussian reconstruction
  and nonlinear spatial averaging. Do not substitute peak-depth vorticity.
- Figures 13/14 retain their current depth-by-depth |Omega| and axis-ratio
  comparisons, including the surface-Rc-normalised tilt version.
- Retain peak absolute vorticity as an additional variable for possible future
  use, with its signed value and selected depth; retain the existing 0–1000 m
  search interval. It does not drive the agreed primary PV calculation.
- No temporal smoothing of PV gradients or dominance ratios. Use each eddy-day's
  calculated magnitudes directly.
- Dominance is purely component magnitude: planetary if |g_plan| >= 2|g_topo|;
  topographic if |g_topo| >= 2|g_plan|; otherwise mixed among valid inputs.
  Invalid/missing inputs remain unclassified; zero/zero must not satisfy both
  dominance definitions. No bathymetry in classification.
- After reviewing the all-depth results, user restored a separate shared 3000 m
  depth selection for targeted populations: deep planetary h>=3000 m, shallow
  topographic h<=3000 m. Keep dominance magnitude-only and all-depth Figures 7–9.
  Excluded deep topographic/shallow planetary days retain their original category;
  they are not mixed. Record retained/excluded days and eddies explicitly.
- Retain TiltDis >= 5 km for directional analyses. Magnitude analyses retain all
  valid tilt distances, including valid values below 5 km. This threshold does
  not affect the PV dominance category itself. Apply alongside the relevant
  validity/QC requirements; it is not a replacement for them.
- Primary population summaries use equal weight per eddy-day, rather than equal
  total weight per eddy. User prefers to leave per-eddy analysis for future work.
  Whole-eddy bootstrap resampling remains appropriate for confidence intervals;
  keep the pooled eddy-day estimator inside each resample. IQR remains population
  spread, not a confidence interval. Retain the existing explicit per-eddy panels
  in Figures 4/5 as a comparison with eddy-day summaries; user considers their
  similarity useful support for the primary weighting choice. Reassess that
  similarity after regeneration rather than assuming it. Per-eddy depth-ranking logic in Figures 13/14
  also needs review separately from weighting the plotted curves.
- Figure 5 per-eddy direction retains the current equal-day circular mean of
  unit directions (mean sin/cos of daily bearing), not a magnitude-weighted mean
  tilt vector. Correct its caption/description during implementation. Which
  days enter this summary still follows the agreed directional threshold/QC.
- Figures 13/14 use fixed existing fitted levels corresponding approximately
  to 0, 181 and 515 m; remove automatic correlation-based depth ranking from
  the new primary workflow. Verify and record exact stored depth values rather
  than rounding them for equality matching. Keep matched Eddy-Day samples
  across the selected depths and relevant properties.
- User explicitly accepts changes to both code and manuscript where needed
  to improve analysis. Neither current code nor prose is automatically the
  scientific authority; agreed definitions govern reconciliation. Do not force
  regenerated results to preserve existing conclusions.
- Figure 15 ellipse-depth comparison uses actual fitted levels, with matched
  Eddy-Day observations across compared depths; no interpolation of ellipse
  geometry. User approved five existing levels: approximately 0, 49, 106, 181,
  327 m; exact values are recorded in controls.json and the current status above.
- Figure 15 directional geometry requires 1.1 <= axis ratio <= 5 alongside
  TiltDis >= 5 km and valid geometry. Apply the axis-ratio requirement at each
  compared depth before forming the matched Eddy-Day sample.
- Figure 17 composites retain actual fitted centres at 0–1000 m inclusive,
  without interpolation. Reference each day's centres to its own surface centre
  and rotate into along-total-environmental-PV and +CCW perpendicular components
  before averaging. Equal eddy-day weights at each depth; whole-eddy bootstrap
  confidence intervals; report contributing days/eddies by depth. Directional
  selection uses saved TiltDis >= 5 km. Regimes follow the newly agreed raw
  magnitude-only classification. Future targeted planetary/topographic comparisons
  should use the agreed deep/shallow selection separately; mixed remains purely
  the magnitude-defined mixed category, not a complement of depth selections.
- Figure 11 background: retain centred 91-day seasonal climatology as primary,
  using monthly climatological model fields; retain surface, 0–200 m and
  200–500 m velocities and primary surface-minus-200–500 m zonal contrast.
  Rebuild in fresh scratch and use the deep-planetary selection (h>=3000 m).
  Validate layer coverage and report exclusions separately.
- Agreed small robustness comparison: seasonal windows 31, 61 and 91 days,
  holding depth ranges and Eddy-Day sample fixed. All use the same monthly
  climatological fields; shorter windows do not provide daily temporal resolution.
  Inspect vertical coverage before finalising retained depth ranges. Keep this
  comparison compact, for internal checking only; paper inclusion is not required.
  Primary window remains 91 days. Discuss any material sensitivity before finalising
  interpretation; do not tune settings for stronger results.
  No additional depth-boundary sweep or broader background study is agreed.
- Scope: do not recreate hypothesis_schematic.png (Figure 1) or
  dtheta_dom_cat.png (after composites). Include supporting manuscript statistics
  and tables. This workflow must support ongoing manuscript writing: user can
  update analysis here, generate figures and obtain the information needed to
  revise the manuscript, including outputs suitable for handing to ChatGPT.
- Each analysis/test should eventually produce a compact manuscript handoff:
  question/definition, settings and sample/exclusion counts, numerical estimates
  with units and uncertainty where appropriate, figure/table references, and
  interpretation limits. Preserve machine-readable results alongside readable
  summaries and source/run provenance. Exact format is a plan detail, not yet
  implemented. Distinguish descriptive summaries from inferential tests.
- Apply TiltDis >= 5 km to Figure 5 direction distributions and Figure 6 regional
  tilt roses too. Magnitude census retains all valid tilt distances. Figure 5
  per-eddy direction uses the qualified daily directions, followed by equal-day
  circular averaging. PV-gradient-only distributions do not automatically acquire
  a tilt-distance filter: they describe environmental gradients, not tilt bearings.
- Figure 3: user agrees to investigate which of the three tagged schematic
  examples produced the current paper image, then pin its Eddy-Day. No example
  identity is established by this agreement alone. Subsequent read-only visual
  inspection of GitHub Figures/delta_schematic.png shows tilt 3.7 km and maximum
  projected extent 21.1 km. These match the first saved tagged example:
  Eddy 2103, Day 7923 (saved 3.69 km and 21.070723 km), uniquely among the
  three outputs. That pair is now pinned and freshly rerun in `tilt_method/`;
  see the current status above. The >=5 km directional sample
  threshold does not prohibit this methods illustration of a smaller tilt.
- These decisions authorize planning/documentation only, not implementation.
  Remaining population/QC/weighting decisions are still open.

## What was inspected / evidence limits

Read local AGENTS.md, figure-producing notebook source, relevant helper functions,
saved outputs and selected workflow READMEs. Parsed notebooks as JSON; cell
indices below are **zero-based**. Figure markers are source comments, not notebook
metadata tags. A tagged cell with execution_count=null is not proof it has never
run; tags/source edits and saved outputs can have different histories.

GitHub access works via connected GitHub tools. Read default-branch
`agujournaltemplate.tex` from:
https://github.com/regdowse/Tilting-of-Mesoscale-Ocean-Eddies-in-The-Tasman-Sea
The manuscript has active figures outside info.md: hypothesis_schematic.png
(Figure 1) and dtheta_dom_cat.png (after composites). Their local provenance is
not traced in this pass. Do not silently include/exclude them in final scope.
No image-by-image comparison or LaTeX compilation was performed.

Scratch input files are visible on host `k180`. System python3 lacks the science
packages; existing `/home/z5297792/myenv310/bin/python` has NumPy, pandas,
Matplotlib, SciPy, netCDF4, pyarrow and nbclient. Kernelspec
`~/.local/share/jupyter/kernels/myenv310/kernel.json` points there. The empty
`~/.venv` directory is not a usable environment. No packages were installed.

Read-only pandas checks were run using that environment; no original notebooks,
figures, scratch data or caches were written. Permissions currently allow scratch
reads, but new scratch writes would require appropriate approval/access later.
Do not infer run provenance from mtimes alone.

## Base inputs (keep as inputs)

Let `D = /srv/scratch/z5297792/SEACOFS_26yr_eddy_dataset_modular`.

| Input | Current location | Role |
|---|---|---|
| Surface eddy-days | `D/processed/eddy_dataset_processed.parquet` | Geometry, position, polarity, track attributes |
| Saved tilt | `D/tilt/tilt_dataset.parquet` | Authoritative TiltDis/TiltDir |
| Confirmed vertical fits | `D/vertical_profiles_confirmed/profiles.parquet` | Fitted centres, ellipses, vorticity, Omega at depth |
| Physical sigma depths | `D/z_r.npy` | 3-D depth geometry for model layer averages |
| Grid reference | `/srv/scratch/z3533156/26year_BRAN2020/outer_avg_01461.nc` | Coordinates, mask, h, f, angle |
| Model archive | Same directory, `outer_avg_*.nc` | Required for background climatology rebuild |

Existing loader: `seacofs_tilt_tools.Paths`, `load_grid`, `load_tilt_tables`,
`load_vert`, `add_region_labels`. `load_tilt_tables` left-merges saved tilt onto
surface by `(Eddy, Day)` but does not itself validate merge cardinality.
Grid reader transposes horizontal dimensions; constructs x/y in km and uses
scalar `angle[0,0]`. Preserve/verify those conventions during extraction.

Numerical baseline checked on 2026-11-16:

- Surface: 127,426 rows, 2,982 unique Eddy IDs; unique Eddy-Day keys.
- Tilt: 115,156 rows, 2,982 IDs; unique keys; **101,656 finite TiltDis**,
  matching manuscript abstract. Remaining rows can have invalid/missing tilt.
- Shared PV: 127,426 rows, 2,982 IDs; unique keys.
- Outer joining tilt and PV keys: 115,156 shared, 12,270 PV-only, zero tilt-only.
  TiltDis and TiltDir match including NaNs on the union (zero differences using
  np.isclose). Cache staleness is a design risk, **not an observed tilt mismatch**.

## Figure-to-code map

Paths relative to `seacofs_eddy_tilt_analysis/`. Dependency names identify the
relevant existing code, not files to copy wholesale.

| Figure / file | Notebook : cell | Necessary upstream work / helpers |
|---|---|---|
| 2 pv_measurement_demo.png | `surface_pv_esp_gaussian/07_pv_gradient_example_selector.ipynb`: 7 | `esp_pv_tools`, shared PV geometry kernels; separate comparison cache; grid and vertical centres; explicit pairs (613,3370), (76,1675) |
| 3 delta_schematic.png | `delta_tilt_method/delta_tilt_schematic.ipynb`: 6 | `delta_sensitivity_tools.increments`, `fit_snapshot`; plotting is inline cell 4; profile input, example selection cells 1–2 |
| 4 tilt_census.png | `census/tilt_census.ipynb`: 11 | Cells 1,5: surface+tilt+regions, equal-eddy weighted density/boxes, regional per-eddy summaries |
| 5 tilt_dir_census.png | same: 12 | Same base; windrose and mirrored histogram helpers; circular per-eddy mean |
| 6 sectioned_rose.png | `rose_plots.ipynb`: 8 | Base+grid, `rose_plot`, windrose and region geometry |
| 7 pv_gradient_rose.png | `pv_tilt.ipynb`: 17 | Shared surface PV; natural log magnitude, geographic bearing, rose helper |
| 8 beta_dom_map.png | same: 22 | Shared PV; spatial binned mean of log topo/plan ratio, dominance contours, bathymetry |
| 9 topo_plan_ratio.png | `plan_topo_dom_eddies/plan_topo_mixed.ipynb`: 5 | Shared PV and smoothed ratio; histograms, regime proportions, magnitude boxes. Cell comment still calls file beta_dom_map.png; user's info.md corrects name |
| 10 plan_dom_rose.png | `plan_topo_dom_eddies/planetary_dominated_eddies.ipynb`: 6 | Shared PV; `case_study_tools.PVAlignmentConfig`, `add_pv_alignment_diagnostics`; actual population set in cell 3 |
| 11 plan_back_flow.png | `beta_effect_background_flow/planetary_eddy_vertical_flow_shear.ipynb`: 15 | Shared PV, alignment diagnostics, background cache, layer means cell 5, whole-track bootstrap summaries cell 7 |
| 12 topo_dom_dtheta.png | `plan_topo_dom_eddies/topographic_dominated_eddies.ipynb`: 8 | Shared PV, smoothed dominance, <=3000 m and tilt>=5; longitude shelf split, histograms and PV boxplots |
| 13 eddy_tilt_control.png | `eddy_env_tilt_controls/paper_eddy_controls.ipynb`: 5 | `paper_control_tools`: fit QC, candidate depths, matched days, descriptive depth ranking, binned medians/IQR |
| 14 eddy_tilt_control_Rc.png | same: 7 | Same matched observations, outcome TiltDis / surface Rc |
| 15 ellipse_tilt2.png | `ellipse_tilt_analysis/ellipse_orientation_tilt.ipynb`: 34 | `ellipse_tilt_tools`, `publication_tools`; all_shapes setup and PUB settings cell 29; six-panel direction figure |
| 16 env_tilt_control.png | `eddy_env_tilt_controls/paper_environment_controls.ipynb`: 12 | `paper_control_tools`, current base tilt joined to PV covariates; latitude/beta/log10 PV; raw and Rc-normalised outcomes |
| 17 composites.png | `esp_population_composites/pv_dominance_tilt.ipynb`: 5 | `pv_dominance_tilt_tools` -> planetary/topographic composite helpers -> `composite_breakdown_tools.centre_statistics`; fitted centres, no velocity reconstruction |

Companion paper statistics: `census/tilt_stats.ipynb` and `tilt_stats_tools.py`.
Include supporting numbers/tables in eventual scope, not just raster figures.

## Derived products and actual build paths

### A. Shared surface PV (most PV figures)

Existing: `D/pv_gradient_surface/surface_pv.parquet` (~85.6 MB).
Builder is `seacofs_tilt_tools.add_pv_gradient_terms`; explicit settings used by
these notebooks and confirmed in current Parquet attrs:

```
version=2; core_mean=True; frac=1; inner_frac=0;
averaging='nonlinear'; surface_method='esp_gaussian';
use_max_abs_w=True; max_depth_m=1000
```

Signed largest-|w| amplitude is selected per Eddy-Day from 0–1000 m (shallowest
tie); surface ellipse geometry is retained. Reconstruct local zeta as the
Laplacian of the ESP Gaussian streamfunction, evaluate local environmental PV
gradient, then take Gaussian-weighted vector means over the core. Standard
PV_grad excludes internal grad(zeta)/h; full/eddy terms are separate diagnostics.
New workflow needs explicit source profiles, avoiding implicit load_vert()
fallbacks, and must preserve source identity plus peak-depth diagnostics.

`use_cache=True` returns the **entire saved dataframe**, ignoring supplied df/grid.
It can therefore replace current tilt/region/geometry with cached copies.
Settings comparison does not establish current input identity. False calculates
and atomically replaces one configured file; None calculates without saving.
Do not toggle False on existing notebooks for this migration.

Current cache attrs timestamp: 2026-09-25T03:03:46.143204+00:00.
Stored ratio smoothing matches current centred **3 observations, minimum 3**
exactly, for both polarities. AGENTS.md's historical 7/min5 note is stale for
these sources: alignment config also defaults to 3 and clips min_periods=5 to 3.
Observation smoothing can bridge calendar gaps; not equivalent to 3 consecutive days.

### B. Figure 2 PV comparison/examples

Existing: `D/pv_gradient_surface_esp_gaussian/surface_pv_esp_gaussian_laplacian_v2.parquet`
(~143.5 MB across both methods). Builder:
`surface_pv_esp_gaussian/00_build_esp_gaussian_cache.ipynb` ->
`esp_pv_tools.build_cache`. Builds uniform_1 and esp_gaussian_1 on FRAC=1, but
does **not** pass use_max_abs_w: default surface w is used. Candidate selector
uses esp_gaussian_1. Current comparison cache has no peak-w diagnostic columns
and attrs are empty. Local example maps reconstruct from row.w too.

Read-only confirmed example values:

| Eddy,Day | Polarity | Example cache PV magnitude | Shared peak-w PV magnitude | Bearing example -> shared |
|---|---|---|---|---|
| 613,3370 | CE | 1.14667432895e-11 | same | 88.8993 -> same degrees |
| 76,1675 | AE | 7.48291583233e-14 | 5.13632511250e-14 | 70.3022 -> 58.6784 degrees |

AE surface w=3.11699198976e-5 s^-1; shared selected w=4.35994698030e-5 at
126.288126687 m. CE selected depth is 0. The manuscript Figure 2 caption values
match the **example cache**. Its PV methods text describes surface-fit vorticity;
later workflows use peak-w amplitude. Discuss a common definition or explicitly
document the distinct demonstration. This is a real numerical discontinuity.
Explicit example selection currently searches a screened candidate pool; final
reproduction should select fixed Eddy-Day directly from a specified product so
unrelated candidate ranking changes cannot remove the example.

### C. Background flow (Figure 11)

Existing root is under the **non-modular** dataset parent:
`/srv/scratch/z5297792/SEACOFS_26yr_eddy_dataset/background_flow_cache_all_eddies_v4/`.
Table `eddy_day_background.parquet` (~11.8 MB); intermediates
`file_partitions/`, `monthly_climatology.npz`.
`01_build_background_cache.ipynb` / `build_background_cache.py` launch
`background_flow_tools.build_background_cache`.

Reads whole archive, parallel file partitions, monthly/full sums and counts;
builds centred 91-day seasonal climatology by month-overlap day weighting.
This is not the event year's daily running mean. Local sigma-layer thickness
weights yield surface, 0–200 and 0–500 m velocities; derive 200–500 only with
valid coverage. Figure uses surface minus 200–500 m zonal contrast in m/s.
Rebuild in a fresh root including fresh partitions; do not reuse old partitions
as evidence of independent regeneration. This is likely the expensive stage;
no runtime benchmark or archive completeness check performed here.

Figure 11 needs only the selected population, layer contrasts and cluster-median
summary, not all exploratory track propagation/regression/depth-motion cells.

### D. Other derived analysis

Persist useful small products in the new scratch space: selection/exclusion
tables, matched-depth inventory, bins/medians/IQR, bootstrap summaries, composite
members and support counts. No N2, wind, model vertical-velocity or depth-following
PV build was found necessary for the listed primary figures. Avoid automatically
porting their unrelated workflows. Figure 17 needs existing fit centres only;
it does not need the external ESP velocity backend or reconstruction caches.

## Discontinuities / decisions to preserve

1. **PV amplitude:** Figure 2 surface-w vs later peak-w, confirmed above. Highest
   priority definition decision; changing it changes existing caption values.
2. **Population consistency:** Figure 10 cell 3 uses ratio_smooth<=-ln2 and
   TiltDis>=5, without active h>=3000 filter; caption/printout say depth filter.
   Figure 11 uses open_ocean_planetary & direction_valid, which does apply it.
   Read-only current check: Figure 10 selects AE 4053, CE 1395; zero selected
   rows have h<3000 in either polarity. Thus missing depth condition is real in
   source but currently has no depth-exclusion effect. Full Fig10/11 sample
   equality still requires checking direction validity/background coverage.
3. **Smoothing/regimes:** Figures 9/10/11/12 use smoothed ratio (current 3/min3).
   Figure 17 uses unsmoothed component magnitudes and strict depth boundaries
   (>3000 planetary, <3000 topographic); exactly 3000 is mixed. Its mixed class
   is complement of both complete selections, so ratio-dominant days in the
   opposite depth range can be mixed. Figure 9 uses ratio-only categories, no
   depth restriction; initialise-as-mixed can count invalid smoothed rows as
   mixed. Document intended distinctions or harmonize after discussion.
4. **Weighting:** Figure 4 daily distributions use inverse per-eddy day counts
   within the regional sample (equal total weight per eddy), not unweighted
   eddy-days. Figure 5 daily histogram uses pooled days; per-eddy direction is
   circular mean of unit directions, not magnitude-weighted mean tilt vectors.
   Caption describes mean zonal/meridional tilt components: clarify intended
   estimator. Figure 17 equal-day mean with whole-track bootstrap differs from
   Figure 15 equal-eddy estimates. Do not impose universal weighting blindly.
5. **Figure 15 tag is authoritative:** cell 29 deletes depth panels to form a
   2x2 plot, whereas tagged cell 34 produces the six-panel plot matching the
   paper caption. Its name says observation_weighted but equal_eddy=True.
   Use actual call, not name. Depth geometry interpolates Q entries with <=100 m
   brackets at nominal 0/50/100/200/300, unlike exact fitted-level controls.
6. **Rotation metric caption:** Figure 13 caption says absolute rotational
   velocity; code is abs(Omega), units s^-1 (rotation parameter), not tangential
   speed. Confirm precise scientific terminology before carrying caption over.
7. **Delta example identity:** tagged Figure 3 cell loops over three random
   examples (saved pairs: 2103/7923, 2000/7636, 1864/7169; seed 5). Which image
   was exported is not pinned by tag alone. Pin exact pair after visual check.
   Notebook prose says production unchanged, but upstream core/tilt.py now
   implements reference-day support restriction; production defaults are
   centred 5-day, sigma=1, bearing offset=20 degrees. Base tilt regeneration
   is outside current scope; provenance must still be recorded.
8. **Exports:** tagged final panels sometimes differ from objects saved by
   old optional export blocks (notably combined Figure 16). New workflow needs
   explicit filename-to-object mapping. Figure 9 old comment filename differs
   from updated info.md; do not reproduce that naming mistake.

## Proposed minimal structure (discussion only)

Flat modules + topic notebooks initially; avoid a new framework or large copied
monolith. Candidate topics:

- `methods.py/.ipynb`: PV illustration + delta schematic (Figures 2–3).
- `census.py/.ipynb`: census, directions, regional roses, paper counts (4–6).
- `pv.py/.ipynb`: PV build/calculation and spatial/regime figures (7–10,12).
- `background.py/.ipynb`: background build, layers and Figure 11.
- `controls.py/.ipynb`: eddy/environment/ellipse relationships (13–16).
- `composites.py/.ipynb`: PV-relative centrelines (17).
- Small `data.py` for source paths/readers/key checks/geography shared across
  topics; split focused geometry helpers only if size warrants it. Figure
  methods call the same PV kernel as population figures.
- `figures/` for regenerated paper exports, README/index and this review.

Exact filenames/groupings are provisional. A build script can be short Python;
no need for a notebook per cache or automatic expensive builds on import.
Notebooks should run top-to-bottom without exploratory cells or hidden globals.

Proposed dedicated scratch root (not created):
`/srv/scratch/z5297792/eddy_tilt_manuscript/`
with `surface_pv/`, `background_flow/`, `summaries/`. Add an explicit run folder
when a different input/settings version is intentionally rebuilt. Keep a small
JSON manifest: source paths/identity, code revision, settings, units/conventions,
row counts, creation time. Avoid timestamp-only output discovery.

New PV product should contain keys and derived covariates; join onto current
base tilt explicitly with cardinality/coverage checks. This prevents cache loads
from silently defining source outcomes. Preserve source and derived h/w meanings.

## Suggested implementation sequence after discussion

1. Resolve PV amplitude and regime/smoothing/weighting choices; pin schematic
   identity and confirm scope for Figure 1 / post-composite figure / paper tables.
2. Record immutable baseline metadata, old figure files and key summaries for
   comparison; source dataset is accepted as input, not rebuilt.
3. Implement minimal independent loaders and common geometry; validate joins,
   units and bearings against current sources.
4. Extract/review PV kernels and fresh builder; rebuild only into new scratch;
   compare on identical Eddy-Day keys. Classify discrepancies as intended
   definition change, input change, numerical difference or defect.
5. Extract cheap census/PV/controls/composite topics and persist selection/support
   tables. Preserve meaningful bootstrap tests/invariants rather than exploratory
   plotting machinery. Smoke-run clean notebooks in the existing environment.
6. Extract/rebuild expensive background pipeline separately, then Figure 11.
7. Compare figures/captions/numbers; verify no imports or file dependencies on
   old analysis directories. Only then consider user's later archive step.

No implementation, cache rebuild, archiving, commit or push performed in this pass.

## Figure outputs (2026-11-18)

Final PNGs now go to `manuscript/figures/` via each topic's `figure_dir` setting.
Moved 16 existing PNGs unchanged, added the existing conceptual Figure 1, and
removed the 16 generated figure PDFs. All seven generators now save PNG only.
Topic outputs retain numerical tables, reports and run manifests. Historic
manifests still describe their original run; no scientific calculation or scratch
cache was rebuilt. Manuscript/latexdiff PDFs and assets inside the independent
text repositories were retained. The shared figures folder is Git-ignored.
Verified all seven configured export destinations, unchanged hashes for the moved
PNGs, no remaining topic figure PNG/PDF copies, and 49 passing existing tests.

## Land-cell region fallback (2026-11-18)

`census.prepare` retains its nearest-grid-cell lookup. Only when the selected
cell is land, it uses the wet cell closest to that land cell in model x/y, rather
than searching from the fitted eddy centre. Missing centres remain unassigned.
A regression case distinguishes those two search origins. All 51 tests pass.

Full-source comparison: only AE Eddy 2990 / Day 10265 changes, from unassigned to
U2/Upstream. Source U2 AE-days become 11,686 (total 26,280); unassigned days become
zero. Every valid-tilt and directional row, all tilt tables, and all other regional
assignments are identical. This day has valid PV direction, so Figure 7's regional
AE sample increases from 67,247 to 67,248 and its U2 rose gains one observation.
Only `pv_gradient_rose.png` changes; all other PNGs are byte-identical.

Census and PV-dominance notebooks were re-executed on compute node k001 by PBS
job 9821282 (exit 0). Reports, tables and manifests are refreshed; scratch products
were reused. Other workflows do not use the changed regional membership in their
calculations and were not rerun. Verification and before-run outputs are retained
in `context/revision/outputs/region-fix/`.

The paper review branch updates its source-population table and Figure 7 asset;
no prose or scientific interpretation changes. Cumulative diff is in
`context/revision/outputs/review/main-paper.diff`, against the latest-pull baseline
`0f173bb`. Changes remain local for review; no PDF was requested or generated.

Region correction publication: user approved the manuscript push on 2026-11-18.
Committed table/Figure 7 updates as `cbd1c3e37d160708063fbfeed64da543d8154eed` and pushed
to GitHub `main`. Local checkout is clean; review baseline remains unchanged
until Overleaf synchronisation is confirmed.
