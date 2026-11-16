# Implementation plan

Layout: each analysis has its own folder; census code/config/notebook are in
`census/`, with generated files in `census/outputs/`. Paths in this plan are
relative to `manuscript/`. Planning Markdown is stored in `context/`. Automated
checks share `tests/`, which is Git-ignored as requested.

Build one topic at a time. Keep Python calculations separate from short notebooks.
Use JSON settings, CSV results, figures and concise Markdown reports.
Preserve the original tagged figure styles; change presentation only after discussion.
Adapt needed plotting routines locally and compare against saved original figures.

| Order | Topic | Outputs |
|---|---|---|
| 1 | Census | Figures 4–6, counts, distributions and direction summaries |
| 2 | Surface PV | Fresh surface-fit PV products; method illustration |
| 3 | PV dominance | Unsmoothed magnitude-only categories and comparisons |
| 4 | Controls and composites | Fixed fitted depths, matched samples, Figures 13–17 |
| 5 | Background flow | Fresh seasonal climatology; internal 31/61/91-day comparison |

Agreed definitions are in `workflow_review.md`. No legacy imports. Source datasets
remain inputs; future derived products use a new scratch root. Do not recreate
the hypothesis schematic or the additional post-composite figure.

Census: `census/census.json` → `census/census.py` → `census/census.ipynb`.
Outputs go to `census/outputs/`: figures, CSV tables, `report.md`, and a manifest
with the effective settings, source/code hashes, environment and Git revision.
Rerunning replaces this topic's outputs. Copy an output folder to retain a run.
No PV or background cache is needed for census.

Checks: unique join keys and coverage, finite values, circular bearings, weighting,
threshold boundaries, shared region definitions, full-data execution and figure
inspection. Descriptive census summaries do not require significance tests.

Surface PV: `pv/pv.json` → `pv/pv.py` with local scientific kernels and original
Figure 2 plotting helpers. Build explicitly with `python pv/pv.py --build` from
manuscript; then run `pv/pv.ipynb` to validate/load and export. Existing products
require explicit `--build --rebuild` to replace. Default fresh product root:
`/srv/scratch/z5297792/eddy_tilt_manuscript/surface_pv/surface_fit_v1/`.
Local figure/report/CSV outputs: `pv/outputs/`. No legacy runtime imports.
Only the surface Gaussian method is built; unrelated uniform-control comparisons
and PV-dominance figures are outside this step.

PV dominance: `pv_dominance/dominance.json`, Python analysis/plotting modules and
short notebook. Figures 7–10/12, writing report and supporting CSVs are in
`pv_dominance/outputs/`. Reads and validates the fresh surface-PV product; no new
scratch computation. Uses shared manuscript census geography/rose styling.
Approved display corrections: rotate Figure 7 bearings into model-grid axes;
expand Figure 9 distribution limits and Figure 10 radial limits only as needed.
Figure 11 is implemented separately in the background-flow topic.
Keep magnitude-only categories distinct from target selection: Figure 10 uses
deep planetary days (footprint-mean h>=3000 m), Figure 12 shallow topographic
days (h<=3000 m). Figures 7–9 show all depths. A single configurable threshold
controls both target filters; excluded days keep their dominance category.
Reports and depth_selection_counts.csv distinguish depth exclusions from other QC.

Controls: `controls/controls.json`, short notebook and local Python helpers retain
the tagged Figures 13–16 styles. Structure depths are exactly 0, 181.26479327,
515.41648935 m; ellipse depths (approved) are 0, 48.89763802, 105.61553821,
181.26479327, 327.44085148 m. Fixed fits only; no interpolation/ranking. Plotted
statistics weight eddy-days equally; Figure 15 intervals resample whole tracks.
Figure 16 raw y-limit is 50 km (approved to display full IQR); normalised limit
remains 0.6. Local outputs include CSV estimates/support/audits and report.md.
Figure 13 shared y-axis expansion is also approved to show full IQR bands;
`figure13_expand_to_iqr` records this presentation correction in JSON.

Composites: `composites/composites.json`, short notebook and local Python helpers
produce Figure 17 in its original six-panel style. Select deep planetary and
shallow topographic days using the shared 3000 m threshold; mixed remains purely
magnitude-defined. Use valid actual fitted centres within 0–1000 m, referenced to
each day's own fitted surface centre, rotated into that day's total signed-PV
frame before averaging. Positive perpendicular is counterclockwise. Equal-day
means and pointwise whole-track bootstrap CIs (500 draws, seed 731); counts vary
with depth and are exported. No interpolation or velocity reconstruction.
The local output folder contains figures, statistics/support, selection audit,
members, writing report and manifest. Background flow (Figure 11) is implemented
below; Figure 3's tilt-method illustration is implemented in `tilt_method/`.

Background: `background/background.json`, short notebook, local monthly/layer and
bootstrap kernels, and original Figure 11 plotting style. Build the full archive
into fresh `background_flow/monthly_v1/` scratch space, including new file partitions,
monthly sums/counts and all-source eddy-day tables for windows 31/61/91 days.
Primary 91-day centred seasonal climatology; fixed surface, 0–200 and 200–500 m
layers. Preserve the original physical sigma-centre midpoint-interface approximation
and require full finite layer coverage. Use the same complete deep-planetary
directional sample for every window and displayed layer. Pooled eddy-day medians,
whole-track bootstrap (5000 draws, seed 42 with exported offsets). Internal window
checks stay in compact tables, with paired differences from 91 days; no extra paper
panels. Preserve legacy products and reject changed product inputs/code/content.
Executed against all 307 archive files (9,189 consecutive days); all 6,103 target
days have complete layer/month/window coverage. The 31/61/91-day population contrast
is consistently positive; detailed checks are in the topic's local outputs.
Fresh-kernel execution and original-style figure comparison passed.

Tilt method: `tilt_method/tilt_method.json`, short notebook, local delta kernels
and original Figure 3 plotting function. Fixed Eddy 2103 / Day 7923; five calendar
days, Gaussian sigma 1 day, 10 m increments, cap 1000 m, reference-supported
intervals only, inverse unweighted sample-variance depth weights and the existing
weighted 3-D SVD convention. This method retains within-profile interpolation,
unlike the exact-depth control/composite comparisons. No extrapolation. Export
source profiles, temporal increments, reconstructed trace/weights/support, line
endpoints, projected profiles, example summary, report and manifest. Saved source
tilt and scratch products stay unchanged. Recomputed example agrees with saved
tilt and the tagged original; Figures 2–17 are now recreated independently.
