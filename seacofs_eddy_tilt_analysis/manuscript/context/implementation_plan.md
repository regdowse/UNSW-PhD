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
