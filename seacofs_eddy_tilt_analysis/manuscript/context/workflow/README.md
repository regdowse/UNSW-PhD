# Manuscript analysis

Usage guide. Commands and inline code paths below are relative to `manuscript/`.
For manuscript revision, start at [context/README.md](../README.md).

Start with [census/census.ipynb](../../census/census.ipynb). Edit settings in
[census/census.json](../../census/census.json), then run all
cells. From the manuscript directory, the same analysis runs with:

```bash
python census/census.py --config census/census.json
```

Use a Python environment with NumPy, pandas, SciPy, Matplotlib, netCDF4 and pyarrow;
Jupyter/nbformat/nbclient are needed for notebook use. The existing Katana kernel
`myenv310` has these dependencies.

`census/outputs/` contains Figures 4–6 (PNG/PDF), CSV tables, a concise `report.md`
for manuscript updates and `manifest.json` with the effective configuration and
source/code hashes. Each run replaces these files. Copy the folder to retain a run.
Generated outputs are local artifacts, excluded from Git; source and config are tracked.

All summaries use one shared region definition. Primary distributions weight
eddy-days equally; per-eddy panels remain comparisons. Directions require tilt
≥5 km. No PV cache or legacy imports are needed.
`census/census_plots.py` preserves the original figures' presentation locally;
future topics should preserve their tagged legacy styles too.

See [implementation_plan.md](implementation_plan.md) for sequence;
[workflow_review.md](workflow_review.md) records decisions
and historical differences. The census report is descriptive, not a set of
independent-observation significance tests.

Layout: one folder per analysis topic (`census/` first); planning and navigation
Markdown files in `context/workflow/`; automated checks in shared `tests/` (Git-ignored).
Generated reports remain with their topic's outputs. Run checks from manuscript:
`python -m unittest discover -s tests -p 'test_*.py' -v`.

Surface PV: settings `pv/pv.json`, notebook `pv/pv.ipynb`. From manuscript, build
once with `python pv/pv.py --build`; subsequent runs use `python pv/pv.py` or the
notebook. Input/config/code hashes are checked before loading the scratch product.
Only scientific dependencies invalidate that product; plot/report edits do not.
Figure 2 rotates geographic arrows into model-grid coordinates for display.
An intentional replacement requires `--build --rebuild`. Figure 2 and its writing
report are in `pv/outputs/`. Both surface PV and auxiliary peak-vorticity fields
are saved; only surface vorticity enters the PV calculation.

PV dominance: run `pv_dominance/dominance.ipynb`, or
`python pv_dominance/dominance.py` from manuscript. Settings are in the matching
JSON; outputs include Figures 7–10/12, `report.md`, counts, directional summaries,
map-bin support and an audit. Shared census/PV modules are inside manuscript;
no legacy functions or old PV caches are used.
Dominance is magnitude-only. Figures 10/12 additionally select deep planetary
(h>=3000 m) / shallow topographic (h<=3000 m) days, respectively; Figures 7–9
remain all-depth. `water_depth_threshold_m` sets the shared cutoff. The report,
audit and `depth_selection_counts.csv` show retained/excluded populations without
reclassifying excluded days as mixed.

Controls: run `controls/controls.ipynb`, or `python controls/controls.py` from
manuscript. Exact depths, QC, bootstrap and plotting settings are in `controls.json`.
Figures 13–16 and CSV estimates, support/audits and writing report are in
`controls/outputs/`. Reads existing source fits and validated PV; no legacy imports
or new scratch builds. Axis alignment uses day weighting and whole-track CIs.

Composites: run `composites/composites.ipynb`, or `python composites/composites.py`
from manuscript. Settings are in `composites.json`. Figure 17, depthwise means/CIs
and support, selection audit, constituent centres and writing report are in
`composites/outputs/`. Uses shared dominance/depth selections and actual fitted
centres within 0–1000 m, referenced to each day's own surface centre. Means weight
eddy-days equally; confidence intervals resample whole tracks. No interpolation,
legacy imports, ESP backend or new scratch product is required.

Background flow: settings `background/background.json`, notebook
`background/background.ipynb`. Build once with `python background/background.py --build`
from manuscript; the notebook or `python background/background.py` validates and
loads the independent scratch product. A matching interrupted build resumes its
new partitions; changed scientific inputs/code require a new `product_dir`.
The primary Figure 11, layer summaries, matched internal 31/61/91-day checks,
selection audit and writing report go to `background/outputs/`. This topic uses
monthly climatological model fields and the shared deep-planetary selection.
The archive build additionally uses joblib for file-level parallel processing.

Tilt method: run `tilt_method/tilt_method.ipynb`, or
`python tilt_method/tilt_method.py` from manuscript. Settings are in the matching
JSON. Figure 3 pins Eddy 2191 / Day 7923, recomputes the five-day delta fit from
confirmed centres and compares it with saved tilt. Original plotting style is
retained. Figure, concise writing report and calculation-step tables are in
`tilt_method/outputs/`. No scratch build or legacy/production imports are needed.
All requested figure workflows (Figures 2–17) are now available by topic.

Dataset refresh 2026-11-17: all seven notebooks are saved with executed inline
figures and updated reports/tables. See workflow_review.md for cache reuse and
example-ID remapping after the source dataset changed. Use the latest topic
reports rather than historical numerical values in the context notes.
