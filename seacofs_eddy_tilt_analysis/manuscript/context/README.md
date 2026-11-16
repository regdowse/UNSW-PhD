# Manuscript analysis

Start with `../census/census.ipynb`. Edit settings in `../census/census.json`, then run all
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

See `implementation_plan.md` for sequence; `workflow_review.md` records decisions
and historical differences. The census report is descriptive, not a set of
independent-observation significance tests.

Layout: one folder per analysis topic (`census/` first); planning and navigation
Markdown files in `context/`; automated checks in shared `tests/` (Git-ignored).
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
