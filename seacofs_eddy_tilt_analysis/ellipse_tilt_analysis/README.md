# Ellipse orientation, deformation and tilt

Run `ellipse_orientation_tilt.ipynb` on Katana, from this directory, the
parent analysis directory, or the `UNSW-PhD` repository root. It uses the
existing `seacofs_tilt_tools.Paths` defaults. Override those paths in the
settings cell if necessary. No new model-data cache is required.

The surface ellipse is the primary analysis; 50, 100, 200 and 300 m are
secondary shape comparisons against the same existing whole-column tilt.
Set `RUN_DEPTHS=False` for a surface-only first pass.

The notebook includes:

- Q-eigenvalue axis ratios and major-axis bearings aligned with the tilt
  pipeline's bearing convention, with a synthetic convention check.
- Separate AE/CE surface magnitude, alignment and elongation-class figures.
- Between-eddy correlations and equal-eddy within-eddy slopes with
  whole-eddy bootstrap intervals.
- Regional and PV-regime summaries, plus consecutive-day rotation diagnostics.
- Available-depth and matched-eddy-day depth comparisons.
- Directional threshold and nearest-depth versus interpolated-Q sensitivity.
- Optional CSV, PDF and 600 dpi PNG export with analysis settings.

`ellipse_tilt_tools.py` contains workflow-specific helpers. Shared loaders and
PV/region calculations remain in `../seacofs_tilt_tools.py`; the existing ML
feature definitions are not changed. The code uses NumPy, pandas, SciPy and
Matplotlib, plus the shared loader's netCDF4 and a pandas Parquet engine
(typically pyarrow). IPython is used only for notebook display.

The bearing offset defaults to 20 degrees, matching the modular pipeline's
`config/local.yaml` setting `tilt.bearing_offset_deg`. This is not an additional
rotation to apply to already geographic Q matrices. Verify the setting if
using a dataset generated with different coordinate conventions.

Only finite positive-definite Q matrices are accepted. Depth interpolation
uses adjacent valid matrices, never extrapolates, and limits bracket width.
The stored `DepthLower`, `DepthUpper` and `DepthMethod` fields record sampling.
Matched samples are constructed after the outcome-specific quality filters.

These are exploratory associations. Near-circular axes, small tilts, shared
fit errors, temporal smoothing, variable measured vertical span, repeated
eddy observations and geographic covariates require care. The notebook
documents these limitations explicitly; its uniform-angle line is a
descriptive reference, not a formal test of independence.

Run the small synthetic regression suite from this folder:

```bash
MPLBACKEND=Agg python3 -m unittest test_ellipse_tilt_tools.py
```

Production data are not bundled. Locally validated notebook execution uses
synthetic input and cannot establish the scientific result on Katana data.

## Publication section

Section 10 appends a six-panel alignment figure and a separate magnitude
figure, plus observation-weighted and stricter-QC alignment sensitivity
figures. `publication_tools.py` contains the associated plotting and bootstrap
helpers. The section needs only the earlier imports and data loading through
`all_shapes`; it does not depend on intermediate exploratory plots.

The primary direction sample requires tilt >=5 km and AR >=1.1; magnitude
plots retain small tilts and nearly circular ellipses. Estimates require
five valid days per eddy within each actual sample/AR class and at least
20 eddies. Unsupported bins are left blank and retain counts in the tables.
The depth panels use identical eddy-days at all depths. Paired
alignment-score differences relative to the surface are exported as well.
Both weighting schemes resample whole eddies for uncertainty.

`PUB_SAVE=True` writes vector PDFs (rasterized scatter points), 600 dpi PNGs,
CSV statistics and coverage, draft captions and settings to
`outputs/publication/` when run. Set it to `False` for a display-only pass.
These ignored output files are generated on Katana and are not committed.
Existing notebook cells and saved outputs were preserved when this section
was added. New publication cells are delivered unexecuted.

Section 11 adds separate 1x2 surface-angle distribution figures for all
regions, shelf, and off-shelf. Each AE panel uses red shades and each CE
panel blue shades, with darker lines for higher AR. The fixed 5-degree
histograms use equal-eddy weights within each region/class, pointwise
bootstrap bands, common y-limits, and explicit per-line sample counts.
Set `ARH_SHOW_CI=False` to hide bands or `ARH_SAVE=True` to export the three
figures, probability/coverage CSVs, captions and settings. No existing
publication figures need to be rerun for this section.
