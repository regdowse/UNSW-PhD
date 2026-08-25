# SEACOFS Eddy Tilt Analysis

This folder is a cleaned working set for the core vertically checked eddy-tilt
analysis. It keeps the original exploratory notebooks intact in
`SEACOFS_dataset/detect_track_tilt_dataset/vert_checked_dataset/tilt_analysis_vert_checked/`
and collects the selected analysis notebooks here with shared utilities in
`seacofs_tilt_tools.py`. The helper module is self-contained and does not
depend on the old `clim_functions.py` module.

## Structure

- `seacofs_tilt_tools.py` contains shared loading, grid, PV-gradient,
  propagation, binning, histogram, windrose, regional rose-plot, and
  sample-matching helpers.
- `pv_tilt.ipynb` analyses alignment between tilt direction and PV-gradient
  terms.
- `secondary_tilt_factors.ipynb` compares tilt distance with vorticity,
  shape, stratification, and propagation diagnostics.
- `long_eddy_case_study.ipynb` selects long-lived, strongly propagating AE/CE
  examples and plots their tilt histories.
- `tilt_analysis_vc.ipynb` summarises the main vertically checked tilt-distance
  relationships.
- `tilt_census_vc.ipynb` gives counts and distribution summaries by cycle and
  region.
- `windrose_plots_vc.ipynb` builds direction/magnitude windrose summaries.
- `sample/` contains the sample-eddy notebooks and the old-to-new sample ID
  matcher.

The default paths in `seacofs_tilt_tools.py` match the current HPC scratch
locations used by the original notebooks. Override the `Paths` dataclass in a
notebook if the data move.
