# Eddy cross-sections and 3-D ESP reconstruction

Open `eddy_sections_and_esp.ipynb` from this directory on Katana. Change
`EDDY_ID` and, optionally, `DAY` in the controls cell. When `DAY=None`, the
notebook chooses the requested eddy's deepest, best-resolved fitted profile.

The notebook provides:

- a plan view showing the two section locations and depth-dependent fitted centres;
- true east–west and north–south sections of N², temperature, salinity,
  potential density, eastward velocity, northward velocity and speed, using
  horizontal interpolation but retaining native sigma levels vertically;
- the matching v4 N² cache summary for the selected eddy-day;
- depth-dependent ESP velocity reconstruction;
- original-versus-ESP vertical sections and 3-D velocity-vector views;
- a selected horizontal slice comparing original and ESP velocity fields and
  their relative-vorticity fields on common colour scales. The original
  vorticity is the numerical curl of the interpolated velocity, while the ESP
  vorticity is evaluated analytically from the fitted non-axisymmetric
  Gaussian streamfunction.

Set `HORIZONTAL_SLICE_DEPTH_M` in the controls cell to select the horizontal
comparison. The nearest depth with an ESP fit is used and reported above the
figure. `HORIZONTAL_QUIVER_STEP` controls the velocity-arrow density.

The N² cache stores reduced eddy-core statistics, not spatial fields. Section
N² is therefore calculated from the source `temp` and `salt` fields using the
same surface-referenced `xroms.potential_density(..., z=0)` and vertical
buoyancy-gradient method used to build the v4 cache.

ESP requires `/home/z5297792/ESP_zonodo/functions.py` by default. Change
`ESP_ROOT` in the notebook if that checkout moves.

`eddy_esp_composite.ipynb` reconstructs every selected day of one eddy and
then averages the evaluated ESP velocity fields. Its default onshore-aligned
frame recentres every day and rotates the local core-mean bathymetric gradient
to a common direction. This is intended for testing persistent shelfward tilt
without smearing the eddy as it translates southward. `esp_composite_tools.py`
contains the reusable compositing functions and also supports radius-normalised
coordinates for later multi-eddy comparisons.
