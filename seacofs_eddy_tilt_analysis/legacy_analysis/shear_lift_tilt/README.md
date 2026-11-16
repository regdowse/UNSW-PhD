# Rotation-dependent shear/lift test

`rotation_dependent_shear_lift.ipynb` tests whether environmental vertical
shear has a common along-shear effect and an opposing cross-shear effect for
EAC anticyclones and cyclones, as suggested by the Kutta-Zhukhovski
interpretation in Chao and Shaw (1998).

The notebook treats `TiltDir` and `TiltDis` as the authoritative coherent tilt
measurements. It never reconstructs tilt from raw daily centre positions.

It reuses:

- `seacofs_tilt_tools.py` for paths, grid data, tilt tables and regions;
- `tilt_mechanisms/mechanism_tools.py` for tilt-vector conversion;
- `beta_effect_background_flow/background_flow_tools.py` for the existing
  centred 91-day climatological and full-archive background-flow cache.

The primary comparison uses 0-200 m minus 200-500 m climatological flow so
that energetic surface variability is downweighted. Surface minus 200-500 m,
and full-archive means, are predefined sensitivity tests. Run
`beta_effect_background_flow/01_build_background_cache.ipynb` first if the v4
background cache is not already available on Katana.

The primary mechanism subset is open-ocean planetary dominated: the centred
seven-day median planetary PV-gradient magnitude must exceed the topographic
term by at least a factor of two and bathymetry must be at least 3000 m. Strongly
topographic days form the explicit complication/comparison group; mixed days
are audited but not pooled with either extreme. Additional models distinguish
within-eddy from between-eddy shear, allow a nonlinear shear response, test
scaling with fitted rotation strength, and control lagged tests for the current
tilt state, latitude and region.

`shear_lift_tools.py` contains only the new shear-coordinate geometry, exact-lag
matching and compact binned summaries. Positive cross-shear displacement means
left of the upper-minus-deep shear vector. Polarity is not embedded in that
definition, allowing the AE/CE reversal to be tested rather than assumed.
