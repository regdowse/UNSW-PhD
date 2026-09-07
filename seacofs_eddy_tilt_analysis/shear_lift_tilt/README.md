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

`shear_lift_tools.py` contains only the new shear-coordinate geometry, exact-lag
matching and compact binned summaries. Positive cross-shear displacement means
left of the upper-minus-deep shear vector. Polarity is not embedded in that
definition, allowing the AE/CE reversal to be tested rather than assumed.
