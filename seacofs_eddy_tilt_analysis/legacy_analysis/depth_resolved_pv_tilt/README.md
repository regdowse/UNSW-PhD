# Depth-resolved PV-gradient and eddy-tilt analysis

This visual-first notebook series tests whether the effective topographic
PV-gradient environment sampled along a displaced eddy spine differs from the
environment inferred beneath its surface centre, and whether the depth-following
estimate better describes measured whole-column tilt.

The analysed quantity is the shallow-water/topographic proxy
`grad[(f + zeta) / h]`, sampled using the detected eddy ellipse and relative
vorticity at each cached depth. It is **not** local three-dimensional Ertel PV.
`TiltDis` and `TiltDir` are existing whole-column measurements and are not
re-estimated at each depth.

Run on Katana after building the two caches in `depth_following_pv/`. The
notebooks select the cached levels nearest 0, 200, 500, 700 and 1000 m, always
report the exact selected depths, and use matched Eddy-Day samples for formal
cross-depth comparisons.

1. `01_depth_data_coverage.ipynb` — coverage, sample attrition and spine displacement.
2. `02_surface_depth_pv_comparison.ipynb` — paired changes in magnitude, direction and bathymetry.
3. `03_regime_switching_with_depth.ipynb` — planetary/mixed/topographic reclassification.
4. `04_depth_resolved_tilt_direction.ipynb` — polarity-aware directional agreement.
5. `05_depth_resolved_tilt_magnitude.ipynb` — tilt magnitude versus depth-dependent mismatch.
6. `06_best_pv_representation.ipynb` — fixed-depth versus 0–1000 m vector-mean scorecard.
7. `07_topographic_mismatch_case_studies.ipynb` — maps of surface-missed deep topography.

The 2:1 regime threshold is used throughout. AEs are expected along signed
`grad(PV)` and CEs opposite it. Summaries are eddy-equal where repeated daily
observations would otherwise allow long-lived eddies to dominate.
