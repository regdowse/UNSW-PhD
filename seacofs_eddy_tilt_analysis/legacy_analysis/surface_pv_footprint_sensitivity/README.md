# Surface PV-gradient footprint sensitivity

This visual-first series determines how the surface-centred shallow-water PV
gradient depends on the area sampled around each eddy. It deliberately leaves
the depth-following calculation unchanged.

The corrected default calculation averages the completed nonlinear local
gradient components. Existing `PV_grad_*_mag` values remain magnitudes of the
net spatially averaged vectors. New diagnostics retain strong but cancelling
topographic gradients:

- `PV_grad_topo_mean_local_mag` — mean local magnitude;
- `PV_grad_topo_rms_local_mag` — RMS local magnitude;
- `PV_grad_topo_p90_local_mag` — upper-tail local magnitude;
- `PV_grad_topo_coherence` — net magnitude / mean local magnitude.

Low coherence with high local magnitude describes a strong, heterogeneous
feature such as a seamount without falsely assigning it one coherent forcing
direction. `averaging="legacy"` preserves the historical standard columns.

Run `00_build_surface_footprint_cache.ipynb` first on Katana, followed by:

1. `01_legacy_vs_nonlinear.ipynb`
2. `02_fraction_sensitivity.ipynb`
3. `03_annulus_sensitivity.ipynb`
4. `04_seamount_detection.ipynb`
5. `05_tilt_relationship_stability.ipynb`
6. `06_footprint_recommendation.ipynb`

The prespecified filled-ellipse scales are 0.25, 0.5, 0.75, 1.0 and 1.25.
Because `frac` is a linear scale, the corresponding areas scale as `frac**2`.
Annuli isolate 0.25–0.5, 0.5–0.75, 0.75–1.0 and 0.5–1.0 of the original
ellipse. Do not choose a footprint only because it maximises a tilt
relationship; prefer stability, feature detection and physical interpretability.
