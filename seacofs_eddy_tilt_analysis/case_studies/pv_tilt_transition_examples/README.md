# Selected eddy PV-gradient and tilt transitions

This visual-first case-study series brings together the existing long-eddy,
tilt, bathymetry and fixed-core ESP-Gaussian PV-gradient analyses for the
following user-selected eddies:

- AE356, AE484, AE604, AE739, AE1015, AE1260, AE1268 and AE2255.
- CE1105, CE1194, CE1552 and CE1927.

Run the notebooks on Katana in order:

1. `00_build_selected_eddy_cache.ipynb` calculates the surface PV-gradient
   terms only for the selected tracks using `FRAC=1`, Gaussian vorticity and
   Gaussian spatial weighting.
2. `01_eddy_overviews.ipynb` creates one compact time-series-and-map overview
   per eddy. Green, orange and grey denote planetary, topographic and mixed
   regimes, respectively. Blue arrows show tilt and magenta arrows show the
   measured environmental PV-gradient direction.
3. `02_selected_day_core_maps.ipynb` suggests planetary, transition and
   topographic days, lets the user override them, and plots the local
   fixed-core gradients, their Gaussian-weighted mean and the depth-resolved
   eddy spine.

The overview uses `PV_grad_mag` and `PV_grad_theta` as the primary
environmental measure. Low-coherence directions are drawn with reduced
opacity because they are residual directions formed by cancellation among
strong local gradients.
