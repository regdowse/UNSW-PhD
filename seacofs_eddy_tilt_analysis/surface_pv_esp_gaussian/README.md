# ESP-Gaussian surface PV-gradient reconstruction

This visual-first series showcases a single-vector surface PV-gradient method
that uses the complete fitted ESP parameterisation. For every Eddy-Day and
grid cell,

`rho**2 = [x-xc, y-yc] Q [x-xc, y-yc].T`

and

`zeta(x,y) = w * exp(-rho**2 / Rc**2)`.

The same dimensionless Gaussian shape supplies smooth spatial weights. The
primary one-row-per-snapshot vector remains the environmental gradient:

`grad(PV)_environment = grad(f)/h - (f + zeta) grad(h)/h**2`.

The internal `grad(zeta)/h` contribution and complete reconstructed PV
gradient are stored separately as `PV_grad_eddy_*` and `PV_grad_full_*`. They
are not silently mixed into the steering/topographic vector.

`surface_method="esp_gaussian"` is an option in `add_pv_gradient_terms`, not a
separate public calculation. This preserves the established API and output
schema while keeping the reconstruction in a tested private helper. The
depth-following pathway remains unchanged.

Run `00_build_esp_gaussian_cache.ipynb` first on Katana, then notebooks 01–06:

1. Gaussian reconstruction and geometry.
2. Truncation/support convergence.
3. Uniform versus Gaussian-weighted environmental gradients.
4. Environmental versus full reconstructed PV gradients.
5. Seamount encounter examples through time.
6. Tilt relationships and final method assessment.

The Gaussian is evaluated to ellipse scale 2 for the primary candidate. At
that boundary `rho**2/Rc**2 = 2`, so its weight is approximately 0.135. Runs
at scales 1 and 1.5 establish whether truncating the remaining tail affects
the result materially.
