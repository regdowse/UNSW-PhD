# Observed cyclones versus model climatology

Run `obs_clim_check.ipynb` on Katana in the usual scientific Python environment.
It reads existing processed, saved-tilt and confirmed-profile tables plus the
reference grid and z_r; no production cache regeneration is required. It does
not change those files. Defaults and paper measurement provenance are in the
notebook. `obs_clim_tools.py` contains only the reusable fitting/matching/reporting
operations. The original measurement notes remain at the end of the notebook.

The main comparison uses daily upper-column OLS tilt and peak signed relative-
vorticity Rossby number over 0–600 m. Observed paper Ro is converted using its
stated solid-body velocity relation, zeta=2*Omega. Actual model fit depths are
reported and must reach at least 600 m. Paper SW bearings are deep-to-shallow
(see Figure 10); model daily centres rotate with grid.angle. Saved production
delta tilt is displayed separately, already in its geographic convention.

Matching is CE, within 1.5 observed radii, ±25% peak Ro. Samples are weighted by
eddy-day, with track counts and selection attrition. No tests of statistical
significance or confidence intervals are implied. Optional PV values use the
shared point-proxy implementation without reading/writing the shared PV cache.
This does not reproduce the Gaussian core-mean workflow.

Set SAVE_RESULTS=True for a timestamped local output directory containing
observations, membership, profile metrics, tables, the main figure and input
path/size/mtime provenance. Input metadata are not full data checksums.

Local checks (from repository root):

```bash
MPLBACKEND=Agg python -m unittest discover -s seacofs_eddy_tilt_analysis/real_obs -p 'test_*.py' -v
```

Production climatological results require a Katana run. Local synthetic tests
verify signs, physical conversion, depth support, peak selection, key joins,
circular differences and empty matches; they do not validate real results.
