# Codex implementation brief — RC3 sensitivity experiments for ESP

## Mission

Implement **two independent, reproducible numerical experiments** supporting the revised Elliptical Streamfunction Parameterisation (ESP) methods paper:

1. **Noise and sampling sensitivity:** quantify how accurately ESP recovers *known* eddy parameters under degraded velocity sampling and measurement noise.
2. **Background-flow sensitivity:** quantify biases in recovered parameters when a spatially uniform background velocity is added to a known eddy and not explicitly represented by the fitted ESP model.

Deliver **exactly two concise Jupyter notebooks**, one for each experiment, and **one or two reusable `.py` modules** containing all simulation, inversion-adapter, metrics, and plotting helpers. Each notebook must produce **exactly one publication-ready multi-panel figure** (one `Figure` object; subplots are allowed). Save the figures to disk. Avoid unrelated analyses, diagnostic figures, or a large framework.

These experiments are controlled *sensitivity tests*, not a comprehensive observational uncertainty study. Report precisely what the results demonstrate and what they do not.

## Repository-first instructions (mandatory)

Existing ESP repository, which was used for preivous anlaysis of the paper, can be found here /srv/scratch/z5297792/UNSW-MRes/MRes/MOCA_PAPER/final_ESP_paper_notebooks with idealiesed tests for SOLO (solo_idealised_tests_updated.ipynb), DOPPIO (doppio_idealised_tests_updated.ipynb), LATTE (latte_idealised_tests.ipynb) and also real applications of the methods in [elliptical_streamfunction_parametrisation](https://github.com/regdowse?tab=repositories#:~:text=elliptical_streamfunction_parametrisation).

1. Inspect the repository before writing code. Identify the canonical elliptical Gaussian streamfunction, velocity evaluation, parameter representation, and existing SOLO, DOPPIO, and LATTE inversion implementations. Search for existing idealised experiments, test cases, geometry builders, and plot styling.
2. **Reuse existing inversion code and established parameter conventions.** Do not reimplement the inversion algorithms or silently modify them to improve the results. If an adapter is necessary, isolate it in a helper module and document any conversion.
3. Read the current manuscript / method documentation in the repository, if available, to establish the definition of `Rc`, `Q`, aspect ratio, orientation, vorticity, and coordinate signs. Do not infer definitions from variable names alone.
4. Keep changes limited to the two notebooks, one or two new `.py` files, and their generated figures/results. Do not modify the core ESP library unless a genuine bug prevents the experiments; report the issue first.
5. Before implementing, briefly document in a notebook Markdown cell which existing functions are reused, their module paths, and any assumptions necessary to run the tests.
6. If a required interface is unavailable or ambiguous, stop and ask a focused question instead of inventing an ESP implementation.

## Target files

Use existing repository naming/layout conventions where possible; otherwise create:

```text
rc3_noise_sampling.ipynb
rc3_background_flow.ipynb
rc3_validation.py        # synthetic fields, geometry, inversion adapters, metrics
rc3_plotting.py          # optional; only if it improves clarity
```
Produce the plots in the notebooks. Do not commit large intermediate arrays. 

## Common scientific setup

### Truth model

Use the manuscript's existing elliptical Gaussian ESP velocity model with a known parameter vector. In generic notation,

\[
\psi(\mathbf{x})=\psi_0\exp\left[-\frac{(\mathbf{x}-\mathbf{x}_0)^TQ(\mathbf{x}-\mathbf{x}_0)}{R_c^2}\right],
\quad Q=Q^T\succ0,\quad\det Q=1.
\]

Obtain velocities from the existing library and respect its rotation/sign convention. Choose a representative, noncircular eddy based on an existing manuscript idealised case, preferably its published baseline (e.g., `Rc ≈ 85 km` **only if that is the actual reference case**). Use the same truth in both experiments. Explicitly record truth parameters and their units.

Use a noncircular eddy so orientation and aspect-ratio recovery are meaningful. Do not report orientation error for nearly circular fitted eddies without a clear convention.

### Inversion methods and sampling

Prefer testing **SOLO, DOPPIO, and LATTE**, using their existing geometries and fitting functions. Where a method is mathematically underdetermined for a chosen sampling geometry, mark it as unavailable rather than returning a misleading estimate. The notebook should remain short by delegating all experiment loops to the helper module.

Sampling density and spatial coverage are **different variables**. Hold geometry and coverage fixed when varying the number of samples; do not conflate fewer observations with shorter transects. Use the same sampling geometry for clean/noisy comparisons. Keep inversion initialisation, bounds, search radii, stopping conditions, and other tuning choices **fixed across all experimental conditions**, unless a method requires an explicit documented exception. Do not tune on the evaluation outcomes.

### Metrics

Compute parameter-recovery errors against known truth, including at least:

- Centre position error, `||x0_hat - x0_true||` in km.
- Core radius error, `|Rc_hat - Rc_true| / Rc_true` in percent.
- Aspect-ratio error, `|alpha_hat - alpha_true| / alpha_true` in percent, if alpha is part of the existing fit.
- Optionally orientation error in degrees (axial/circular distance modulo 180°), only if reliably identifiable.

Record convergence/failure rate separately. **Do not discard failed fits silently**; report success counts and use `NaN` only with explicit failure status. In-sample `R²` may be retained as a diagnostic but must not replace truth-based errors. If methods estimate different subsets of parameters, show only comparable metrics and explain omissions.

Use a fixed random seed (e.g. `20261027`) and `numpy.random.default_rng`. Keep all physical and numerical quantities in consistent units, with conversion at input/output boundaries.

## Experiment 1 — noise and sampling sensitivity

### Scientific question

How do known-parameter recovery errors and inversion success change with velocity noise amplitude and observation count for the supported ESP sampling configurations?

### Experimental design

1. Generate noise-free velocity observations from the reference ESP eddy at existing method-specific sampling geometries.
2. Define a modest grid of observation counts, e.g. `N = [10, 20, 50, 100]`, **subject to each method's minimum required sample count**. For SOLO/DOPPIO, specify whether `N` is per transect or total and keep that definition consistent. LATTE should use a fixed spatial domain and reproducible scattered-point sampling.
3. Define dimensionless noise levels relative to a clearly stated reference eddy speed, e.g. `sigma / U_ref = [0, 0.025, 0.05, 0.10]`. Choose `U_ref` from the true eddy field (for example the maximum model speed over the full reference domain), not separately from each degraded sample. **Label these as controlled perturbations, not measured observational error statistics.**
4. Add zero-mean Gaussian perturbations to observed velocity components. Make component covariance assumptions explicit. For the primary test, independent errors are sufficient; if easy to implement correctly, add one **clearly distinguished** spatially correlated error condition, using a positive-semidefinite covariance kernel and fixed correlation length, without increasing the figure's complexity excessively. Do not claim correlated-error robustness unless actually tested.
5. For each method, `(N, noise)` condition and random replicate, generate observations, fit with the existing method, and calculate recovery errors and success/failure. Use approximately 100 replicates per stochastic condition if computationally feasible. A noise-free baseline need not be repeated unnecessarily. Precompute or cache static geometry where useful.
6. If fitting all three methods at full scale is prohibitively expensive, prioritise a scientifically coherent subset and document the reason; **never silently reduce the study**. The preferred outcome is all three methods.
7. Keep sampling coverage constant across `N`. Use common random seeds / deterministic random-number management so reruns reproduce the figure exactly. Avoid conditioning the analysis on successful fits without displaying the corresponding success fraction.

### One paper figure

Create one clear multi-panel figure showing recovery error against noise amplitude, with curves or grouped points for `N`, faceted by inversion method and/or metric. Prefer median with 10th–90th percentile (or 5th–95th), consistently labelled; include sample counts and failed-fit information in the legend, a compact annotation, or a separate panel **within the same figure**. If a three-method × three-metric grid is too dense, select the two most informative metrics for the main figure and save the remaining metrics to a compact dataframe (saved in the notebooks memory).

The figure must support a restrained conclusion about tested sampling densities and perturbation amplitudes. Avoid a plot dominated by goodness-of-fit scores.

## Experiment 2 — uniform background-flow sensitivity

### Scientific question

How does a **spatially uniform, unmodelled** background velocity bias ESP parameter estimates, and does the sensitivity differ between SOLO, DOPPIO, and LATTE?

### Experimental design

1. Use the **same truth eddy and sampling geometries** as Experiment 1, with a fixed adequately resolved observation count (e.g. `N=100` if appropriate).
2. Construct synthetic observations as

   \[
   \mathbf{u}_{\rm obs}(\mathbf{x})=\mathbf{u}_{\rm ESP}(\mathbf{x})+\mathbf{U}_{\rm bg},
   \]

   with a constant vector `U_bg` across all sample locations. The fit must use the **unmodified existing ESP inversion** unless that inversion already includes a background term by default; if so, explicitly identify this and provide a defensible no-background-fit comparison. Do not subtract the imposed background before fitting the primary experiment.
3. Vary background speed as a fraction of the same `U_ref`, e.g. `|U_bg| / U_ref = [0, 0.05, 0.10, 0.20, 0.30, 0.50]`.
4. Test **background direction**, not just magnitude: at minimum eastward and northward vectors in the model's Cartesian coordinates. These directions matter because SOLO transect orientation can cause strongly anisotropic bias.
5. Fit each method at each speed and direction. Record centre displacement vector and magnitude, relative core-radius error, aspect-ratio error, convergence/failure, and any other metric justified by the existing implementation.
6. This is a **deterministic sensitivity test**; do not add random noise to the primary analysis. No Monte Carlo replication is required unless a method itself uses stochastic sampling. For LATTE, fix one reproducible point cloud across background conditions.
7. Include the zero-background baseline and check that baseline errors are small enough to interpret incremental biases. Where baseline bias is nonzero, distinguish **total error relative to truth** from **incremental change relative to zero background**.
8. Do not assume DOPPIO or LATTE are immune to background flow. Measure the bias and describe only what the results support. A uniform flow is not the same as background shear: **do not claim to have tested shear**.

### One paper figure

Create one multi-panel figure showing parameter-recovery error versus `|U_bg| / U_ref`, separating methods and distinguishing background direction by line style or colour. Prioritise **centre-position error (km)** and **core-radius relative error (%)**; aspect ratio can be a third panel or saved as a table if needed. Include the zero-background baseline. Use consistent scales and readable legends. A direction-sensitive centre-displacement plot is an acceptable alternative if it more clearly reveals SOLO biases, but the figure must show the dependence on imposed background magnitude.

## Notebook structure (both notebooks)

Keep each notebook short, approximately 4–7 cells:

1. **Markdown:** purpose, precise scientific question, reused library functions and assumptions.
2. **Imports/config:** local module imports, immutable experiment settings, seed, output path.
3. **Run:** call one high-level helper function returning tidy `pandas.DataFrame` results.
4. **Summary:** print a compact parameter/trial/success summary and optionally display a small table. No extra figures.
5. **Figure:** call one plotting helper, create **one** `matplotlib.figure.Figure` and display it in the notebook.
6. **Markdown:** 2–4 factual interpretation bullets automatically grounded in the results where feasible, plus a one-sentence limitation. Avoid hard-coded claims before running the experiment.

Use standard NumPy, pandas and Matplotlib and the repository's existing numerical dependencies. Do not introduce heavy dependencies without need.

## Quality control and acceptance criteria

- Both notebooks execute from a clean kernel in repository context, with no hidden state or manual edits.
- Exactly **two** notebooks, **one or two** supporting `.py` modules, and **one** figure per notebook.
- Each figure is publication-ready: Labelled units, panel labels `(a)`, `(b)` etc., consistent terminology, no titles unless needed, no excessive decimals. No need to make any adjustements to things like fontsizes, this is something I can tune myself. 
- Baseline (`sigma=0`, `U_bg=0`) recovers the reference parameters to expected numerical accuracy for an adequately sampled, identifiable configuration; if not, investigate and **report** rather than hide it.
- `Q` remains symmetric positive definite with determinant 1 when applicable. Preserve the established definition of aspect ratio and `Rc`.
- Report inversion failures, number of successful trials, and randomness/covariance assumptions. Do not fabricate uncertainty bars from deterministic conditions.
- No fitted-model `R²` presented as independent validation. Avoid universal claims about real-data robustness, arbitrary sampling, background shear, or 3D reconstruction.

## Suggested public helper API (adapt to existing code)

```python
run_noise_sampling(config) -> pd.DataFrame
run_background_flow(config) -> pd.DataFrame
plot_noise_sampling(results, config) -> matplotlib.figure.Figure
plot_background_flow(results, config) -> matplotlib.figure.Figure
```

Results should be **tidy**, one row per method / condition / replicate, with stable fields such as `method`, `n_samples`, `noise_fraction`, `background_fraction`, `background_direction`, `replicate`, `success`, `centre_error_km`, `rc_error_pct`, `alpha_error_pct`, and optional `orientation_error_deg`. Keep the original fitted parameters as additional columns where practical so errors can be audited.

## Final execution order

1. Locate and verify existing ESP model/inversion APIs.
2. Implement the smallest reusable adapter and metrics layer.
3. Run zero-perturbation sanity checks for each method.
4. Implement and smoke-test noise/sampling experiment.
5. Implement and smoke-test background-flow experiment.
6. Run full experiments and generate both final figures.
7. Inspect figures for interpretability and verify files exist.
8. Summarise outcomes and limitations without modifying manuscript text.
