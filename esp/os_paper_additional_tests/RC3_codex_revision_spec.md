# Codex task: revise ESP RC3 sensitivity experiments

## Objective

Modify the **existing** RC3 sensitivity analysis, rather than creating independent alternative implementations. Produce exactly two concise executable notebooks and at most two supporting Python modules, with one publication-ready results figure from each notebook (we are sticking to striclty just the circular idealised eddy i.e. no elliptical tests). Produce one small supplementary sampling-geometry figure and one background flow figure. The scientific objective is to quantify recovery of known ESP parameters under (1) independent velocity noise at three controlled sampling densities and (2) unmodelled uniform background flow at high sampling density. Do not interpret the results as a universal ranking of SOLO, DOPPIO and LATTE.

## Read before editing

1. Inspect the repository, especially `rc3_noise_sampling.ipynb`, `rc3_background_flow.ipynb`, `rc3_validation.py`, `rc3_plotting.py`, `modules/utils.py`, and the source implementations in `final_ESP_paper_notebooks/solo_idealised_tests_updated.ipynb`, `doppio_idealised_tests_updated.ipynb`, `latte_idealised_tests.ipynb` (adapt paths if necessary).
2. Read the current local README explaining the previous figures and parameter conventions. Preserve the existing method implementations, numerical conventions, optimizer settings, parameter meanings and accepted-fit checks unless a necessary change is explicitly documented.
3. Preserve the manuscript's existing idealised **circular Gaussian eddy** exactly: centre `(0, 0)`, `Rc = 85 km`, `psi0 = 250000 m^2/s`, `Q = I`, velocity/sign conventions and reference-speed calculation. Do **not** alter the eddy strength or truth to improve results. All three methods use this identical truth.
4. Both inner-core and outer-core fitting stages receive **all sampled synthetic observations**, as in the existing idealised tests. Do not introduce the inner/outer observation-selection procedure used for real observational datasets.
5. Reuse the current fitting adapters and **final fitted** parameter `Omega` (the value reported with the final recovered parameters), **not** the preliminary inner-core `Omega_0`. Verify the final output field and units against the existing notebooks before computing error.

## Shared sampling design

Use the same centred geometry and high-density observation locations in both experiments.

- SOLO: east–west transect through the true centre: `x/Rc` equally spaced from `-1.2` to `+1.2`, `y=0`.
- DOPPIO: two orthogonal transects through the true centre, each extending from `-1.2 Rc` to `+1.2 Rc` along its axis. The intersection `(0,0)` is one unique observation, shared by both transects. Preserve the expected DOPPIO input structure.
- LATTE: reproducible scattered positions within a disk of radius `1.2 Rc`, centred at `(0,0)`. Keep the point set fixed across all noise levels and Monte Carlo trials at a given density. Use reproducible nested sets where practical, but do not compromise the required outer-core coverage to achieve nesting.
- For all geometries, verify and record the count of **unique** observations satisfying `rho > Rc/sqrt(2)` (for this circular truth, `rho = hypot(x,y)`). Require LATTE to have **at least 20%** of its points outside that contour and **at least four** such points at LOW density. For transects, ensure the fixed endpoints and equally spaced points provide observations beyond that contour. Fail explicitly if constraints are not met. This is an experimental sampling constraint, **not** a general threshold recommended for ESP.
- Do not vary sampling offset, transect orientation, eddy shape, or spatial coverage in these tests. The manuscript's earlier idealised experiments already examine sampling geometries; this work isolates the stated perturbations.

### Density regimes

| Regime | SOLO (total unique) | DOPPIO (per transect; unique total) | LATTE (total unique) |
|---|---:|---:|---:|
| LOW | 11 | 11; 21 | 21 |
| MEDIUM | 31 | 31; 61 | 61 |
| HIGH | 101 | 101; 201 | 201 |

Do not imply the methods have equal observation budgets or equivalent information content. Report both per-transect and unique DOPPIO counts and outer-core counts in the notebook and saved metadata.

## Experiment 1: noise and sampling

Notebook: `rc3_noise_sampling.ipynb`.

- Evaluate LOW, MEDIUM and HIGH for SOLO, DOPPIO and LATTE using the common circular truth.
- Independent additive Gaussian noise to each velocity component at each **unique** location: `epsilon_u, epsilon_v ~ N(0, sigma^2)`. At DOPPIO's shared intersection, draw noise **once** and reuse it on both transects.
- Levels: `100*sigma/U_ref = [0, 2.5, 5, 10]` (%), where `U_ref` is the maximum true vector speed of the unchanged circular eddy, excluding perturbations.
- Run 100 independent Monte Carlo noise realisations for every nonzero-noise method/density/level. Fit the clean (`0%`) configuration once per method/density. Use a fixed recorded RNG seed (retain existing seed if feasible), paired standard-normal draws across noise amplitudes within a method/density/replicate, and fixed observation locations within a density regime.
- Do not add background current. Do not vary positions between noise replicates. Do not tune the inversion per noise level.
- Calculate errors **against known truth**, using final fitted parameters:
  - centre displacement: `norm(center_hat - center_true)` in km;
  - angular velocity: `100 * abs((Omega_hat - Omega_true) / Omega_true)`;
  - aspect ratio: `100 * abs((alpha_hat - alpha_true) / alpha_true)`;
  - Gaussian radius: `100 * abs((Rc_hat - Rc_true) / Rc_true)`.
- SOLO fixes `alpha=1`; show its aspect-ratio panel as **"Not estimated"**, not zero error. For DOPPIO/LATTE, recovered `alpha > 1` for circular truth represents spurious ellipticity. Orientation is undefined for circular truth; do not score it.
- Plot the **median** error with **10th–90th percentile** bands over accepted trials. The bands show Monte Carlo spread, not confidence intervals. Zero-noise points are deterministic baselines.
- Retain every trial record, including optimizer failures, invalid estimates, and extreme-but-valid estimates. Compute error distributions conditional on accepted fits and report `accepted/attempted` separately. Do not trim, cap, silently omit or replace large valid errors. If failure rates are nonzero, disclose them prominently in the notebook, caption and summary.

### Noise results figure

Exactly **one** main results figure, **3 rows × 4 columns**:
- rows: SOLO, DOPPIO, LATTE;
- columns **in this order**: **Centre**, **final angular velocity Ω**, **aspect ratio α**, **Gaussian radius Rc**;
- horizontal axis: component noise SD / `U_ref` (%);
- colours: LOW, MEDIUM, HIGH, consistent in all panels;
- SOLO aspect-ratio panel: `Not estimated`;
- clear units, readable legend, lettered panels, consistent metric-axis scaling where defensible, and no misleading clipping. Use log/symlog only where justified and label it clearly. Save publication-resolution PNG. No separate acceptance-rate column; give rates in caption/text or supplementary output, especially if below 100%.

## Experiment 2: uniform background flow

Notebook: `rc3_background_flow.ipynb`.

- Reuse **exactly** the HIGH-density geometries and observation positions from Experiment 1. This must be verified programmatically, not merely assumed.
- No random velocity noise; each condition is one deterministic fit.
- Add the **same constant vector** to every sampled velocity: `u_obs(x) = u_eddy(x) + U_bg`. Do not subtract the background or jointly fit a background parameter; the point is to test unmodelled flow.
- Background-speed levels: `100*|U_bg|/U_ref = [0, 5, 10, 20, 30, 50]` (%).
- SOLO: test both **eastward** `(B,0)` (parallel to transect) and **northward** `(0,B)` (perpendicular).
- DOPPIO: test **northward only**.
- LATTE: test **northward only**.
- Include one shared zero-background baseline per method (do not double-count the duplicated SOLO zero level).
- Compute the same four **total absolute errors against truth** as in Experiment 1, using final `Omega`, not `Omega_0`. Do not plot incremental changes from the zero-background fit in the main figure.
- Treat this as sensitivity to **one representative flow direction** for DOPPIO and LATTE. Do not claim those methods are direction-independent. SOLO's two directions are intentionally retained to expose direction-dependent biases.
- Retain invalid fits and large finite estimates in outputs; show or annotate numerical failures and disclose acceptance rates. A successful optimizer termination is not proof of accurate recovery. Do not silently hide large `Rc` errors.

### Background results figure

Exactly **one** main results figure, **1 row × 4 columns**, in the same order: **Centre**, **final angular velocity Ω**, **aspect ratio α**, **Gaussian radius Rc**. Use consistent colours for SOLO/DOPPIO/LATTE, and distinct solid/dashed line styles for SOLO's eastward/northward flows. DOPPIO/LATTE each have northward curves only. Label speed as `|U_bg|/U_ref (%)`. Show total errors against truth, not differences from baseline. Save PNG. If the radius axis spans orders of magnitude, choose a legible justified scale and annotate extreme values rather than silently cropping.

## Supplementary Figures — Sampling Geometries

Create two additional publication-quality supplementary figures illustrating the observation geometries used in the sensitivity experiments. These figures must be generated from the actual sampling coordinates used by the experiments, not independently constructed schematic coordinates.

### Supplementary Figure S1 — Noise and Sampling Experiment

**Output:** `results/rc3_noise_sampling_geometry.pdf` and `.png`

- Create a 3 × 3 figure.
- Rows: SOLO, DOPPIO, LATTE.
- Columns: LOW, MEDIUM, HIGH sampling density.
- Plot all observation locations as scatter points.
- Overlay the true eddy centre at (0, 0), marked with a black cross.
- Overlay two concentric circles:
  - Dashed circle: maximum tangential-velocity contour, \(r=R_c/\sqrt{2}\).
  - Solid circle: Gaussian characteristic radius, \(r=R_c\).
- Use identical axis limits and equal aspect ratios across all nine panels.
- Label panels (a)–(i).
- Annotate each panel with its number of **unique velocity observations**, \(N\), and the number beyond the maximum tangential-velocity contour, \(N_{\mathrm{outer}}\).
- For DOPPIO, count the shared transect intersection only once.
- For LATTE, show the actual reproducible scattered observation locations.
- Include a concise legend identifying the eddy centre and the two reference contours.

### Supplementary Figure S2 — Background-Flow Experiment

**Output:** `results/rc3_background_flow_geometry.pdf` and `.png`

- Create a 1 × 3 figure.
- Panels: (a) SOLO, (b) DOPPIO, (c) LATTE.
- Use the exact HIGH-density sampling coordinates from Supplementary Figure S1.
- Show observation locations, true eddy centre, \(r=R_c/\sqrt{2}\), and \(r=R_c\), using identical formatting to S1.
- Overlay arrows indicating the directions of the imposed uniform background currents:
  - SOLO: eastward and northward.
  - DOPPIO: northward.
  - LATTE: northward.
- Use arrows of equal display length to represent **direction only**, not background-current magnitude. State this explicitly in the caption.
- Position arrows so they do not obscure observation points or reference contours.
- Annotate each panel with \(N\) and \(N_{\mathrm{outer}}\).
- Use identical spatial axis limits and equal aspect ratios across all three panels.

### Implementation Requirements

1. Add reusable plotting functions to the existing `rc3_plotting.py` module. Do not duplicate geometry-generation logic.
2. Reuse the exact sampling arrays from the two experiments. Assert that the background experiment's coordinates match the HIGH-density noise experiment's coordinates.
3. Calculate \(N_{\mathrm{outer}}\) directly from the coordinates and the true eddy centre.
4. Do not regenerate random LATTE locations when creating the figures.
5. Export both figures at publication quality (300 dpi PNG).
6. Keep the figures supplementary: they must not change the primary sensitivity experiments or their numerical results.
7. Include descriptive captions explaining the sampling configurations, reference contours and background-current directions.

## Code and deliverables

- Keep only the two requested notebooks as primary analysis entry points. Put all nontrivial reusable simulation, sampling, fitting-adapter, aggregation and plotting code in **one or two** Python modules (prefer editing existing `rc3_validation.py` and `rc3_plotting.py`). Do not duplicate model or fitting implementations.
- Put all configurations used for each test within a config file. This includes the shared idealised eddy parameters. 
- Notebook structure: run experiment; print compact sampling/acceptance/summary tables; produce **one** main paper figure; save outputs. Supplementary sampling diagram  generated by a shared helper but should not become another main notebook results figure.
- Write machine-readable trial-level CSV, summary CSV, and provenance JSON for **each** experiment. Include seed, truth parameters, `U_ref`, noise/background levels, sampling geometry, unique counts, outer-core counts, optimizer/acceptance settings, definitions/units, code provenance if accessible, and figure paths.
- Ensure deterministic reruns. Use non-destructive filenames for revised outputs; avoid overwriting earlier results until verified.
- Do not introduce arbitrary new rejection rules based on parameter error, goodness of fit or proximity to truth. Preserve and report numerical validity checks.
- Include a small self-test/assertion suite (in the helper module or a minimal test file only if needed): geometry extents, outer-core counts, shared DOPPIO intersection, matching HIGH geometry between notebooks, correct error definitions, final `Omega` field, consistent `U_ref`, and finite clean-baseline results.
- Do not change the underlying idealised eddy, source fitting methods, observational workflows, or manuscript prose. If a source method cannot operate on the centred geometry, investigate and report the reason instead of silently modifying the method or moving the transect.


## Interpretation guardrails

These are controlled sensitivity tests with an exactly specified Gaussian truth. Noise is independent between observations; correlated errors, evolving eddies, shear, non-Gaussian structures and platform-specific observational errors are **not** assessed. A uniform background current is a specific model mismatch, not random noise. Sampling density is varied at fixed geometry and spatial extent; this is not a comprehensive survey-geometry study. Preserve any baseline bias instead of forcing clean fits to have zero error. The results should support qualified statements about when parameter recovery degrades, not blanket claims of threshold-free operation or general robustness.
