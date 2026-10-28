# Context for the revised circular ESP sensitivity experiments

This document provides experimental context for interpreting two main results figures and one supplementary geometry figure. It describes what the figures show, not their observed trends, recovered errors or success rates. The current revision follows `RC3_LATTE_symmetric_sampling_codex.md` with the subsequent clarification requiring exact 90° rotational symmetry as well as both coordinate-axis reflections. The preceding asymmetric LATTE run and its stale background-geometry figure are preserved in `archive/before_latte_d4_sampling/`; current main figures and numerical outputs are in **`results/revised/`**.

## Files and configuration

Exactly two notebooks are the analysis entry points:

| Entry point | Main results figure | Saved supplementary figure |
| --- | --- | --- |
| `rc3_noise_sampling.ipynb` | `results/revised/rc3_noise_sampling.png` | S1: `results/rc3_noise_sampling_geometry.png` |
| `rc3_background_flow.ipynb` | `results/revised/rc3_background_flow.png` | Shares S1; no separate background geometry figure |

Each notebook displays one main results figure. The noise notebook saves the sole geometry supplement S1 separately; the background notebook reuses its HIGH configurations and does not generate another supplement. Figures are exported as 300 dpi PNG only. All test configurations—including truth, sampling counts, seeds, noise/current levels, original numerical settings, output directories, colours and scale conventions—are in **`rc3_config.json`**. Restart the kernel before executing either notebook. Both work from this directory or the UNSW-PhD repository root.

The shared helpers are `rc3_validation.py` and `rc3_plotting.py`. No alternate inversion implementation is introduced. The adapter selects the unchanged fitting functions from `solo_idealised_tests_updated.ipynb`, `doppio_idealised_tests_updated.ipynb` and `latte_idealised_tests.ipynb` under the configured source root's `MOCA_PAPER/final_ESP_paper_notebooks/`, and reuses `modules/utils.py` for Gaussian velocities and projection/shape utilities.

## Common known truth and parameter conventions

All three methods receive the **identical stationary circular Gaussian eddy**:

\[
\psi(\mathbf{x})=\psi_0\exp(-\rho^2/R_c^2),\qquad
\rho^2=(\mathbf{x}-\mathbf{x}_0)^TQ(\mathbf{x}-\mathbf{x}_0).
\]

The true centre is `(0,0)` m, `Rc = 85000 m`, `psi0 = 250000 m²/s`, and `Q = I`. Thus `alpha = 1`, where alpha is major/minor semiaxis length. Orientation is undefined and is not scored. There are no elliptical-truth tests in this revision. DOPPIO and LATTE still estimate unrestricted elliptical shapes; an estimated alpha greater than one represents spurious ellipticity. SOLO prescribes alpha=1 and its shape error is **not estimated**, rather than a measured zero error.

The actual source velocity code uses **u = −∂ψ/∂y and v = ∂ψ/∂x**, with x eastward and y northward. Its signed amplitude is `Omega = −2 psi0/Rc²` in s⁻¹. Centre vorticity is `Omega*trace(Q)`; Omega and vorticity are different quantities. The derivative signs stated in the previous README were reversed; this documentation correction does not alter the original numerical model. An independent finite-difference self-test checks these signs.

Rc is the Gaussian characteristic scale. The true maximum tangential-speed contour is **r = Rc/√2**, not r=Rc. The circular truth has rho=r. The fixed reference speed `U_ref` is the maximum **true vector-speed magnitude**, evaluated with the original velocity function at its analytic peak. It is common to all methods and conditions and excludes imposed noise/current. It is not recomputed from degraded sample sets or fitted velocities.

Public fields use SI units. The adapter converts coordinates to the original notebooks' km/m·s⁻¹ fitting convention and converts fitted parameters back. **`omega_s` is the final outer-fit Omega in s⁻¹ and is the quantity scored.** `omega_inner_s` is the preliminary inner-core amplitude, retained only for audit.

## Centred sampling and density regimes

All positions are fixed within each regime and shared across perturbation levels and replicates. Spatial extent does not change with density.

| Regime | SOLO unique total | DOPPIO per transect | DOPPIO unique total | LATTE unique total |
| --- | ---: | ---: | ---: | ---: |
| LOW | 11 | 11 | 21 | 21 |
| MEDIUM | 31 | 31 | 61 | 61 |
| HIGH | 101 | 101 | 201 | 201 |

- **SOLO:** equally spaced east–west transect through the true centre, x/Rc from −1.2 to +1.2, y=0.
- **DOPPIO:** two equally spaced orthogonal transects through the true centre, each from −1.2 Rc to +1.2 Rc. The intersection is one unique vector observation reused by both transects, including the same noisy u/v draw.
- **LATTE:** a stratified, D4-symmetric random point set strictly inside a disk of radius 1.2 Rc. There is exactly one centre point. Generic base points are sampled with an area-uniform radial law and angle uniform strictly between 0 and π/4; all eight sign/swap counterparts are included. Four-point groups on the coordinate axes use random radii with the same area-uniform radial law in their stratum. Axis groups are constrained locations, not independent uniform 2D disk observations. There are no boundary anchors. All positions are unique and invariant under 90° rotations and reflections across both axes. Inner and outer strata are separated strictly by r=Rc/√2. No point lies on that partition or the disk's outer boundary. The 1.2 Rc extent describes the fixed sampling domain, not an anchor radius.

The exact LATTE composition is:

| Regime | Inner eight-point groups | Inner axis groups | Outer eight-point groups | Outer axis groups | Inner points including centre | Outer points |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| LOW | 1 | 0 | 1 | 1 | 9 | 12 |
| MEDIUM | 4 | 0 | 3 | 1 | 33 | 28 |
| HIGH | 14 | 1 | 10 | 1 | 117 | 84 |

The fixed seed is 20261027 with geometry stream key 17 and separate inner/outer substreams. Seed pools are shared across densities; lower-density point sets are nested subsets, not necessarily prefixes of the complete coordinate array. Ordering is centre, inner eight-point groups, inner axis groups, outer eight-point groups, outer axis groups. The stored order is reused elementwise in the background experiment. No geometry is selected by fitted errors.

`N_outer_true` counts **unique observed positions with true rho > Rc/√2**, computed directly from coordinates. LATTE outer counts are exactly 12, 28 and 84 at LOW/MEDIUM/HIGH, equal to DOPPIO outer counts. LATTE inner counts, including the centre, are exactly 9, 33 and 117. This matches one aspect of radial coverage, not the full geometry or information content. Transects must include observations beyond it. Geometry validation fails explicitly if these requirements are not met. These are experimental design constraints, not general recommended ESP thresholds.

The unique counts and DOPPIO per-transect counts are printed and saved. They are not equal observation budgets or equivalent information content. The manuscript's earlier geometry tests, rather than this experiment, address offsets, orientation and other sampling configurations.

Both fitting stages receive **all** synthetic observations. No observational-data inner/outer selection workflow is added. The source's internal finite-value and zero-projection filtering remains. In the CSV, `n_outer_fit` counts points eligible for that original outer fit; it is distinct from the true-coordinate coverage count `n_outer_true`. The original SOLO `core_thresh=30000` is retained in backend km units and covers the complete transect. Original data-derived initialisation, stopping rules and post-fit radius caps are preserved and recorded, not tuned by perturbation level.

### Exact coordinate reuse

One immutable, sampling-signature-named coordinate manifest under `results/revised/` supplies both experiments and all supplementary plots. Plotting does not generate new LATTE positions. Background execution asserts array equality against that manifest and, when available, the saved noise experiment's coordinate provenance. Running the noise notebook first provides both checks. A standalone background run can use the shared manifest. A conflicting saved revised noise geometry raises an error instead of being silently reused.

## Four total-error metrics and numerical acceptance

All plotted errors are against known truth:

| Metric | Definition | Units |
| --- | --- | --- |
| Centre | norm(fitted centre − true centre) | km |
| Final angular velocity | 100 × abs((final fitted Omega − true Omega)/true Omega) | % |
| Aspect ratio | 100 × abs((fitted alpha − true alpha)/true alpha) | %; unavailable for SOLO |
| Gaussian radius | 100 × abs((fitted Rc − true Rc)/true Rc) | % |

The three parameter errors are absolute, so they do not indicate the sign of over/underestimation. Centre error is a displacement magnitude, not its direction. No main panel plots a change relative to a fitted zero-perturbation baseline. In-sample R² is not a validation metric in these figures.

Accepted fits require successful source optimiser completion, finite estimates, positive radius, and symmetric positive-definite Q with determinant one, without a source radius-cap fallback. These are numerical validity checks, **not accuracy or identifiability guarantees**. There is no rejection threshold on parameter error, goodness of fit, distance from truth, or closeness to the clean fit. Extreme but valid estimates remain in trial outputs and summaries. Rejected trials retain statuses, warning messages and attempted parameters where available. Conditional error summaries must be read alongside accepted/attempted counts.

The cubic/affine inner-core estimators approximate the Gaussian field. Clean fits may have baseline bias, which is retained rather than subtracted or forced to zero. The self-tests require finite accepted clean HIGH fits, not exact agreement with truth.

## Main figure: independent noise and sampling

**Image:** `results/revised/rc3_noise_sampling.png`.

The observations are the true eddy velocities plus independent zero-mean Gaussian errors in **each component at each unique location**, with variance sigma². The two components are independent; DOPPIO's shared intersection is reused rather than independently perturbed twice. LATTE velocity perturbations are drawn independently at every unique location, including the centre, axis-group points and every reflected/rotated point; the noise vectors are not transformed or symmetrised. There is no background current or randomisation of observation positions.

Horizontal axis: **100 sigma/U_ref = 0, 2.5, 5, 10%**. Each nonzero method/density/level has 100 independent Monte Carlo realisations. Each zero-noise method/density baseline is fitted once. Seed 20261027 is retained. Standard-normal draws are paired across noise amplitudes within a method/density/replicate; physical noise amplitude changes by the specified fraction of the common U_ref.

The layout is **3 method rows × 4 error columns**, in this exact order:

| Row | Centre | Final Omega | Alpha | Rc |
| --- | --- | --- | --- | --- |
| SOLO | (a) | (b) | (c): Not estimated | (d) |
| DOPPIO | (e) | (f) | (g) | (h) |
| LATTE | (i) | (j) | (k) | (l) |

Colours encode LOW, MEDIUM and HIGH consistently across panels. Lines show the median of individual accepted-fit errors. Shading spans the 10th–90th percentiles of those errors: **Monte Carlo spread, not confidence intervals**. A deterministic clean baseline has no estimated replicate spread. Quantiles use every accepted trial, including extreme errors; no error-based trimming is performed. Connecting lines do not imply simulations between tested levels.

There is no acceptance-rate column. Aggregate accepted/attempted counts are in the figure caption; any failures and the lowest condition acceptance are disclosed there, and full per-condition counts/statuses are in notebook tables and saved CSVs. If an entire condition has no accepted error estimate, an x is placed at a fixed annotation height rather than presented as an error value. Corresponding metric columns share vertical scales across methods. Symlog axes are used only when the plotted range warrants them; axis labels and captions identify the scale and its linear threshold. All actual quantile bounds determine the axis limits.

## Main figure: unmodelled uniform background flow

**Image:** `results/revised/rc3_background_flow.png`.

Use exactly the HIGH coordinates above. Add the **same constant vector** to every observed velocity, without subtraction or a jointly fitted background parameter. There is no random velocity noise or Monte Carlo replication. Horizontal axis: **100 |U_bg|/U_ref = 0, 5, 10, 20, 30, 50%**.

SOLO tests eastward (parallel to its transect) and northward (perpendicular). DOPPIO and LATTE test **northward only**, one representative direction. Each method has one shared zero-current baseline; the two SOLO plotted series use that same record rather than counting it twice.

The layout is **1 row × 4 columns**: (a) centre, (b) final Omega, (c) alpha, (d) Rc. All are total absolute errors against truth. SOLO has no alpha curve. Method colours are consistent across panels. Solid/dashed styles distinguish eastward/northward; DOPPIO/LATTE have only northward series. There are no uncertainty bands. Caption and tables disclose acceptance; failed conditions are marked with an x at an annotation height and omitted as error values.

Large valid errors are not cropped or discarded. If a metric spans the configured logarithmic-range criterion, its axis uses a clearly labelled symlog scale with headroom above the plotted curves; maximum-value annotations are disabled. Inspect numeric ticks rather than comparing slopes across differently scaled axes.

## Supplement S1: density-specific sampling geometry

**Image:** `results/rc3_noise_sampling_geometry.png`.

A **3×3** layout has SOLO/DOPPIO/LATTE rows and LOW/MEDIUM/HIGH columns. Points are the actual experiment locations; all panels share limits and equal spatial aspect. The black cross is the true centre, the dashed circle is r=Rc/√2, and the solid circle is r=Rc. Each panel reports its unique N and N_outer. DOPPIO's intersection is counted once. HIGH panels (c), (f) and (i) are reused in the uniform-background-flow experiment: SOLO eastward and northward; DOPPIO/LATTE northward only. There are no background arrows or extra panels. These are sampling diagrams, not additional recovery experiments or inferred eddy boundaries.

## Outputs, verification and interpretation scope

Each experiment saves trial-level CSV, summary CSV, sampling-count CSV and provenance JSON under `results/revised/`. Trial records contain the original fitted parameters, final/preliminary Omega, acceptance/status fields, geometry identifiers and count definitions. Noise summaries contain medians and 10th/90th percentiles; deterministic background summaries do not invent uncertainty intervals. Provenance records the full config and its hash, truth/U_ref, exact coordinates, coordinate-manifest hash, source/helper hashes, package versions, self-test results, numerical settings, definitions/units and figure paths. Machine-readable JSON uses null for unavailable quantities.

Self-tests check geometry extents, unique/outer counts, the single DOPPIO intersection, nested LATTE sets, exact inner/outer counts, 90° rotational and both reflection invariances, unchanged preceding-revision transects, exact HIGH-coordinate reuse, error definitions, final Omega and units, streamfunction signs, common U_ref and finite accepted clean fits. The two notebooks are also executed from clean kernels. Older result files are retained separately; they describe the previous designs and should not be combined with revised trials.

These tests assess fixed-geometry circular-Gaussian recovery under independent component noise or a uniform unmodelled current. They do not assess correlated errors, variable survey coverage, temporal evolution, background shear, non-Gaussian structure, 3D recovery or platform-specific observational errors. One tested direction for DOPPIO/LATTE does not establish direction independence. Matching truth does not equalise geometry, observation budgets or free shape parameters. Derive proposed paper statements from the supplied current images/tables, with these limits, rather than treating the setup description as a result or a universal method ranking.
