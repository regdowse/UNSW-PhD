# Context for interpreting the ESP sensitivity figures and sampling overviews

This document accompanies four sensitivity-result images and the sampling overview placed first in each notebook, generated for potential inclusion in a revised Elliptical Streamfunction Parameterisation (ESP) methods paper. It is self-contained context for a ChatGPT project or another reader receiving the images without access to the analysis code. It describes the experimental design and figure encodings, **not the observed results**. No trends, recovered parameter values, observed success rates, or method rankings are asserted here.

The figures address two questions:

1. **Noise and sampling:** how does parameter recovery respond to controlled velocity noise and observation count, with spatial coverage held fixed?
2. **Uniform background flow:** how does adding an unmodelled constant velocity vector affect parameter recovery, including dependence on current direction?

Each question is examined with two truth assignments: the original method-specific shapes, and an additional circular eddy shared by all three methods. These are controlled synthetic sensitivity tests, not observational error estimates or a comprehensive uncertainty assessment.

## Identify the images

The numbering below is a reference for this document, not a proposed manuscript figure numbering. PNG and PDF versions have the same content. Filenames are relative to this directory.

| Figure | Image | Notebook | Synthetic truth assignment |
| --- | --- | --- | --- |
| 1 | `results/rc3_noise_sampling.png` | Second figure in `rc3_noise_sampling.ipynb` | SOLO: circular; DOPPIO and LATTE: elliptical |
| 2 | `results/rc3_noise_sampling_circular.png` | Third figure in `rc3_noise_sampling.ipynb` | All three methods: the same circular eddy |
| 3 | `results/rc3_background_flow.png` | Second figure in `rc3_background_flow.ipynb` | SOLO: circular; DOPPIO and LATTE: elliptical |
| 4 | `results/rc3_background_flow_circular.png` | Third figure in `rc3_background_flow.ipynb` | All three methods: the same circular eddy |

Figures 2 and 4 are additional comparisons. They do not replace Figures 1 and 3. SOLO's truth and observations are identical between Figures 1 and 2, and between Figures 3 and 4. For DOPPIO and LATTE, the added comparisons change the truth shape and its associated reference speed; the sampling and inversion settings are retained.

### Sampling overview: first figure in each notebook

Files: `results/rc3_noise_sampling_layout.png` and `results/rc3_background_flow_layout.png`, with matching PDF versions. These two overviews show the same geometry and truth assignments. They are setup diagrams, not additional recovery experiments.

Each has two rows and three method columns, left to right SOLO, DOPPIO and LATTE. Panels (a–c), the top row, show the original assignments: circular SOLO and elliptical DOPPIO/LATTE. Panels (d–f), the bottom row, show all methods on circular truth. All axes use kilometres, equal spatial aspect and the same −100 to 100 km limits.

Faint, small grey arrows show the unperturbed true eddy velocity, without noise or a background current. A thin grey outline marks the elliptical-radius contour \(\rho=R_c/\sqrt{2}\); the small cross marks the true centre. The outline is a reference contour of the Gaussian eddy, not a hard eddy boundary or sampling-selection boundary. There is no velocity colour scale.

Coloured hollow circles show the actual observation locations for every noise-experiment count, \(N=11,21,51,101\), using the same colour mapping as the noise-result figures. Marker size decreases with increasing count so coincident observations remain visible as concentric rings; locations are not jittered or moved to separate colours. Each colour represents a complete alternative sampling set, not an additional platform or a disjoint group of observations. In particular, nested LATTE sets and shared transect locations can have multiple rings at one position.

The count definitions below still apply: \(N\) is per transect for DOPPIO and total for SOLO/LATTE. The background notebook's overview deliberately includes all these counts for context, but its background-flow fits use **only \(N=101\)**. The overview does not imply that background-flow sensitivity was tested at every displayed count.

## Shared physical model and parameter definitions

The known, stationary, two-dimensional eddy is an elliptical Gaussian streamfunction:

\[
\psi(\mathbf{x})=\psi_0\exp\left[-\frac{\mathbf{r}^{T}Q\mathbf{r}}{R_c^2}\right],
\qquad \mathbf{r}=\mathbf{x}-\mathbf{x}_0.
\]

Here \(x\) increases eastward and \(y\) northward. The velocity convention is

\[
u=\frac{\partial\psi}{\partial y},\qquad
v=-\frac{\partial\psi}{\partial x},\qquad
\Omega=-\frac{2\psi_0}{R_c^2}.
\]

The true centre is \(\mathbf{x}_0=(0,0)\), the Gaussian radius parameter is \(R_c=85\) km, and the streamfunction amplitude is \(\psi_0=250{,}000\) m²/s for both shapes. The original idealised notebooks express the same amplitude as 250 km·m/s. Public calculations use SI units, with conversion to the original notebooks' km-coordinate fitting conventions at the inversion boundary.

\(Q\) is dimensionless, symmetric positive definite and has determinant one. The two truth shapes are:

| Truth shape | \(Q\) | True aspect ratio \(\alpha\) | Major-axis orientation |
| --- | --- | --- | --- |
| Circular | \(I\) | 1 | Undefined |
| Elliptical | \(\begin{bmatrix}1&0.4\\0.4&1\end{bmatrix}/\sqrt{0.84}\) | \(\sqrt{7/3}\approx1.528\) | 135° anticlockwise from east, modulo 180° |

The aspect ratio is major semiaxis divided by minor semiaxis: \(\alpha=\sqrt{\lambda_{\max}(Q)/\lambda_{\min}(Q)}\). The elliptical radius is \(\rho=\sqrt{\mathbf{r}^{T}Q\mathbf{r}}\). **\(R_c\) is the Gaussian scale**, not the radius of maximum normalised tangential speed; that peak is at \(\rho=R_c/\sqrt{2}\). Centre vorticity is \(\Omega\,\mathrm{tr}(Q)\), although vorticity is not a plotted metric. Orientation error is not shown in any of the four figures.

### Reference speed used on the horizontal axes

\(U_{\mathrm{ref}}\) is the maximum **vector-speed magnitude of the true full eddy field**, evaluated using the original velocity model at its analytic speed maximum. It is not estimated from each sample set or fitted field, and it does not include the imposed noise or background current.

The circular truth uses \(U_{\mathrm{ref}}\approx2.523\) m/s; the elliptical truth uses \(U_{\mathrm{ref}}\approx3.118\) m/s. These are setup normalisations, not recovered parameters. Each reference stays fixed across sampling counts and perturbation levels and is the same in the noise and background experiments for that truth shape.

Consequently, the same horizontal-axis percentage corresponds to different physical perturbation amplitudes for circular versus elliptical truth in Figures 1 and 3. In Figures 2 and 4, all methods share the same reference speed and physical perturbation amplitude at a given percentage.

## Methods and sampling geometries

The experiments reuse the original manuscript implementations of SOLO, DOPPIO and LATTE. Each estimates inner-core parameters before an outer Gaussian-profile fit estimates radius and amplitude. The cubic/affine inner-core approximations are not a joint exact fit of every Gaussian parameter. Zero perturbation therefore provides an empirical baseline to inspect; exact recovery should not be presumed from the synthetic truth alone.

SOLO prescribes \(Q=I\). DOPPIO and LATTE estimate \(Q\), including when their input eddy is circular. Their shape errors for circular truth measure departure from a true aspect ratio of one: **spurious ellipticity**. SOLO's shape is prescribed, so its missing shape-error curves are an intentional omission, not a zero error estimate or a failed fit.

| Method | Fixed observation geometry | Meaning of \(N\) |
| --- | --- | --- |
| SOLO | Straight eastward transect, \(x=-70\) to 70 km at \(y=20\) km | Total number of velocity-vector observations |
| DOPPIO | Orthogonal transects intersecting at (10,10) km. Horizontal: \(x=-20\) to 40 km at \(y=10\) km. Vertical: \(y=-20\) to 40 km at \(x=10\) km | Points **per transect**; \(2N-1\) unique observations because the intersection is shared |
| LATTE | Scattered points in the fixed square \([-50,50]\times[-50,50]\) km. Four fixed corner anchors; the remaining points are reproducible uniform random locations | Total number of velocity-vector observations |

For SOLO and DOPPIO, varying \(N\) changes spacing while retaining transect endpoints. DOPPIO uses odd counts to retain the exact shared intersection. LATTE uses nested point sets with fixed corners, so the domain and convex hull stay fixed as \(N\) changes. Observation locations do not change between replicates or perturbation levels at a given method/count. The background experiment uses exactly the corresponding noise-experiment geometry at its selected count.

All sampled observations are supplied to both fitting stages. No extra observations or truth-centred selection window are supplied. Initialisation rules, tuning and stopping rules remain fixed across conditions. In particular, the original SOLO call retains `core_thresh=30000` in its km-coordinate implementation; this covers its full 140 km transect rather than imposing a new 30 km fitting window. Original post-fit radius caps are retained: \(10^6\) km for SOLO/DOPPIO and \(10^5\) km for LATTE. These are numerical rejection limits, not physical truth values or optimisation bounds.

The geometries and numbers of unique observations differ between methods. Even the common circular-truth figures are **not a comparison at matched observation budget and geometry**. The original figures additionally use different truth shapes for SOLO versus the other methods.

## Error definitions and numerical acceptance

A hat denotes an estimate; a subscript \(t\) denotes known truth. The plotted total errors are:

\[
E_{\mathrm{centre}}=\|\widehat{\mathbf{x}}_0-\mathbf{x}_{0,t}\|
\quad\text{in km},
\]

\[
E_{R_c}=100\frac{|\widehat R_c-R_{c,t}|}{R_{c,t}},\qquad
E_\alpha=100\frac{|\widehat\alpha-\alpha_t|}{\alpha_t}.
\]

These are errors against truth, not residuals against the fitted observations. Total radius and aspect-ratio errors are absolute percentages and do not indicate the sign of overestimation or underestimation. Centre error is a displacement magnitude, not its direction. In-sample goodness of fit is not a plotted validation metric.

A full fit is accepted only when its optimisers complete successfully, parameters are finite, radius is positive, and \(Q\) passes symmetry, positive-definiteness and unit-determinant checks, without a radius-cap fallback. Acceptance is a numerical/physical validity check, **not an accuracy criterion or evidence that a parameter is well constrained**. No threshold on recovery error removes a fit. Failed fits retain explicit statuses in the trial tables; their full-fit metrics are excluded from plotted error summaries. Large finite errors remain in the figures when the fit passes acceptance checks.

## Figure 1 — noise and sampling, original truth assignments

**File:** `results/rc3_noise_sampling.png`.

SOLO receives the circular truth. DOPPIO and LATTE receive the elliptical truth. The figure asks how recovery and fit acceptance vary with noise amplitude and sampling density under those method-specific cases.

The observations are

\[
\mathbf{u}_{\mathrm{obs}}(\mathbf{x}_i)
=\mathbf{u}_{\mathrm{ESP}}(\mathbf{x}_i)+\boldsymbol\epsilon_i,
\qquad \boldsymbol\epsilon_i\sim\mathcal N(\mathbf{0},\sigma^2I_2).
\]

Errors are independent between unique locations and between the two velocity components. A shared DOPPIO intersection receives one vector-noise draw, reused on both transects. \(\sigma\) is the standard deviation of **each component**, not the standard deviation of vector-speed magnitude. There is no imposed background flow.

The horizontal axis is \(100\sigma/U_{\mathrm{ref}}\), with levels **0, 2.5, 5 and 10%**. Colours identify **\(N=11,21,51,101\)**, using the method-specific count definitions above. Each nonzero-noise condition has 100 independent noise replicates; each clean condition is fitted once. Randomness uses seed 20261027 and `numpy.random.default_rng`. Standard-normal draws are paired across noise amplitudes within a method/count/replicate; different method/count/replicate combinations use separate reproducible streams. Sampling-location variability is not replicated.

The layout is three method rows by four metric columns:

| Row | Centre error | Radius error | Aspect-ratio error | Accepted fits |
| --- | --- | --- | --- | --- |
| SOLO, circular | (a) | (b) | (c): intentionally not scored | (d) |
| DOPPIO, elliptical | (e) | (f) | (g) | (h) |
| LATTE, elliptical | (i) | (j) | (k) | (l) |

Error curves show the **median of individual accepted-fit errors**. Shading spans the **10th–90th percentiles** of those errors. It is simulation spread conditional on full-fit acceptance, not a confidence interval, standard error, or observational uncertainty band. At zero noise there is one fit, so no replicate spread is estimated. The accepted-fit panels show \(100\times\text{accepted}/\text{attempted}\) using all attempted fits in each condition; read these alongside the conditional error summaries. Their vertical axes can be zoomed near 100%, so inspect tick labels.

Corresponding metric columns share vertical scales across methods within this figure. Radius-error axes use a symmetric-log scale: linear from 0 to 1% error and logarithmic above 1%. Other metric axes are linear. Points mark tested conditions; connecting lines and shading between points do not represent additional simulated conditions.

## Figure 2 — noise and sampling, common circular truth

**File:** `results/rc3_noise_sampling_circular.png`.

This repeats Figure 1's design, panel layout, count definitions, colours, statistics and axis conventions, but **all three methods receive the same circular eddy**. Rows (a–d), (e–h) and (i–l) are therefore SOLO/circular, DOPPIO/circular and LATTE/circular, respectively. All methods use the same \(U_{\mathrm{ref}}\).

DOPPIO and LATTE still estimate unrestricted elliptical shapes. Their panels (g) and (k) measure spurious ellipticity relative to \(\alpha_t=1\); they are not orientation-error panels. SOLO's panel (c) remains intentionally unscored.

Standard-normal noise draws are paired with Figure 1. DOPPIO/LATTE physical noise amplitudes are rescaled by the circular reference speed; SOLO observations and fits are identical to Figure 1. This is a paired truth-shape comparison at fixed dimensionless perturbation levels, not a new independent SOLO ensemble. Vertical limits are chosen within each figure, so do not assume identical limits between Figures 1 and 2.

## Figure 3 — uniform background flow, original truth assignments

**File:** `results/rc3_background_flow.png`.

SOLO receives circular truth; DOPPIO and LATTE receive elliptical truth. Each method uses **\(N=101\)**, meaning 101 total observations for SOLO/LATTE and 201 unique observations for DOPPIO. This is a deterministic experiment with fixed locations, no added noise and no Monte Carlo replication.

The observations are

\[
\mathbf{u}_{\mathrm{obs}}(\mathbf{x}_i)
=\mathbf{u}_{\mathrm{ESP}}(\mathbf{x}_i)+\mathbf U_{\mathrm{bg}}.
\]

The same \(\mathbf U_{\mathrm{bg}}\) is added at every location. The inversion does **not** subtract that vector or fit an explicit background-current term. The horizontal axis is \(100|\mathbf U_{\mathrm{bg}}|/U_{\mathrm{ref}}\), at **0, 5, 10, 20, 30 and 50%**. Direction is either eastward \((+U,0)\) or northward \((0,+U)\) in the model coordinates. Eastward is parallel to the SOLO transect; northward is perpendicular. There is no imposed background shear.

Colours distinguish methods: SOLO blue, DOPPIO orange and LATTE green. Solid lines indicate eastward flow; dashed lines indicate northward flow. Legends also identify each method's truth shape. Both direction sequences include the same zero-current control for their method; these overlapping controls are not independent replicates.

The six panels distinguish errors relative to truth from changes relative to the method's own zero-current fit:

| Panel | Quantity |
| --- | --- |
| (a) | Total centre error relative to true centre, km |
| (b) | Total absolute radius error relative to true radius, % |
| (c) | Total absolute aspect-ratio error relative to true aspect ratio, %; DOPPIO/LATTE only |
| (d) | Magnitude of fitted-centre change from the zero-current estimate, km |
| (e) | Signed fitted-radius change from the zero-current estimate, % of that estimate |
| (f) | Signed fitted-aspect-ratio change from the zero-current estimate, % of that estimate; DOPPIO/LATTE only |

Writing \(\widehat p(B)\) for a fitted parameter under a background vector \(B\), the lower row uses

\[
\Delta_{\mathrm{centre}}(B)
=\|\widehat{\mathbf{x}}_0(B)-\widehat{\mathbf{x}}_0(0)\|
\quad\text{in km},
\]

\[
\Delta_{R_c}(B)=100\left[\frac{\widehat R_c(B)}{\widehat R_c(0)}-1\right],
\qquad
\Delta_\alpha(B)=100\left[\frac{\widehat\alpha(B)}{\widehat\alpha(0)}-1\right].
\]

Panel (d) is **not** the difference between the centre-error magnitudes in panel (a). Panels (e–f) use fitted baseline denominators, whereas panels (b–c) use true denominators. Negative lower-row radius/shape changes mean a decrease from the fitted baseline, not a negative total error. Zero-current total errors need not be zero; incremental changes are zero by definition when that baseline fit is accepted.

SOLO has no curves in (c) or (f), because it does not estimate shape. There are no uncertainty bands or error bars: each point is one deterministic fit. Lines connect tested magnitudes. The figure footer reports numerical acceptance; if a fit fails, its error value is omitted and an `x` is placed near the top of the panel at that horizontal position. Such a marker's vertical position is an annotation, not an error measurement. All vertical axes in Figure 3 are linear.

## Figure 4 — uniform background flow, common circular truth

**File:** `results/rc3_background_flow_circular.png`.

This repeats Figure 3 with **the same circular eddy and reference speed for all three methods**. The observation geometries, counts, current percentages and directions, inversion settings, colours, line styles, six panel definitions, and total-versus-incremental distinctions are retained. SOLO's observations and fits are identical to Figure 3.

DOPPIO/LATTE panel (c) measures spurious ellipticity relative to a circular truth. Panel (f) measures change from each method's own fitted zero-current aspect ratio, which need not equal the true value of one. Orientation remains unscored. SOLO shape curves remain absent.

In this figure, radius panels **(b) and (e) use symmetric-log vertical scales**: linear within ±1% and logarithmic beyond. This retains zero, signed incremental changes and large finite errors on the same axes. Other panels remain linear. A different radius-axis scale from Figure 3 is a display choice, not a different radius definition, fit constraint or acceptance rule. Read axis ticks rather than comparing apparent curve slopes or heights between the two background figures.

## Scope for subsequent interpretation

Use this document to establish what the images encode, then derive any results or proposed manuscript statements from the supplied images or associated tables. Do not treat the figure definitions as findings. In particular:

- Distinguish a baseline recovery error from the incremental response to an imposed perturbation.
- Distinguish numerical acceptance from accurate or identifiable parameter recovery; read noise-error summaries with their acceptance panels.
- The noise bands describe independent-component perturbations at fixed locations. Spatially correlated errors, randomised survey coverage and real measurement-error distributions are not tested.
- The background experiment tests spatially uniform currents, not shear, temporal evolution, vertical structure or background-flow removal.
- Shared circular truth improves comparability but does not equalise observation budgets, geometry or the number of free shape parameters. Avoid an unconditional ranking of algorithms from these plots.
- Do not infer orientation recovery from circular truth, or describe SOLO's prescribed aspect ratio as an estimated quantity.
- Do not report precise values that cannot be read from the images. Each figure has corresponding `*_trials.csv`, `*_summary.csv` and `*_provenance.json` files in `results/` if numerical support is needed. Provenance records the settings, truth parameters and source hashes; these files need to be supplied separately if the receiving assistant cannot access this directory.

The source methods are selected without modification from `solo_idealised_tests_updated.ipynb`, `doppio_idealised_tests_updated.ipynb` and `latte_idealised_tests.ipynb` in the manuscript's original `final_ESP_paper_notebooks` directory. Synthetic field evaluation and geometry conventions reuse that analysis's `modules/utils.py`. The local adapters and plotting definitions are in `rc3_validation.py` and `rc3_plotting.py`. This document describes the four result figures and their sampling overviews; it does not propose a selection of figures for the paper or supply a results interpretation.
