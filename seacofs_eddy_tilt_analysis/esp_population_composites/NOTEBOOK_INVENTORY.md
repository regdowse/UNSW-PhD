# Repository notebook inventory

Design-stage review of the local checkout, 16 September 2026. This is a source inventory, not a rerun or validation of saved scientific results. Checkpoints are excluded.

## aviso_eddy_dataset/notebooks/Untitled.ipynb

Exploratory/code-only notebook; inspect source for purpose.

## aviso_eddy_dataset/notebooks/analyse_processed_dataset.ipynb

# AVISO processed eddy dataset: quick analysis This notebook opens the final processed Parquet table and provides simple quality-control summaries and plots. It does not modify the dataset. ## Dataset overview ## One row per eddy ## Basic distributions ## Spatial coverage Plot one mean position per tracked eddy to keep the figure lightweight. ## Temporal coverage ## Inspect an individual eddy Change `EDDY_ID` to inspect its track and fitted properties. #### General plots

## aviso_eddy_dataset/notebooks/compare_original_esp_reconstruction.ipynb

# AVISO velocity versus ESP reconstruction Compare the native AVISO velocity field with the velocity reconstructed from the fitted ESP/DOPPIO parameters. The notebook automatically selects three long-lived anticyclonic eddies (AEs) and three long-lived cyclonic eddies (CEs), then shows early, middle, and late stages of each life. Only the displayed AVISO field may have a constant far-field background removed. DOPPIO reconstruction always represents the fitted eddy component. ## Controls ## Load

## aviso_eddy_dataset/notebooks/eddy_velocity_overlay.ipynb

# AVISO eddy velocity overlay Plot native AVISO absolute geostrophic velocity and overlay the processed ESP/DOPPIO eddy centres, IDs, polarity, and fitted maximum-tangential-velocity contours. ## Load the processed dataset ## Plotting function Pass an integer pipeline `Day` or a date such as `"1994-06-01"`. The source NetCDF is taken from `source_file`, with an annual-file fallback. ## Plot a day The default selects the day containing the most processed eddies. Replace `DAY` with another integer

## aviso_eddy_dataset/notebooks/generate_daily_overlay_frames.ipynb

# Generate daily AVISO eddy-overlay frames Generate one PNG per available AVISO day for later assembly into an animation. Each frame follows `eddy_velocity_overlay.ipynb`: native `ugos`/`vgos` speed and arrows, processed eddy centres and IDs, and the fitted $R_c^2/2$ ESP contour. The NetCDF files are processed one year at a time, figures are closed immediately, existing frames can be skipped, and a CSV manifest records the animation order. ## Controls ## Load eddies and source files ## Determine

## aviso_eddy_dataset/notebooks/run_pipeline.ipynb

# AVISO Eddy Dataset Pipeline This notebook is a concise control panel for the AVISO surface eddy pipeline. Scientific and orchestration code lives in `src/aviso_eddy_dataset/`; this notebook loads configuration, runs stages, and checks outputs. ## Load configuration The notebook uses `config/local.yaml` when present and otherwise falls back to the committed `config/example.yaml`. ## Run one stage Uncomment a stage while developing or resuming the workflow. ## Run selected stages Uncomment the s

## seacofs_eddy_dataset_modular/notebooks/Untitled.ipynb

Exploratory/code-only notebook; inspect source for purpose.

## seacofs_eddy_dataset_modular/notebooks/dataset_analysis.ipynb

Exploratory/code-only notebook; inspect source for purpose.

## seacofs_eddy_dataset_modular/notebooks/eddy_velocity_overlay.ipynb

# SEACOFS eddy velocity overlay Plot rotated SEACOFS surface velocity from the source `outer_avg_*.nc` files and overlay the processed ESP/DOPPIO eddy centres, IDs, polarity, and fitted maximum-tangential-velocity contours. ## Load configuration and processed eddies ## Source-file helpers Pass an integer pipeline `Day` or a date such as `"2018-04-01"`. The notebook first uses `fname` from the processed table, then falls back to scanning model-file `ocean_time` values. ## Plotting function ## Plo

## seacofs_eddy_dataset_modular/notebooks/esp_eddy_vel_reco.ipynb

# ESP 3-D velocity-field reconstruction This notebook demonstrates how well the elliptical streamfunction parameterisation (ESP) reconstructs one SEACOFS eddy from its depth-dependent fitted parameters. It resolves the original monthly NetCDF file from `fname`, applies the same depth interpolation and grid rotation used by the pipeline, evaluates ESP on the native grid, and compares maps, vertical sections, depth-resolved skill, 3-D velocity vectors, and normalised Okubo–Weiss isosurfaces. By de

## seacofs_eddy_dataset_modular/notebooks/generate_daily_overlay_frames.ipynb

# Generate daily SEACOFS eddy-overlay frames Generate one PNG per available model day for later assembly into an animation. Each frame follows `eddy_velocity_overlay.ipynb`: rotated surface speed and arrows, processed eddy centres and IDs, and the fitted $R_c^2/2$ ESP contour. The NetCDF files are processed one file at a time, figures are closed immediately, existing frames can be skipped, and a CSV manifest records the animation order. ## Controls ## Load eddies and source files ## Determine on

## seacofs_eddy_dataset_modular/notebooks/run_pipeline.ipynb

# SEACOFS Eddy Dataset Pipeline This notebook is a concise control panel for the modular SEACOFS eddy dataset pipeline. The scientific functions live in `src/seacofs_eddy_dataset/`; this notebook only loads config, runs stages, and checks outputs. ## Load Config Copy `config/example.yaml` to `config/local.yaml` and edit paths before running the pipeline. ## Run Individual Stages Run any single stage while developing or resuming part of the workflow. ## Run The Workflow Edit this list if you want

## seacofs_eddy_tilt_analysis/Untitled.ipynb

Exploratory/code-only notebook; inspect source for purpose.

## seacofs_eddy_tilt_analysis/beta_effect_background_flow/01_build_background_cache.ipynb

# Build the background-flow cache Run this notebook once before the background-flow analyses. It uses the current Jupyter kernel, avoiding differences between the notebook environment and the shell Python. The builder is restartable: completed monthly-file partitions are skipped if the notebook is interrupted and rerun. ## Load and select eddies The cache is built for all eddy observations. The topographic PV-gradient calculation retains `(w + f)`, but planetary/topographic dominance is delibera

## seacofs_eddy_tilt_analysis/beta_effect_background_flow/background_relative_beta_effect.ipynb

# Background-relative beta-effect analysis This notebook tests whether EAC anticyclonic eddies (AEs) have a northward/equatorward residual drift and cyclonic eddies (CEs) a southward/poleward residual drift after subtracting the strong southwestward background flow. It uses the cache produced by `build_background_cache.py`. The primary sample has planetary PV-gradient magnitude greater than topographic PV-gradient magnitude. The topographic term uses $(w+f)\nabla h/h^2$, so the dominance mask is

## seacofs_eddy_tilt_analysis/beta_effect_background_flow/planetary_eddy_vertical_flow_shear.ipynb

# Planetary-dominated eddies: vertical background flow and tilt This notebook asks whether depth-dependent environmental flow can explain the shared zonal component of planetary-dominated eddy tilt, while the beta effect explains the opposite meridional components. It uses exactly the same 2:1 planetary-dominance, 3000 m minimum depth, and 5 km minimum tilt criteria as `plan_topo_dom_eddies/planetary_dominated_eddies.ipynb`. `TiltDir` is a compass bearing **from the deep centre to the shallow ce

## seacofs_eddy_tilt_analysis/beta_effect_background_flow/vorticity_loss_during_beta_drift.ipynb

# Does beta drift weaken eddy relative vorticity? This notebook tests whether relative-vorticity magnitude decreases while AEs drift equatorward and CEs drift poleward. With the dataset convention `w > 0` for AEs and `w < 0` for CEs, both polarities should have $d|w|/dt<0$ during beta-aligned drift. The primary population is composed of contiguous off-shelf, planetary-PV-gradient-dominant trajectory segments. Track velocity and vorticity tendency are calculated from complete trajectories **befor

## seacofs_eddy_tilt_analysis/beta_effect_background_relative_propagation.ipynb

# Beta effect after removing background EAC advection This notebook asks whether EAC eddies have a polarity-dependent **residual drift** after subtracting the strong southwestward background current. It extends `beta_effect_vorticity_tilt_analysis.ipynb`. The primary diagnostic is $$\mathbf c_{res}=\mathbf c_{track}-\mathbf U_{background},$$ with the Southern Hemisphere beta-effect prediction $$c_{res,N}>0\quad\mathrm{for\ AEs},\qquad c_{res,N}<0\quad\mathrm{for\ CEs}.$$ Four background estimate

## seacofs_eddy_tilt_analysis/beta_effect_vorticity_tilt_analysis.ipynb

# Beta effect, relative vorticity, and eddy tilt This notebook tests whether the EAC eddy dataset contains a linked sequence of signatures expected from the planetary beta effect: 1. anticyclonic eddies (AEs) propagate equatorward and cyclonic eddies (CEs) propagate poleward; 2. the observed relative-vorticity tendency is compatible with $d\zeta/dt \approx -\beta v_y$; 3. within a given eddy, weaker $|\zeta|$ is associated with greater tilt magnitude; 4. vorticity weakening precedes tilt growth;

## seacofs_eddy_tilt_analysis/case_studies/01_planetary_reference_cases.ipynb

# Story 1: planetary-dominated reference eddies Select stable deep-water AE and CE examples that establish the background polarity-dependent tilt. Under the corrected signed-PV convention, AEs are expected to tilt along signed $\nabla PV$ and CEs opposite it. These are reference states, not examples of every offshore eddy. ## Load and classify Core-mean bathymetry and gradients are used. The main regime requires planetary dominance by at least 2:1 and depth of at least 3,000 m. Bearings with til

## seacofs_eddy_tilt_analysis/case_studies/02_topographic_response_cases.ipynb

# Story 2: polarity-dependent topographic response Select strongly topographic eddies that illustrate the asymmetric response found in the population analysis. CEs are ranked for coherent opposition to the signed total PV gradient; AEs are ranked for loss of a single coherent PV-relative direction. The AE category is therefore a disrupted-response example, not a claim that every topographic AE is random. ## Load and classify The primary topographic regime requires the smoothed topographic PV-gra

## seacofs_eddy_tilt_analysis/case_studies/03_regime_transition_cases.ipynb

# Story 3: within-eddy PV-regime transitions Find long-lived eddies that spend sustained periods in both deep planetary-dominated water and shallower topographic-dominated water. These are the strongest explanatory cases because each eddy provides its own reference. CE ranking rewards following the rotating opposite-PV target; AE ranking rewards a clean planetary reference followed by weaker topographic directional coherence. ## Load and classify A transition must contain sustained qualifying ob

## seacofs_eddy_tilt_analysis/case_studies/eddy_cross_sections_3d/eddy_esp_composite.ipynb

# Daily ESP fields and a mean three-dimensional eddy This notebook reconstructs the ESP velocity field independently for every selected day of one eddy, recentres each field on its shallowest fitted centre, and then composites the evaluated velocity fields. It does not average ESP parameters before reconstruction. The default onshore frame rotates each day using the local core-mean bathymetric gradient: positive axis 1 is onshore and positive axis 2 is alongshore. This prevents southward transla

## seacofs_eddy_tilt_analysis/case_studies/eddy_cross_sections_3d/eddy_sections_and_esp-Copy1.ipynb

# Eddy hydrographic sections and 3-D ESP reconstruction Select an eddy and optionally a day, then inspect native-level zonal and meridional sections of temperature, salinity, surface-referenced potential density, $N^2$, eastward velocity, northward velocity and speed. The final cells compare the background-removed SEACOFS velocity with the depth-dependent ESP reconstruction. The saved N² cache contains eddy-core summary diagnostics rather than a spatial field. It is loaded for case-level validat

## seacofs_eddy_tilt_analysis/case_studies/eddy_cross_sections_3d/eddy_sections_and_esp.ipynb

# Eddy hydrographic sections and 3-D ESP reconstruction Select an eddy and optionally a day, then inspect native-level zonal and meridional sections of temperature, salinity, surface-referenced potential density, $N^2$, eastward velocity, northward velocity and speed. The final cells compare the background-removed SEACOFS velocity with the depth-dependent ESP reconstruction. The saved N² cache contains eddy-core summary diagnostics rather than a spatial field. It is loaded for case-level validat

## seacofs_eddy_tilt_analysis/case_studies/long_eddies.ipynb

Exploratory/code-only notebook; inspect source for purpose.

## seacofs_eddy_tilt_analysis/case_studies/long_eddies_depth.ipynb

# Long-lived eddies: depth-following PV-gradient case studies This is the depth-resolved counterpart of `long_eddies.ipynb`. It loads the cached depth-following tables rather than recalculating ellipse means. Each colour represents the same fixed cached depth in every time-series panel and on the map. Values are selected directly at those model levels: no vertical interpolation or extrapolation is performed. Blue and orange background shading use the depth-averaged snapshot classification so eac

## seacofs_eddy_tilt_analysis/case_studies/long_eddies_depth_snapshot.ipynb

# Long-lived eddies: depth-averaged snapshot case studies This is the depth-averaged companion to `long_eddies_depth.ipynb`. It loads the cached `source='depth_snapshot'` table, which contains one thickness-weighted, column-averaged row per eddy-day. The panels therefore show one curve per variable rather than separate fixed-depth curves. Blue and orange background shading identify planetary- and topographic-beta-dominated snapshots. ## Settings and cached data ## Select the longest-lived eddies

## seacofs_eddy_tilt_analysis/case_studies/plan.ipynb

Exploratory/code-only notebook; inspect source for purpose.

## seacofs_eddy_tilt_analysis/case_studies/pv_tilt_transition_examples/00_build_selected_eddy_cache.ipynb

# Build the selected-eddy fixed-core cache Calculate the ESP-Gaussian surface PV-gradient terms for the supplied eight AEs and four CEs only. Run this first after changing the selected eddy list.

## seacofs_eddy_tilt_analysis/case_studies/pv_tilt_transition_examples/01_eddy_overviews.ipynb

# Selected eddy overviews For each supplied eddy, show compact tilt/PV time series beside its trajectory over bathymetry. Regime colours: green planetary, orange topographic and grey mixed. Blue arrows are surface-to-deep tilt; magenta arrows are the fixed-core mean environmental PV-gradient direction. PV-gradient arrows are normalised to show direction; their magnitude is shown in the time series. Arrow opacity increases with vector coherence, so cancellation-dominated directions are visually d

## seacofs_eddy_tilt_analysis/case_studies/pv_tilt_transition_examples/02_selected_day_core_maps.ipynb

# Selected-day core PV-gradient and eddy-spine maps Create individual maps for representative planetary, transition and topographic days. Edit `DAY_OVERRIDES` to replace the automatically suggested days for any eddy. White arrows show local environmental PV-gradient directions inside the `FRAC=1` ellipse. The magenta arrow is their Gaussian-weighted mean. Connected depth-coloured points are the fitted eddy centres down the water column on that day.

## seacofs_eddy_tilt_analysis/case_studies/topo.ipynb

Exploratory/code-only notebook; inspect source for purpose.

## seacofs_eddy_tilt_analysis/case_studies/trans.ipynb

Exploratory/code-only notebook; inspect source for purpose.

## seacofs_eddy_tilt_analysis/dep_or_indep_tilt_factors.ipynb

# Secondary Tilt Factors Screen tilt distance against eddy structure, stratification, and propagation diagnostics.

## seacofs_eddy_tilt_analysis/depth_following_pv/build_depth_following_pv_cache.ipynb

# Build depth-following PV-gradient tables This notebook follows each fitted eddy ellipse through the upper 1000 m. It calculates horizontal ellipse means at every fitted depth, then forms a thickness-weighted eddy-day table by averaging vector components before reconstructing magnitudes and bearings. The calculation is intentionally cached because it is substantially more expensive than the surface-centred method. ## Settings The output directory is separate from the source eddy tables. Existin

## seacofs_eddy_tilt_analysis/depth_resolved_pv_tilt/01_depth_data_coverage.ipynb

# 1. Depth-cache coverage and eddy-spine displacement Before comparing physics, this notebook visualises which Eddy–Day observations survive at each cached depth. Formal comparisons use only observations available at every selected level, preventing depth-dependent sample attrition from masquerading as a physical signal. ## Reading this notebook Use the coverage figures to decide whether the deepest selected level is sufficiently represented. The paired notebooks deliberately use the matched sam

## seacofs_eddy_tilt_analysis/depth_resolved_pv_tilt/02_surface_depth_pv_comparison.ipynb

# 2. Surface versus depth-following PV-gradient environment This notebook asks how far the inferred bathymetry and PV-gradient vector move away from their shallow-level values as the detected eddy spine shifts horizontally. Large rotations with increasing displacement support the geometric hypothesis: the deeper core samples a different bathymetric setting rather than merely a stronger version of the surface gradient.

## seacofs_eddy_tilt_analysis/depth_resolved_pv_tilt/03_regime_switching_with_depth.ipynb

# 3. Planetary–topographic regime switching with depth Regimes use a 2:1 threshold. The primary result is the fraction of matched snapshots whose depth-following classification differs from the shallow classification.

## seacofs_eddy_tilt_analysis/depth_resolved_pv_tilt/04_depth_resolved_tilt_direction.ipynb

# 4. Which depth's PV-gradient direction best matches tilt? `TiltDir` is the existing whole-column tilt direction. Agreement is polarity-aware: AEs are expected along signed grad(PV), CEs opposite it. Each eddy contributes equally to the summary. A depth-dependent minimum in directional error would identify a better explanatory sampling level. Interpret effect size and consistency across polarities, not only significance.

## seacofs_eddy_tilt_analysis/depth_resolved_pv_tilt/05_depth_resolved_tilt_magnitude.ipynb

# 5. Tilt magnitude and depth-dependent environmental mismatch This notebook tests whether large tilt occurs where deeper parts of the eddy encounter bathymetry or PV conditions missed by the shallow centre.

## seacofs_eddy_tilt_analysis/depth_resolved_pv_tilt/06_best_pv_representation.ipynb

# 6. Fixed-depth versus column-vector PV representation This notebook places the fixed-level estimates and the cached thickness-weighted 0–1000 m vector mean on one visual scorecard. Fixed levels use the same matched Eddy–Day sample. Select a preferred representation from consistent effect size, sample coverage and physical interpretation. The column-vector mean is a distinct whole-column hypothesis, not a substitute for examining where along the spine the signal changes.

## seacofs_eddy_tilt_analysis/depth_resolved_pv_tilt/07_topographic_mismatch_case_studies.ipynb

# 7. Case studies where the shallow centre misses deep topography Candidates begin planetary dominated at the shallowest level but become topographic dominated at depth. Ranking is only a screening device; every track must be visually inspected for boundary effects, coverage and physical coherence. ## Visual audit Prefer cases with sustained rather than single-day switching, continuous depth coverage, coherent tracks, and a visibly plausible bathymetric feature beneath the displaced deep centre.

## seacofs_eddy_tilt_analysis/eddy_categories/categorising.ipynb

#### off shelf #### Shear

## seacofs_eddy_tilt_analysis/eddy_census.ipynb

# Eddy Census Counts and distributions eddies in the EAC region.

## seacofs_eddy_tilt_analysis/ellipse_tilt_analysis/ellipse_orientation_tilt.ipynb

# Eddy deformation, ellipse orientation and tilt **Primary question:** does surface elongation relate to measured tilt magnitude, and does tilt have a preferred direction relative to the ellipse? Parallel alignment is a possibility, not an assumed outcome. AE and CE are analysed separately. **Secondary question:** do ellipse geometries at 50, 100, 200 and 300 m show stronger relationships? Every depth uses the **same existing whole-column `TiltDis` and `TiltDir`**; this notebook never recomputes

## seacofs_eddy_tilt_analysis/long_eddy_case_study.ipynb

# Long Eddy Case Study Select strongly propagating long-lived AE/CE examples and compare tilt with surface and bottom motion.

## seacofs_eddy_tilt_analysis/ml_subsurface_predictive_modelling.ipynb

# Predicting eddy tilt magnitude and direction from surface information This workflow answers two different questions for **AE and CE separately**: 1. Can surface structure and environmental variables predict tilt in an unseen eddy? 2. Which feature groups have stable descriptive or predictive associations with tilt? Tilt magnitude and circular direction are modelled separately. Model family and feature-set selection occur inside every outer eddy-grouped fold, so no repeatedly inspected final se

## seacofs_eddy_tilt_analysis/plan_topo_dom_eddies/ae_shelf_constraint_hypothesis.ipynb

# AE signed-PV alignment versus shelf constraint ## Updated hypothesis and present evidence The corrected hypothesis is that AEs align with signed $\nabla PV$, whereas CEs oppose it. The completed notebooks currently show: - Planetary-dominated AEs have a broadly equatorward tendency, but alignment is moderate rather than tight (median mismatch about $61.8^\circ$). - As topographic influence increases, AE preference error rises and the aligned fraction falls. This is **loss of PV-aligned directi

## seacofs_eddy_tilt_analysis/plan_topo_dom_eddies/plan_dom_sensitivity.ipynb

Exploratory/code-only notebook; inspect source for purpose.

## seacofs_eddy_tilt_analysis/plan_topo_dom_eddies/planetary_dominated_eddies.ipynb

# Planetary-dominated eddies: background tilt direction ## Corrected signed-PV hypothesis This uses the mathematical gradient of **signed** PV. In the Southern Hemisphere, $\nabla f/h$ points northward (equatorward). AEs are predicted to **align** with signed $\nabla PV$; CEs are predicted to **oppose** it. Both may retain a similar westward component, so this is an angular tendency rather than an exact bearing. Planetary-dominated observations define the background bearing used later. Cached ou

## seacofs_eddy_tilt_analysis/plan_topo_dom_eddies/rossby_number_pv_gradient.ipynb

# Rossby number and the topographic PV gradient ## Aim For shallow-water PV, `PV = (zeta + f) / h`. Neglecting horizontal gradients of relative vorticity gives `grad(PV) = grad(f)/h - (zeta + f) grad(h)/h**2`. With the **signed** Rossby number `Ro = zeta/f`, the topographic contribution is `beta_topo = -f (1 + Ro) grad(h)/h**2`. This notebook tests whether CE relative vorticity amplifies the same environmental slope while AE relative vorticity cancels it, and whether this helps explain the polar

## seacofs_eddy_tilt_analysis/plan_topo_dom_eddies/topographic_dominated_eddies.ipynb

# Topographic-dominated eddies: polarity-dependent response to signed PV gradients ## Corrected hypothesis For signed $PV=(\zeta+f)/h$, AEs are predicted to **align** with $\nabla PV$ ($\Delta\theta\to0^\circ$), whereas CEs are predicted to **oppose** it ($\Delta\theta\to180^\circ$). The CE transition is: natural poleward tilt → departure as total signed $\nabla PV$ rotates → increasing opposition to the signed topographic gradient. The secondary AE hypothesis is that shelf proximity disrupts it

## seacofs_eddy_tilt_analysis/pv_tilt.ipynb

# PV Tilt Compare eddy tilt direction with planetary, topographic, and total shallow-water PV-gradient directions. ## PV-gradient regional rose plot The same regional rose layout can be used with PV-gradient magnitude and direction. #### PV dominance

## seacofs_eddy_tilt_analysis/rose_plots.ipynb

# Windrose Tilt Summaries Direction and magnitude summaries for AE/CE tilt vectors. ## Regional rose plot This uses the shared `rose_plot` helper so the same regional windrose layout can be reused in other notebooks.

## seacofs_eddy_tilt_analysis/shear_lift_tilt/rotation_dependent_shear_lift.ipynb

# Does vertical shear produce rotation-dependent EAC eddy tilt? This notebook tests the mechanism suggested by Chao and Shaw (1998): differential advection produces an along-flow displacement, while interaction between the ambient flow and eddy rotation produces a transverse response whose sign reverses with rotation. The Arctic result is a hypothesis to test, not a sign convention to impose on the Southern Hemisphere EAC. The existing `TiltDir` and `TiltDis` are the authoritative coherent tilt

## seacofs_eddy_tilt_analysis/shear_tilt_climatology/00_prepare_inputs.ipynb

# Prepare the population shear-tilt inputs Run once on Katana, then run **01_shear_tilt_climatology.ipynb**. No individual case selection and no new ROMS velocity/temperature processing. Reuses the v4 background-flow cache and authoritative coherent `TiltDis`/`TiltDir`. The current **ESP-Gaussian, nonlinear, FRAC=1** environmental PV-gradient method supplies the regime diagnostic. This does not use `PV_grad_full` or presume that gradient dominance proves dynamical control. The output retains eve

## seacofs_eddy_tilt_analysis/shear_tilt_climatology/01_shear_tilt_climatology.ipynb

# Climatology of shear-relative eddy tilt and its evolution **Questions:** Do AEs and CEs have a preferred along/against/across-shear orientation? Does it differ between planetary and topographic conditions? Given the current orientation, how does coherent tilt change over the next 1, 3 and 5 days? This is a population analysis; no case studies or hand-selected eddies. Positive parallel = downshear; negative = upshear; positive transverse = left. The primary tilt is the existing measured lower-t

## seacofs_eddy_tilt_analysis/shear_tilt_climatology/02_shear_direction_preference.ipynb

# Does the tilt preference follow the direction of background shear? **Primary question:** do AEs tilt left and CEs right across different geographic shear directions, separately in planetary and topographic eddy-days? We also test whether topographic eddies encounter more variable shear. This is a hypothesis, not a selection assumption. No individual case studies are selected. The mixed population and **all quality-controlled eddy-days**, including unknown PV labels, test generalisation directl

## seacofs_eddy_tilt_analysis/statistical_tilt_direction_associations.ipynb

# Statistical associations with eddy tilt direction This notebook asks which surface and environmental feature groups have stable associations with **circular tilt direction**, independently for AE and CE. Direction is analysed only for tilt magnitude ≥5 km. For scalability, two clustered GEE models estimate the eastward and northward components of the unit tilt vector. Their predictions are recombined into direction and resultant concentration. This respects angular wrap-around and repeated edd

## seacofs_eddy_tilt_analysis/statistical_tilt_magnitude_associations.ipynb

# Statistical associations with eddy tilt magnitude This notebook addresses the explanatory question separately from ML prediction: **which surface and environmental feature groups have stable associations with tilt magnitude?** AE and CE are analysed independently. The primary model is a Gaussian generalized estimating equation (GEE) for `log1p(TiltDis)`. Eddy is the clustering unit, giving cluster-robust uncertainty for repeated eddy-days. Nonlinear spline terms describe effect shape. Time-var

## seacofs_eddy_tilt_analysis/surface_pv_esp_gaussian/00_build_esp_gaussian_cache.ipynb

# Build the fixed-core ESP-Gaussian cache Reconstruct $\zeta=w\exp(-\rho^2/R_c^2)$ and compare Gaussian weighting with an equal-weight control. Both methods use the same physically fixed `FRAC=1` radius-of-maximum-tangential-velocity core. Run this notebook first on Katana.

## seacofs_eddy_tilt_analysis/surface_pv_esp_gaussian/01_gaussian_reconstruction.ipynb

# What the fixed-core ESP-Gaussian reconstruction measures Visualise the fitted vorticity, Gaussian weights, and PV-gradient fields using only seabed cells within the `FRAC=1` dynamical core. The cyan ellipse is both the radius-of-maximum-velocity boundary and the hard sampling cutoff. The Gaussian remains 0.607 at this boundary; it prioritises the inner core without discarding dynamically relevant outer-core cells.

## seacofs_eddy_tilt_analysis/surface_pv_esp_gaussian/02_fixed_core_gaussian_weighting.ipynb

# Gaussian weighting inside the fixed dynamical core `FRAC=1` is not tuned. This notebook instead shows how the Gaussian redistributes influence among cells inside that fixed boundary and how much effective sample size is retained.

## seacofs_eddy_tilt_analysis/surface_pv_esp_gaussian/03_uniform_vs_gaussian.ipynb

# Uniform versus Gaussian weighting within FRAC=1 Both calculations sample exactly the same seabed cells inside the radius-of-maximum-velocity core. Only their spatial weights and reconstructed vorticity differ.

## seacofs_eddy_tilt_analysis/surface_pv_esp_gaussian/04_environmental_vs_full_pv.ipynb

# Environmental forcing versus total reconstructed PV Keep the environmental vector as the forcing estimate, then test how much adding the eddy's own $\nabla\zeta/h$ changes the result.

## seacofs_eddy_tilt_analysis/surface_pv_esp_gaussian/05_seamount_encounter_examples.ipynb

# Low-coherence bathymetric-edge encounters within FRAC=1 Select snapshots only after requiring strong local topographic exposure and weak vector coherence. Maps show the fixed core, local gradient directions, and the single Gaussian-weighted environmental vector. These are bathymetric-feature candidates; isolated seamounts should be distinguished from the continental slope during interpretation.

## seacofs_eddy_tilt_analysis/surface_pv_esp_gaussian/06_tilt_relationships_method_assessment.ipynb

# Fixed-core weighting assessment against eddy tilt Compare Gaussian and uniform weighting only within the same `FRAC=1` physical core. The tilt association is an evaluation, not a criterion for changing the boundary. ## Interpretation rule `esp_gaussian_1` is the primary method because it uses the complete ESP structure while remaining strictly inside the radius-of-maximum-velocity core. Report its environmental `PV_grad_mag` and `PV_grad_theta` as the single forcing estimate. Use local magnitu

## seacofs_eddy_tilt_analysis/surface_pv_footprint_sensitivity/00_build_surface_footprint_cache.ipynb

# 0. Build the surface-footprint sensitivity cache This is the only expensive notebook. It applies the corrected surface calculation to prespecified filled ellipses and annuli, plus two legacy controls. The depth-following cache is not read or changed.

## seacofs_eddy_tilt_analysis/surface_pv_footprint_sensitivity/01_legacy_vs_nonlinear.ipynb

# 1. Historical versus corrected nonlinear averaging This isolates the mathematical averaging change at the same surface centre and footprint size. The diagonal represents identical results.

## seacofs_eddy_tilt_analysis/surface_pv_footprint_sensitivity/02_fraction_sensitivity.ipynb

# 2. Filled-ellipse fraction sensitivity Shrinking the ellipse changes both feature detection and sample support. These plots retain the net vector and non-cancelling exposure separately.

## seacofs_eddy_tilt_analysis/surface_pv_footprint_sensitivity/03_annulus_sensitivity.ipynb

# 3. Does the outer eddy detect features missed by its centre? Annuli identify where within the surface ellipse strong topographic gradients occur. The 0.75–1.0 annulus represents the outer core, where a nearby seamount or slope may first intersect the eddy.

## seacofs_eddy_tilt_analysis/surface_pv_footprint_sensitivity/04_seamount_detection.ipynb

# 4. Detect strong but cancelling topographic features Candidates combine high local-gradient exposure with low vector coherence. This detects heterogeneous bathymetry without pretending that opposing gradients define a single net direction. Maps are a mandatory audit. Reject land-boundary artefacts, grid-edge cases and isolated single-cell gradients before interpreting a candidate as seamount interaction.

## seacofs_eddy_tilt_analysis/surface_pv_footprint_sensitivity/05_tilt_relationship_stability.ipynb

# 5. Stability of tilt relationships across surface footprints The best footprint should not be selected from the strongest isolated relationship. We look for polarity-consistent behavior that remains interpretable across neighbouring scales.

## seacofs_eddy_tilt_analysis/surface_pv_footprint_sensitivity/06_footprint_recommendation.ipynb

# 6. Evidence synthesis and footprint recommendation This notebook visualises the trade-off between feature detection, directional coherence, sample support and tilt agreement. It does not automatically declare the fraction with the strongest tilt relationship to be correct. ## Decision rules Prefer a footprint that has adequate grid-cell support, detects visually verified features, gives similar conclusions at adjacent fractions, and does not collapse net-vector coherence unnecessarily. Report

## seacofs_eddy_tilt_analysis/tilt_analysis.ipynb

# Vertically Checked Tilt Analysis Main summary plots for vertically checked tilt distance.

## seacofs_eddy_tilt_analysis/tilt_census.ipynb

# Tilt census A tilt-distance-specific census for anticyclonic (AE) and cyclonic (CE) eddies. The notebook retains the original mirrored-histogram idea but uses compact publication dashboards, smooth shading and embedded weighted box summaries. Daily distributions give every eddy equal total weight, preventing long tracks from dominating. Per-eddy metrics contain one observation per eddy. ## Census summary Counts remain available for documentation, while the figures below focus on tilt-distance

## seacofs_eddy_tilt_analysis/tilt_census_trial.ipynb

# Tilt-distance census: publication dashboard gallery Ten alternative **full-page, multi-panel dashboards** for presenting the tilt-distance census. Every comparison separates anticyclonic (AE) and cyclonic (CE) eddies. Direction is intentionally left for later analysis. Following the regional interpretation used here, the upstream sector contains **S1, U1 and U2**, while the downstream sector contains **S2, D1 and D2**. Daily distributions use equal total weight per eddy; per-eddy summaries con

## seacofs_eddy_tilt_analysis/tilt_mechanisms/00_build_stratification_cache.ipynb

# Build the eddy-centre stratification cache This notebook calculates environmental buoyancy frequency and writes one row per existing eddy-day. It **does not calculate or modify tilt**. `xroms.potential_density(..., z=0)` supplies surface-referenced potential density, then $N^2=-g\rho_0^{-1}\partial\sigma_0/\partial z$ is evaluated down each requested column. Using potential rather than pressure-dependent in-situ density prevents adiabatic compression from being misidentified as stratification.

## seacofs_eddy_tilt_analysis/tilt_mechanisms/01_vertical_shear.ipynb

# Does vertically sheared background flow explain measured tilt? Primary hypothesis: measured `TiltDir` aligns with the trailing integral of surface-minus-depth environmental velocity, and `TiltDis` increases with the magnitude of that accumulated displacement. This is tested separately for AE/CE, on/off shelf, background definition, depth, and accumulation window. Instantaneous alignment is descriptive; accumulated shear is the mechanistic test because tilt is a displacement. ## Decision rule S

## seacofs_eddy_tilt_analysis/tilt_mechanisms/02_beta_stratification_burger.ipynb

# Beta, stratification, and Burger-number hypotheses This notebook asks whether the equatorward increase in measured `TiltDis` is better described by beta alone, stratification, `N2/f²`, a transparent deformation-radius/Burger proxy, or background shear. Beta is inseparable from latitude, so spatial robustness and comparisons among competing predictors are central. ## Eddy-level magnitude relationships and robustness tests The following figures use one median observation per eddy. Raw relationsh

## seacofs_eddy_tilt_analysis/tilt_mechanisms/03_pv_topographic_steering.ipynb

# Planetary PV versus topographic steering This notebook tests whether tilt aligns with planetary PV gradients offshore but with topographic PV gradients over the shelf/slope. It reuses the established PV calculations in `seacofs_tilt_tools`; no PV formula is duplicated here. ## Strong topographic-steering evidence would be - tighter topographic-PV alignment as slope and bottom influence increase; - stronger alignment in topographic-dominant than planetary-dominant regimes; - a distinct on-shelf

## seacofs_eddy_tilt_analysis/tilt_mechanisms/04_depth_dependent_propagation.ipynb

# Depth-dependent centre propagation: a modal-dynamics proxy This notebook uses the existing vertical-profile dictionary to test the kinematic prediction \[d\mathbf{T}/dt \approx \mathbf{U}_{surface}-\mathbf{U}_{deep}.\] It is a practical precursor to a full baroclinic-mode decomposition. It does not re-estimate `TiltDis` or `TiltDir`. ## Escalation to full modal decomposition Proceed only if depth-dependent centre motion is coherent. Then solve vertical modes from the density profiles, estimate

## seacofs_eddy_tilt_analysis/tilt_mechanisms/05_wind_ekman_sensitivity.ipynb

# Optional wind/Ekman sensitivity Run this notebook only after creating an eddy-day wind-stress cache with geographic `tau_east_pa` and `tau_north_pa`. Wind is treated as a surface-forcing sensitivity, not assumed to be present in the existing BRAN velocity cache. The prediction is directional: deep-to-surface `TiltDir` should align with recent Ekman transport if wind displaces the shallow structure. ## Required extension The decisive wind test should use trailing 2, 5, 10 and 20-day stress inte

## seacofs_eddy_tilt_analysis/tilt_mechanisms/06_joint_mechanism_comparison.ipynb

# Joint comparison of candidate tilt mechanisms This final notebook is run after the component notebooks pass QC. It asks which mechanisms retain within-eddy and between-eddy associations after adjustment. Direction and magnitude are analysed separately, AE and CE are kept separate, and whole eddies remain the clustering unit. The goal is not to declare causality from a single coefficient. A mechanism is considered supported only when its directional, magnitude, timescale, regime and within-eddy

## seacofs_eddy_tilt_analysis/tilt_mechanisms/07_temporal_tilt_evolution.ipynb

# Temporal evolution of measured eddy tilt This notebook tests the literature prediction that vertical tilt accumulates as the upper and lower eddy centres propagate differently. It uses the existing `TiltDis` and `TiltDir` measurements and never re-estimates tilt. The analysis separates three related questions: 1. Does tilt magnitude grow within individual eddies as they age? 2. Is growth different for eddies persistently exposed to planetary- or topographic-PV dominance? 3. Do observed day-to-

## seacofs_eddy_tilt_analysis/tilt_mechanisms/08_stratification_polarity_tilt.ipynb

# Does stratification differ between AEs and CEs, and does it explain tilt magnitude? This notebook tests two related hypotheses without assuming either result: 1. **CEs are more strongly stratified than AEs.** 2. **More strongly stratified eddies have larger coherent tilt magnitudes (`TiltDis`).** The primary stratification measures are the eddy-core mean $N^2$ over 0–200 m and 0–500 m. The shallower measure should be more sensitive to cyclone/anticyclone thermocline displacement; the 0–500 m m
