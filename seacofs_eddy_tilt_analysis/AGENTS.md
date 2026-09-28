# SEACOFS eddy tilt: agent context and navigation

Last reviewed: **2026-09-28**, against repository commit **09d5dd9** (before this
documentation addition). This is a maintained context snapshot, not a scientific
results archive. Paths below are relative to this directory unless stated otherwise.

## Start here

1. Work in **regdowse/UNSW-PhD**, primarily `seacofs_eddy_tilt_analysis/`.
   The enclosing local project is named **Eddy Dataset** and contains several
   old `UNSW-MRes*` checkouts. Those are historical, not the current destination.
2. Read this file, then the relevant folder README and notebook introduction,
   settings, helper functions and saved outputs. Do not reread every notebook
   or every chat for an ordinary focused task.
3. Check `git status --short --branch`, the remote, and recent commits before
   editing. Current collaboration uses `main`; the user has previously preferred
   direct pushes they can pull on Katana. Follow the current request, preserve
   unrelated edits and notebook outputs, and never force-push over others' work.
4. Identify the exact outcome, population, coordinate frame, weighting and
   cache settings. Similar notebook names do not imply identical definitions.
5. Use current implementation and input provenance to resolve discrepancies.
   This guide, READMEs, chat interpretations and saved outputs can become stale.
   Record a material discrepancy rather than silently changing the science.

## Scientific purpose

The project studies the three-dimensional structure and tilt of East Australian
Current (EAC) eddies in the 26-year SEACOFS model dataset. It supports a tilt
paper: describe the census, distinguish AE/CE behaviour, and test relationships
with planetary/topographic PV gradients, background-flow shear, stratification,
eddy geometry, lifecycle and vertical velocity. Case studies and ESP composites
help interpret population patterns. Associations and mechanistic consistency
are not automatically causal evidence.

The upstream sibling `../seacofs_eddy_dataset_modular/` performs:

`detect_nencioli → fit_doppio_surface → track_eddies → process_tracked_dataset
→ compute_vertical_profiles → qc_vertical_profiles → compute_tilt → analyse_tilt`.

Its `src/seacofs_eddy_dataset/core/` holds science kernels; `stages/` orchestrates
I/O and processing. Read its README, `config/local.yaml` and relevant kernel
when a request changes dataset construction. `../aviso_eddy_dataset/` is a
separate observational workflow; do not interchange its inputs with SEACOFS.

## Data contracts and conventions that must survive edits

| Quantity | Meaning / caution |
|---|---|
| `Eddy`, `Day` | Snapshot join keys. Never join by `Day` alone. Check uniqueness; vertical tables additionally require `Depth`. Do not assume IDs persist across dataset regeneration. |
| `Cyc` | `AE` = anticyclone, `CE` = cyclone. In this Southern Hemisphere dataset fitted relative vorticity is positive for AE and negative for CE. Plotting usually uses AE red, CE blue. |
| `TiltDis`, `TiltDir` | Authoritative saved coherent whole-column delta tilt: km and compass degrees, **deep-to-shallow** direction. Analysis normally loads these rather than refitting them. |
| Bearings | 0° north, 90° east; east component = distance × sin(bearing), north = distance × cos(bearing). Use circular differences, not ordinary subtraction/means across 0°. |
| Model geometry | `xc`, `yc`, `Rc` use km; depths/bathymetry use m. Model-grid axes, geographic axes, and rotated analysis axes differ. The production tilt bearing offset is currently 20°; inspect conversion helpers before applying any rotation. |
| Fitted `w` | Relative vorticity, s⁻¹. It is **not** ROMS vertical water velocity `w`, which is upward m/s. Helpers may call the latter `vv_*`. `Omega` is a separate fitted rotation parameter. |
| `h` | In core-mean PV tables, footprint-averaged bathymetry. This differs from centre-column `water_depth_m`, used in some shear coverage checks. |
| PV | A shallow-water environmental/topographic proxy, not full 3-D Ertel PV. Standard `PV_grad_*` excludes internal `grad(zeta)/h`; `PV_grad_full_*` includes it. |
| Rossby number | Signed `Ro = zeta/f` matters for `1 + Ro` in the topographic term; some descriptive plots explicitly use `abs(Ro)`. Do not exchange them. |
| Depth | Many analyses select nearest **actual fitted levels**, report them and avoid interpolation. Other helpers explicitly interpolate. A requested 500 m need not be an actual 500 m level. |
| Region | Shared labels are `S1,S2,U1,U2,D1,D2`; regional composites pool S1+S2, U1+U2, D1+D2 before averaging. Use `add_region_labels` rather than reinventing boundaries. |

**Different tilt estimands:** raw daily centre offsets, maximum projected extent,
saved delta `TiltDis`, and the displacement of a composite centreline are distinct.
Several newer composites show **surface-to-depth** displacement, opposite to the
saved deep-to-shallow direction. A mean vector can cancel across members and its
magnitude is not the mean of individual magnitudes. Keep labels explicit.

**Statistical support:** eddy-days repeat within tracks. Preserve the workflow's
day weighting versus equal-eddy weighting, use whole-eddy resampling where
implemented, and report both days and eddies. IQR shading is population spread;
bootstrap intervals quantify uncertainty in an estimator. A bootstrap seed does
not subsample the mean population unless a separate sample selection does so.
Matched-depth comparisons must match `(Eddy, Day)` after relevant QC. Core N²
is not background stratification. Radius-normalised tilt shares a denominator
with radius-related predictors. Depth ranking on the same data is descriptive.

## Shared APIs and actual input locations

`seacofs_tilt_tools.py` is the main shared module:

```python
paths = tilt.Paths()                 # override paths for a different system
grid = tilt.load_grid(paths.grid, paths.z_r)
surface, tilt_table = tilt.load_tilt_tables(paths, add_regions=True, grid=grid)
vertical = tilt.load_vert(paths)     # DataFrame
profiles = tilt.load_vert(paths, dic_form=True)  # EddyN -> DayN -> profile
```

Current primary root: `/srv/scratch/z5297792/SEACOFS_26yr_eddy_dataset_modular/`.

| Input | Relative to that scratch root unless absolute |
|---|---|
| Processed eddy-days | `processed/eddy_dataset_processed.parquet` |
| Saved tilt | `tilt/tilt_dataset.parquet` |
| Confirmed profiles | `vertical_profiles_confirmed/profiles.parquet` |
| Physical grid depths | `z_r.npy` |
| Shared surface PV cache | `pv_gradient_surface/surface_pv.parquet` |
| Depth-following PV cache | `pv_gradient_depth_following/pv_gradient_depth_following_{snapshot,depth}_0_1000m.parquet` |
| Model archive/reference grid | `/srv/scratch/z3533156/26year_BRAN2020/outer_avg_*.nc`; reference `outer_avg_01461.nc` |
| Corrected N² cache | `/srv/scratch/z5297792/SEACOFS_26yr_eddy_dataset/tilt_mechanisms/n2_eddy_day_v4_potential_density_core.parquet` — note the **non-modular** parent |
| External ESP backend | `/home/z5297792/ESP_zonodo/functions.py` |

Production data and expensive runs are on **Katana**. A local checkout or saved
notebook output does not establish that a current remote cache exists or matches
the code. Do not rewrite HPC paths to local paths in committed defaults. Core
dependencies include NumPy, pandas, SciPy, Matplotlib, netCDF4 and a Parquet
engine; specialised workflows add statsmodels, xroms or the ESP backend. Use the
relevant environment/configuration rather than assuming system Python suffices.

## PV calculation: current choices and version traps

Read `add_pv_gradient_terms` and the relevant notebook call before using it.
Its **API defaults** remain `core_mean=False`, `surface_method='uniform'`,
`use_max_abs_w=False`; these are not the explicit settings of many newer notebooks.

For the corrected ESP-Gaussian surface method:

```text
r = [x-xc, y-yc];  Q = [[q11,q12],[q12,q22]]
rho² = rᵀ Q r
core FRAC=1: rho² <= Rc²/2
G = exp(-rho²/Rc²)
zeta(r) = w G [1 - 2(rᵀ Q² r)/(Rc² trace(Q))]
g_plan = grad(f)/h
g_topo = -(f + zeta) grad(h)/h²
g_environment = g_plan + g_topo
g_full = g_environment + grad(zeta)/h
```

The Gaussian describes the streamfunction envelope; vorticity is its Laplacian.
Do not restore the superseded `zeta = w*G` reconstruction. Evaluate nonlinear
terms at each wet grid cell first, then average components using positive `G`
inside the core. Magnitude/bearing come from the averaged vector. Mean local
magnitude, RMS and coherence describe heterogeneity/cancellation and are not
interchangeable with the net magnitude. The envelope weight at the core boundary
is exp(-1/2), but vorticity there also contains the Laplacian shape factor.

Many current surface analyses explicitly use `core_mean=True`, `frac=1`,
`averaging='nonlinear'`, `surface_method='esp_gaussian'`, and
`use_max_abs_w=True`. That last option takes the **signed** largest-|w| value
within 0–1000 m for each eddy-day, retains surface ellipse geometry, and records
`w_surface`, `w_selected_depth_m`, `w_profile_n`. It is not a depth-following
footprint or necessarily the surface fitted vorticity. Verify the call before
writing methods text; some older chat descriptions assumed surface w.

Surface-cache semantics are unusually important:

- `use_cache=None`: calculate only, no cache write.
- `use_cache=True`: load the **entire** saved result; `df`/`grid` do not subset it.
  Validate settings and join current outcomes/properties by snapshot keys.
- `use_cache=False`: recalculate and replace the single configured cache file.
  Different settings do not automatically get separate filenames. Settings
  validation is not proof that underlying data are unchanged.
- The `surface_pv_esp_gaussian/` comparison series has its own versioned cache;
  it is distinct from the shared single-file surface cache.
- `source='depth_snapshot'` loads the thickness-weighted column-vector table;
  `source='depth'` loads per-depth rows. `depth_following=True` builds both using
  each level's own ellipse and vorticity, returning `(snapshot_df, depth_df)`.
  Surface-method upgrades must not silently change the depth-following estimand.

## Regimes, direction and shear

Common ratio: `log(|g_topo| / |g_plan|)`. Factor-two dominance corresponds to
thresholds ±log(2). Regimes normally classify **days**, not immutable whole tracks.
The current signed-PV hypothesis for saved deep-to-shallow tilt is **AE along
signed grad(PV), CE opposite**. It is a hypothesis/reference, not a forced sign.

Do not impose one universal selection on every workflow:

- Population dominance notebooks often use a centred 7-observation median
  (minimum 5), planetary depth ≥3000 m and directional tilt ≥5 km.
- Current topographic population/shear notebooks use depth ≤3000 m; paper case
  helpers include a 2000 m threshold. Read actual settings, including overrides.
- `esp_population_composites/pv_dominance_tilt.ipynb` uses **unsmoothed** component
  magnitudes, planetary `plan >= 2*topo` and `h > 3000`, topographic
  `topo >= 2*plan` and `h < 3000`. Mixed is the complement of both complete
  selections among valid days; exactly 3000 m is mixed. Invalid inputs are
  unclassified, not mixed. Selection uses saved `surface.TiltDis >= 5` and
  profile data are capped at 1000 m before averaging.
- Temporal `shear_tilt_climatology/` uses trailing consecutive-calendar-day
  smoothing and exact lag checks, avoiding future data in today's regime.
  Seven observations and seven consecutive days are not equivalent with gaps.

Background flow uses the v4 cache (`background_flow_cache_all_eddies_v4`):
**centred 91-day seasonal climatology**, formed by day-weighting the 26-year
calendar-month climatological fields overlapping the window. It is not a
91-day running mean of the event year's raw daily fields. Full-archive mean
is a sensitivity. The older annulus-background method is retired.

Primary shear/lift comparisons commonly use `U(0–200 m) - U(200–500 m)`;
derive the latter from `(500*U_0_500 - 200*U_0_200)/300` only with valid coverage.
Some control notebooks instead use surface minus 0–500 m: inspect the definition.
These velocity contrasts have m/s units, not s⁻¹ derivatives. Positive transverse
is **left** of upper-minus-deep flow; polarity is not embedded in that axis.
Physical sigma depths vary spatially; do not treat `z_r` as one fixed depth vector.

## Workflow map: read only the branch relevant to the task

| Location | Purpose and entry points |
|---|---|
| `seacofs_tilt_tools.py`, `ml_subsurface_tools.py` | Shared loaders, geography, PV, plotting; ML features/model helpers. Reuse established functions. |
| `census/` | `tilt_census.ipynb` for distributions; `tilt_stats.ipynb` plus `tilt_stats_tools.py` for paper numbers. |
| Root notebooks | Broad census/tilt/rose/PV/statistical/ML exploration. `Untitled.ipynb` is not an entry point. |
| `plan_topo_dom_eddies/` | Planetary/topographic populations, mixed cases, Rossby dependence, shelf-constraint hypothesis and sensitivities. |
| `surface_pv_footprint_sensitivity/` | Build with 00, then 01–06: legacy/nonlinear averaging, ellipse fractions, annuli, seamount cancellation and stability. |
| `surface_pv_esp_gaussian/` | Build corrected comparison cache with 00, then 01–06: reconstruction, weighting, uniform control, environmental/full PV, examples and assessment. |
| `depth_following_pv/` → `depth_resolved_pv_tilt/` | Build snapshot/depth caches, then 01–07: coverage, paired surface/depth comparison, regime switching, direction, magnitude, representation and mismatch cases. |
| `beta_effect_background_flow/` | Build background cache first; residual drift, vorticity-loss and planetary/topographic vertical-flow/shear notebooks. Prefer this maintained cached workflow over the older root background-relative notebook. |
| `tilt_mechanisms/` | 00 stratification cache; 01 shear, 02 beta/Burger, 03 PV/topography, 04 depth propagation, 05 optional wind, 06 joint tests, 07 temporal tilt, 08 stratification/polarity. |
| `shear_lift_tilt/` | `rotation_dependent_shear_lift.ipynb`: along/cross-shear geometry and polarity tests with existing background cache. |
| `shear_tilt_climatology/` | 00 inputs → 01 temporal climatology; 02 geographic shear-direction preference can reuse 00 without rerunning 01. |
| `eddy_env_tilt_controls/` | Broad `eddy_tilt_controls` and `environment_tilt_controls`; compact `paper_eddy_controls` and `paper_environment_controls` use median curves and IQR. Latest paper notebooks also include `TiltDis / surface Rc` figures. |
| `ellipse_tilt_analysis/` | Ellipse orientation/axis ratio versus tilt, depth comparisons and appended publication sections. Orientation is axial (180° symmetry); near-circular axes need QC. |
| `vertical_velocity/` | 01 individual model snapshots; 02 cache/climatology; 03 tilt relationship and 04 broad depth-change proxy reuse that cache. Not a closed stretching budget. |
| `delta_tilt_method/` | Schematic, temporal-kernel sensitivity and saved-tilt versus projected-extent comparison. Changes to production tilt belong upstream too. |
| `case_studies/` | 01 planetary, 02 topographic, 03 transitions; long-eddy/depth companions. Rankings select illustrations, not independent proof. |
| `case_studies/pv_tilt_transition_examples/` | 00 selected-eddy cache → 01 overviews / 02 selected-day core maps. Check arrow direction labels. |
| `case_studies/eddy_cross_sections_3d/` | `eddy_sections_and_esp`: native hydrographic/velocity sections and horizontal model-vs-ESP vorticity; `eddy_esp_composite`: reconstruct daily fields before averaging. |
| `esp_population_composites/` | Current entry points: `composite_breakdown`, `planetary_composite_tilt`, `topographic_composite_tilt`, `regional_composite_tilt`, `pv_relative_tilt`, `pv_dominance_tilt`. Inspect notebook/helper weighting and frames individually. |
| `eddy_categories/` | Exploratory categorisation; retain distinction from PV dominance and shared region labels. |

The population-composite README describes an earlier staged workflow and is not
a reliable inventory of current notebooks. Newer centreline notebooks can use
day-weighted means, while the older population helper supports equal-eddy field
composites. Do not transfer its support thresholds or weighting without checking.

## Recent decisions / issues to retain

- **Production delta depth support:** intervals must be supported on the reference
  day before temporal averaging, variance weighting and reconstruction. Neighbours
  cannot extend that day's depth range. Production remains centred 5 days,
  Gaussian sigma 1 day; the sensitivity notebook explicitly uses 7 days. Saved
  outputs need regeneration (`compute_tilt`, then `analyse_tilt`, with
  `parallel.skip_existing: false`); updating code alone does not change them.
- **N² correction:** v4 uses surface-referenced potential density, replacing the
  in-situ-density v3 calculation. Do not reuse v3. Rebuild with mechanism notebook
  00, then rerun 02, 03, 06, 08 and any other actual consumers. 01/04/05/07 do not
  need rerunning solely for that cache change. Keep version/method validation.
- **Latest paper controls:** commit `2353710` adds surface-radius-normalised
  versions alongside absolute-km figures. Normalisation uses the snapshot's
  **surface** `Rc`, including for subsurface-property comparisons. The subsequent
  `09d5dd9` update includes a 2×3 environmental figure: latitude, beta and log-PV
  versus raw/normalised tilt. Its linear fits/R² describe binned medians, not an
  adjusted individual-observation regression; the older README's no-fit
  description does not cover this addition.
- **Chat interpretation, not revalidated here:** `shear_tilt` reported CE core
  N² remaining higher after the correction. It distinguished absolute core
  stratification from environmental stratification; this contrast alone is not
  a reason to change the formula. Background N / core-minus-background anomalies
  were proposed, not established here as completed work.
- **Chat interpretation, not revalidated here:** `wholistic_pv_grad` found that
  CE right-of-shear preference was more persistent across geographic sectors,
  while AE left-of-shear preference varied/reversed and weakened with equal-sector
  weighting. Do not present universal AE-left/CE-right response as settled.
- **Evidence limits:** saved outputs may predate code/settings. During this review,
  saved error outputs existed in the older root background-relative notebook
  (`KeyboardInterrupt`) and case-study notebook 03 (`KeyError`). These are review
  leads, not newly reproduced bugs or instructions to fix unrelated notebooks.

## Editing, verification and handoff

Keep notebooks readable: question, editable controls, loading, concise analysis
and figures; put reusable calculations in focused Python helpers. Avoid large
unrelated notebook rewrites or clearing outputs the user supplied for review.
For scientific changes, test meaningful invariants with synthetic data (signs,
units, geometry, key joins, gaps, coverage, weighting and cache compatibility).
Use the relevant `test_*.py` suite, e.g. from the repository root:

```bash
MPLBACKEND=Agg python -m unittest discover \
  -s seacofs_eddy_tilt_analysis/eddy_env_tilt_controls -p 'test_*.py' -v
```

Replace the directory for the changed workflow; consult its README for extra
dependencies. Validate notebook JSON/cell syntax without running a multi-year
cache build. Distinguish syntax/unit/synthetic execution checks from real Katana
execution and scientific result validation. Do not fabricate production outputs.
Check diffs for large accidental output changes, then commit only intended files.
When handing off, state the commit/push result and exact rerun/cache requirements.

## How to refresh this guide

When asked to update it: compare commits since the reviewed revision, inspect
changed notebook settings/helpers/outputs and the relevant recent project tasks,
then **replace stale statements**. Update the date/revision and decisions above.
Keep this a compact navigation and decision record; do not append transcripts,
every plot statistic, transient local paths or an ever-growing commit diary.
For a numeric scientific conclusion, identify the notebook/output and cache/run
provenance; otherwise label it provisional. Keep planned work distinct from
implemented code and from executed/validated results.

Initial provenance: scanned current notebook introductions/settings/output-error
status and shared/helper code, read workflow READMEs and recent Git history, and
sampled recent accessible turns from these project tasks (not exhaustive history):
`eddy_dataset_modular`, `eddy_tilt_analysis`, `getting_to_the_bottom_of_eddy_tilt`,
`esp_vel_reco`, `beta_analysis`, `shear_tilt`, `wholistic_pv_grad`,
`plan_topo_dom_eddies`, `pv_grad_measuring`, `case_studies`,
`population_esp_composites`, `ellipse_orientation_tilt`, `vert_velocity`,
`delta_tilt_method`, `tilt_controls`, `text_eddy_dataset`, `aviso_eddy_dataset`.
Use those exact task titles for targeted retrieval if a decision needs history.
