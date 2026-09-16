# ESP population composites

First implementation: **AE/CE × planetary/topographic** composites with
complete-profile selection, equal-eddy weighting and eddy-cluster bootstrap
uncertainty. Full design and future extensions are in `DESIGN.md`; the repository
survey is in `NOTEBOOK_INVENTORY.md`.

## Run on Katana

Pull `UNSW-PhD/main`, then run these notebooks in order:

1. **00_population_audit.ipynb** — choose depths/thresholds, load existing tables,
   inspect exclusions/geography and freeze the selected Eddy-Day manifest.
2. **01_polarity_regime_composites.ipynb** — reconstruct daily ESP fields, save
   per-eddy means, bootstrap and plot maps, centrelines, sections and 3-D vectors.
3. **02_selection_sensitivity.ipynb** — inspect attrition and compare equal-eddy
   versus day-weighted centrelines.

Run from this directory or elsewhere inside `seacofs_eddy_tilt_analysis`.
Notebook 01 loads the last frozen run by default; `RUN_PATH` can select a specific
run. Each parameter/source combination has a separate output directory under:

`/srv/scratch/z5297792/SEACOFS_26yr_eddy_dataset_modular/esp_population_composites/`

Existing dependencies are numpy, pandas, scipy, matplotlib, netCDF4 and a Parquet
engine. The grid loader also needs the existing `z_r.npy`. ESP loads from
`/home/z5297792/ESP_zonodo/functions.py`; change `ESP_ROOT` in notebook 00 if needed.

Inputs are the processed eddy table, confirmed vertical profiles, model grid and
existing `source='depth_snapshot'` PV cache. These notebooks do **not** require
rebuilding the N²/background-flow caches or loading daily velocity NetCDF files.
No population composites have been evaluated on real Katana data locally.

## Defaults and interpretation

- Use cached levels nearest **0, 200 and 500 m**, within a maximum 60 m mismatch;
  actual levels are printed and stored. Adjust target depths to include additional
  levels for finer sections. All levels must be present and valid for every
  selected day; there is no vertical interpolation. All members share the same
  shallow reference depth. A shallower companion run tests deep-profile selection.
- Positive radius, finite required parameters, positive-definite ellipse matrix,
  consistent AE/CE rotation sign and unique keys are required. No tilt magnitude
  or preferred tilt direction is used to select members.
- Regimes use the depth-snapshot log topographic/planetary ratio and a **2:1**
  threshold. Mixed and missing regimes remain in the audit.
- Planetary composites use **east/north**; topographic composites use
  **onshore/alongshore**. Slope orientation uses the surface-core bathymetric
  gradient, with provisional minimum net slope `1e-4 m/m`, coherence `0.5` and
  valid-cell fraction `0.8`. No noisy centre-point fallback is used.
- Whole columns share one reference centre and rotation. Radius normalisation
  uses that reference Rc at all depths. Velocities remain in m/s.
- Ocean support uses nearest native-grid wet/bathymetric values. Land,
  outside-domain and below-seabed locations are missing, not zero. These fields
  are idealised reconstructions evaluated where the ocean exists, not observed
  native velocity composites. The shoreline mask is resolved at the native grid.
- Days are averaged within each eddy/group, then eddies receive equal weight.
  Daily velocity fields are streamed; per-eddy fields are cached. A group of
  per-eddy means is loaded for uncertainty, not every daily field.
- **500 bootstrap draws**, seed 731, resample whole eddy IDs globally. A shared
  eddy retains its covariance across populations. Component and centreline
  pointwise 95% intervals are saved. Means/intervals below **20 contributing
  eddies** are masked. The threshold is a reporting convention, not a power test.
- Grey centrelines show between-eddy variation, while bands describe uncertainty
  in the population mean. Speed plots show the **speed of the mean velocity**,
  not mean daily speed. Deep-offshore centre displacement corresponds to
  deep-to-shallow onshore tilt.

## Outputs and limitations

Each run saves all exclusion reasons, selected rows, exact depths, parameters,
source sizes/timestamps, implementation hash, ESP backend hash, member/day counts,
per-eddy fields and centrelines, bootstrap summary NPZ files and figures.
Source metadata fingerprints are not full checksums of large input datasets.
Reconstruction checks that inputs and the implementation have not changed since
selection. Rerunning reconstruction deliberately rebuilds member files; it does
not silently resume stale caches. It overwrites results inside the same run.

These are **geometry-screened composites**, not a fully skill-validated sample.
Processed surface parameters can contain gap-filling; their original provenance
is not stored in the final table, and upstream IDs are renumbered. Every audit
row records this uncertainty. Reconstruction-error validation needs native model
fields and should first use stratified members in the existing case-study
cross-section notebook. Bootstrap bands do not include fitting error, systematic
selection bias, or dependence between distinct nearby eddies. Do not interpret
pixelwise intervals as a simultaneous significance test across a volume.

Current primary groups have different environmental frames and eligible samples;
a raw difference between their grids is not a controlled mechanism test. The
sensitivity notebook compares weighting, but it does not claim to implement
matched geography/season, spatial-block bootstrap or field-difference intervals.
High/low Rossby, core/background stratification, shear-relative and other groups
remain subsequent extensions after the initial selection audit is reviewed.

## Local validation

```bash
MPLCONFIGDIR=/tmp/population_mpl python3 -m unittest discover \
  -s seacofs_eddy_tilt_analysis/esp_population_composites -p 'test_*.py' -v
```

Synthetic checks cover colliding Day values across eddies, unequal track lengths,
fixed-depth attrition, invalid geometry, polarity, zero tilt retention, coordinate
rotation, wet masks, cluster intervals, shared IDs and all plotting layouts.
The test velocity backend is analytic; it does not validate the external ESP
implementation or the scientific results on Katana.
