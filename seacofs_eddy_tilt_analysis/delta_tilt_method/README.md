# Delta tilt method

Open `delta_tilt_schematic.ipynb` on Katana and run all cells. It loads the actual confirmed vertical profiles from `/srv/scratch/z5297792/SEACOFS_26yr_eddy_dataset_modular/vertical_profiles_confirmed/profiles.parquet` and randomly chooses `N_IMAGES=10` unique real eddy-days with five consecutive usable profiles and a finite, nonzero depth-limited tilt. No synthetic fallback is used.

The schematic retains a centred five-day window and Gaussian temporal sigma of one day. The schematic uses these defaults; depth weights remain inverse unweighted sample variance across the five days.

- Set `N_IMAGES` to the desired batch size (default 10).
- Change `SEED` to try another batch (default 46). Eddies are visited in random order once per round before selecting more days from them.
- Set `EDDY_ID` to choose an eddy, leaving `REFERENCE_DAY=None` for a batch of random eligible days from that eddy.
- Set both controls to reproduce a particular example. This produces one figure regardless of `N_IMAGES`. The notebook prints a selection table and reports each reference profile depth and fitted depth.

Selection is screened for illustration, not an unbiased population sample. Run from this folder or elsewhere inside UNSW-PhD. Dependencies are NumPy, pandas, Matplotlib and a Parquet engine. Defaults match the pipeline tilt configuration.

Panel a shows daily centre profiles; b shows the mean-increment reconstruction, variance-weighted line and deep-to-shallow direction; c compares fitted tilt distance with the reference day's projected horizontal extent. The ribbon represents relative fit weight, not uncertainty. The schematic and sensitivity notebook share the local `fit_snapshot` helper. Production now applies the same reference-day support restriction.

PNG and PDF figures are exported here with the selected Eddy ID and Day in their names when `SAVE=True`, along with `delta_tilt_schematic_selected_examples.csv`. By default figures display inline without saving. If fewer eligible pairs exist than requested, all available pairs are shown with a message. See the notebook for exact depth indexing and coordinate conventions. Real-data execution requires the Katana profile file.

## Temporal-weight sensitivity

Run `delta_temporal_weight_sensitivity.ipynb` on Katana to compare equal, Gaussian 1.5-day and Gaussian 1-day temporal means on `N_SNAPSHOTS` random eligible real eddy-days. Change `SEED` or restrict `EDDY_IDS`. It overlays three-panel snapshot schematics and plots the three delta estimates plus maximum pairwise horizontal centre separation over each sampled eddy's available lifetime. The existing unweighted depth-variance weights are held fixed. All three methods use the reference day’s supported depth intervals.

The small `delta_sensitivity_tools.py` helper keeps the notebook concise. SciPy supplies pairwise distances. Figures/CSV files are optional (`SAVE=False` by default); generated outputs are ignored by Git. This experiment retains its seven-day windows but now limits every estimate to its reference day’s depth support. It does not regenerate the production dataset.

## Reference-day depth limit

Both notebooks now exclude intervals not supported by the reference day **before** temporal averaging, depth-variance weighting and cumulative reconstruction. Neighbouring profiles still contribute inside that range, but cannot extend it above/below the reference profile. No reference profile means no estimate. All temporal kernels use the same depth restriction; lifetime curves apply it separately each day.

The interpolation and upper-interval labels are retained: a profile from 0 to 200 m supplies intervals labelled 0 to 190 m. It therefore fails the unchanged 200 m minimum fit-label range and yields no estimate. The maximum-depth control remains an additional cap. Existing interpolation between fitted levels is retained, with no extrapolation outside the measured range.

The main pipeline in `seacofs_eddy_dataset_modular` now uses this restriction too. Tests compare production against the notebook helper with varying depth ranges, shallow references and missing days. Regenerate the saved climatology by rerunning `compute_tilt` and then `analyse_tilt` with `parallel.skip_existing: false`; existing files do not change merely by updating the code. The production default remains five days with Gaussian sigma = one day, while the sensitivity experiment retains its explicitly configured seven-day windows.

## Largest schematic differences

At the top of `delta_tilt_schematic.ipynb`, set `RUN_DIFFERENCE_SEARCH=True` to enable the expensive bottom section and `TOP_N_DIFFERENCES=10` (or another positive integer). This scans all eddy-days in the loaded profile dataset independently of the random-example filters. It ranks `abs(MaxProjectedExtent - TiltDis)` in kilometres using each day's own fitted horizontal axis and the same reference-depth limit as the schematic. Every valid nonzero fit is eligible, including incomplete neighbour windows when sufficient data remain.

The ranking table includes signed and absolute differences, eddy/day IDs, and reference/fit depths. The final cell renders the requested top N using the existing plotting function. The scan runs once per settings choice; rerun only the plotting cell after changing top N. Progress is printed every 100 eddies. `SAVE=True` also exports the full ranking CSV and the usual PNG/PDF schematics; no additional cache is created.
