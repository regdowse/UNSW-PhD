# Delta tilt method

Open `delta_tilt_schematic.ipynb` on Katana and run all cells. It loads the actual confirmed vertical profiles from `/srv/scratch/z5297792/SEACOFS_26yr_eddy_dataset_modular/vertical_profiles_confirmed/profiles.parquet` and randomly chooses a real eddy-day with five consecutive usable profiles and a finite, nonzero production tilt. No synthetic fallback is used.

The production method now uses a centred five-day window and Gaussian temporal sigma of one day. The schematic uses these defaults; depth weights remain inverse unweighted sample variance across the five days.

- Change `SEED` to try another example (default 731).
- Set `EDDY_ID` to choose an eddy, leaving `REFERENCE_DAY=None` for a random eligible day.
- Set both controls to reproduce a particular example. The notebook prints the selected values.

Selection is screened for illustration, not an unbiased population sample. Run from this folder or elsewhere inside UNSW-PhD. Dependencies are NumPy, pandas, Matplotlib and a Parquet engine. Defaults match the pipeline tilt configuration.

Panel a shows daily centre profiles; b shows the mean-increment reconstruction, variance-weighted line and deep-to-shallow direction; c compares fitted tilt distance with the reference day's projected horizontal extent. The ribbon represents relative fit weight, not uncertainty. The illustrated fit is checked against the production `compute_weighted_tilt` function.

PNG and PDF figures are exported here with the selected Eddy ID and Day in their names when `SAVE=True`. See the notebook for exact depth indexing and coordinate conventions. Real-data execution requires the Katana profile file.

## Temporal-weight sensitivity

Run `delta_temporal_weight_sensitivity.ipynb` on Katana to compare equal, Gaussian 2-day and Gaussian 1.5-day temporal means on 10 random eligible real eddy-days. Change `SEED` or restrict `EDDY_IDS`. It overlays three-panel snapshot schematics and plots the three delta estimates plus maximum pairwise horizontal centre separation over each sampled eddy's available lifetime. The existing unweighted depth-variance weights are held fixed. Equal-weight lifetime estimates are checked against production.

The small `delta_sensitivity_tools.py` helper keeps the notebook concise. SciPy supplies pairwise distances. Figures/CSV files are optional (`SAVE=False` by default); generated outputs are ignored by Git. This experiment retains its seven-day windows and explicitly checks equal weighting against the legacy settings (`num=6, temporal_sigma_days=None`). It does not regenerate the production dataset.
