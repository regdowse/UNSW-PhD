# Delta tilt method

Open `delta_tilt_schematic.ipynb` on Katana and run all cells. It loads the actual confirmed vertical profiles from `/srv/scratch/z5297792/SEACOFS_26yr_eddy_dataset_modular/vertical_profiles_confirmed/profiles.parquet` and randomly chooses `N_IMAGES=10` unique real eddy-days with five consecutive usable profiles and a finite, nonzero production tilt. No synthetic fallback is used.

The production method now uses a centred five-day window and Gaussian temporal sigma of one day. The schematic uses these defaults; depth weights remain inverse unweighted sample variance across the five days.

- Set `N_IMAGES` to the desired batch size (default 10).
- Change `SEED` to try another batch (default 46). Eddies are visited in random order once per round before selecting more days from them.
- Set `EDDY_ID` to choose an eddy, leaving `REFERENCE_DAY=None` for a batch of random eligible days from that eddy.
- Set both controls to reproduce a particular example. This produces one figure regardless of `N_IMAGES`. The notebook prints a selection table and checks each result against production.

Selection is screened for illustration, not an unbiased population sample. Run from this folder or elsewhere inside UNSW-PhD. Dependencies are NumPy, pandas, Matplotlib and a Parquet engine. Defaults match the pipeline tilt configuration.

Panel a shows daily centre profiles; b shows the mean-increment reconstruction, variance-weighted line and deep-to-shallow direction; c compares fitted tilt distance with the reference day's projected horizontal extent. The ribbon represents relative fit weight, not uncertainty. The illustrated fit is checked against the production `compute_weighted_tilt` function.

PNG and PDF figures are exported here with the selected Eddy ID and Day in their names when `SAVE=True`, along with `delta_tilt_schematic_selected_examples.csv`. By default figures display inline without saving. If fewer eligible pairs exist than requested, all available pairs are shown with a message. See the notebook for exact depth indexing and coordinate conventions. Real-data execution requires the Katana profile file.

## Temporal-weight sensitivity

Run `delta_temporal_weight_sensitivity.ipynb` on Katana to compare equal, Gaussian 2-day and Gaussian 1.5-day temporal means on 10 random eligible real eddy-days. Change `SEED` or restrict `EDDY_IDS`. It overlays three-panel snapshot schematics and plots the three delta estimates plus maximum pairwise horizontal centre separation over each sampled eddy's available lifetime. The existing unweighted depth-variance weights are held fixed. Equal-weight lifetime estimates are checked against production.

The small `delta_sensitivity_tools.py` helper keeps the notebook concise. SciPy supplies pairwise distances. Figures/CSV files are optional (`SAVE=False` by default); generated outputs are ignored by Git. This experiment retains its seven-day windows and explicitly checks equal weighting against the legacy settings (`num=6, temporal_sigma_days=None`). It does not regenerate the production dataset.
