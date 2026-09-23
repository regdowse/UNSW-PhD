# Delta tilt method

Open `delta_tilt_schematic.ipynb` on Katana and run all cells. It loads the actual confirmed vertical profiles from `/srv/scratch/z5297792/SEACOFS_26yr_eddy_dataset_modular/vertical_profiles_confirmed/profiles.parquet` and randomly chooses a real eddy-day with seven consecutive usable profiles and a finite, nonzero production tilt. No synthetic fallback is used.

- Change `SEED` to try another example (default 731).
- Set `EDDY_ID` to choose an eddy, leaving `REFERENCE_DAY=None` for a random eligible day.
- Set both controls to reproduce a particular example. The notebook prints the selected values.

Selection is screened for illustration, not an unbiased population sample. Run from this folder or elsewhere inside UNSW-PhD. Dependencies are NumPy, pandas, Matplotlib and a Parquet engine. Defaults match the pipeline tilt configuration.

Panel a shows daily centre profiles; b shows the mean-increment reconstruction, variance-weighted line and deep-to-shallow direction; c compares fitted tilt distance with the reference day's projected horizontal extent. The ribbon represents relative fit weight, not uncertainty. The illustrated fit is checked against the production `compute_weighted_tilt` function.

PNG and PDF figures are exported here with the selected Eddy ID and Day in their names when `SAVE=True`. See the notebook for exact depth indexing and coordinate conventions. Real-data execution requires the Katana profile file.
