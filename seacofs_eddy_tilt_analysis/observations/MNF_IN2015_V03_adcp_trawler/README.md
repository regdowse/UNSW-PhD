# IN2015_V03 ADCP ESP notebooks

Active crossing windows, coordinate origins, and the ADCP CSV path live in [`transects.json`](transects.json). The dated crossing definitions are analysis choices; [`info.md`](info.md) contains only context from Roughan et al. (2017). All application-specific settings are in the JSON file with the same stem as its notebook. Edit JSON, then rerun the corresponding notebook from the top. The shared methods live in `../esp_depth_tools.py`.

| Notebook | Application | Config |
| --- | --- | --- |
| `murphy_solo.ipynb` | Independent SOLO fits of M1–M2 | `murphy_solo.json` |
| `freddy_solo.ipynb` | Independent SOLO fits of selected F1–F9 visits | `freddy_solo.json` |
| `murphy_latte.ipynb` | All ADCP data inside the Murphy regional footprint | `murphy_latte.json` |
| `freddy_latte.ipynb` | All ADCP data inside one Freddy geographic visit footprint | `freddy_latte.json` |
| `murphy_latte_buffered.ipynb` | Independent LATTE fits of buffered M1–M2 windows | `murphy_latte_buffered.json` |
| `freddy_latte_buffered.ipynb` | Independent LATTE fits of buffered F1–F9 windows | `freddy_latte_buffered.json` |

Run SOLO and export its results before LATTE: LATTE uses each selected SOLO crossing's shallowest accepted inner centre as its initial seed. Buffered LATTE uses all ADCP samples within each expanded UTC window, including off-line track positions, without a geographic mask. Edit `buffer_hours` by crossing in the buffered config; optional `fit_overrides` can change a named crossing’s LATTE settings without affecting the others. Set `active_transect` there and rerun only the selected-crossing cells after the fits already exist in memory. For Freddy, `visit` selects F1–F7, F8–F9, or all F1–F9; the regional notebook fits one geographic footprint at a time. The buffered notebook compares every selected crossing before showing detailed diagnostics and reconstruction for `active_transect` at the bottom. No profile is averaged across dates or visits.

Results and provenance go to `/srv/scratch/z5297792/tilt_observation_data/MNF_IN2015_V03_adcp_trawler/outputs/`. The regional and buffered LATTE exports retain their existing distinct filenames. Rerun SOLO exports after editing `transects.json` or SOLO settings, then rerun the relevant LATTE notebook. Notebook outputs copied during the restructure are historical until rerun with the new config files.

At the bottom of both `freddy_solo.ipynb` and `murphy_solo.ipynb`, `along_track_bin_km` in the corresponding JSON config selects distance bins for a sensitivity comparison on the active crossing. It plots the original and binned centre offsets, Gaussian radius, and inner vorticity by depth; this comparison does not change the exported SOLO fits.

The four LATTE notebooks provide the same bin-width comparison using LATTE and their own fit settings (including buffered-crossing overrides). Each config has `along_track_bin_km`. Buffered comparisons use only the selected crossing. Regional comparisons use the existing regional selection and bins projected onto its shallowest-bin principal axis; different ship passes and sampling times can contribute to one bin. These are sampling sensitivities, and do not replace exported fits.
