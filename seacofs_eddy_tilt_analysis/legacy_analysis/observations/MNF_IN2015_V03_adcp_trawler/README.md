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

Run SOLO and export its results before LATTE: LATTE uses each selected SOLO crossing's shallowest accepted inner centre as its initial seed. Buffered LATTE uses all ADCP samples within each expanded UTC window, including off-line track positions, without a geographic mask. Edit `buffer_hours` by crossing in the buffered config; `fit_by_transect` contains a complete independent fit configuration for each crossing. Set `active_transect` there and rerun only the selected-crossing cells after the fits already exist in memory. For Freddy, `visit` selects F1–F7, F8–F9, or all F1–F9; the regional notebook fits one geographic footprint at a time. The buffered notebook compares every selected crossing before showing detailed diagnostics and reconstruction for `active_transect` at the bottom. No profile is averaged across dates or visits.

Results and provenance go to `/srv/scratch/z5297792/tilt_observation_data/MNF_IN2015_V03_adcp_trawler/outputs/`. The regional and buffered LATTE exports retain their existing distinct filenames. Rerun SOLO exports after editing `transects.json` or SOLO settings, then rerun the relevant LATTE notebook. Notebook outputs copied during the restructure are historical until rerun with the new config files.

At the bottom of both `freddy_solo.ipynb` and `murphy_solo.ipynb`, `along_track_bin_km` in the corresponding JSON config selects distance bins for a sensitivity comparison on the active crossing. It plots the original and binned centre offsets, Gaussian radius, and inner vorticity by depth; this comparison does not change the exported SOLO fits.

The four LATTE notebooks provide the same bin-width comparison using LATTE and their own fit settings (including per-crossing settings). Each config has `along_track_bin_km`. Buffered comparisons use only the selected crossing. Regional comparisons use the existing regional selection and bins projected onto its shallowest-bin principal axis; different ship passes and sampling times can contribute to one bin. These are sampling sensitivities, and do not replace exported fits.

## Independent fit settings and adaptive radius thresholds

SOLO and buffered LATTE use `fit_by_transect`, with complete blocks for M1–M2
or F1–F9. Freddy regional LATTE uses `fit_by_visit.visit1` and
`fit_by_visit.visit2`; Murphy regional LATTE uses one `fit` block. No shared
fit defaults or overrides need to be combined when editing these JSON files.
Existing numerical settings were copied into every corresponding block;
`radius_tracking` starts disabled. All distances are metres.

Each block includes these controls alongside the other fit settings:

```json
"core_radius": 30000,
"outer_radius": 70000,
"radius_tracking": true,
"core_radius_factor": null,
"outer_radius_factor": 1.75,
"initial_Rc": null
```

With tracking off, both explicit radii remain fixed. With tracking on, a null
factor keeps its threshold fixed, while a positive factor multiplies the last
accepted Gaussian `Rc`. In this example the inner threshold stays at 30 km and
the outer threshold starts at 70 km, then follows `1.75Rc`. Set the core factor
to `0.2` to make it adaptive too. A supplied `initial_Rc` initialises adaptive
thresholds from the first measured depth; otherwise the explicit radii initialise
them independently until the first accepted Gaussian radius. A fixed threshold
always uses its explicit radius, even when `initial_Rc` is supplied.

Both thresholds use the preceding accepted radius throughout the current depth,
including both LATTE core-selection iterations. A rejected Gaussian fit, rejected
inner fit, or missing depth leaves the radius reference unchanged. Accepted inner
centres and ellipses still update independently of Gaussian acceptance. LATTE
thresholds use ellipse radius `rho`; SOLO retains its along-track core selection.
No extra clipping or smoothing is applied to adaptive thresholds.

Exports include `core_radius_used_m`, `outer_radius_used_m`,
`radius_reference_m`, and `radius_reference_depth_m` for every depth. A reference
with no source depth denotes the configured initial radius; an absent reference
is NaN. Each crossing starts its own tracking state, and each binning-sensitivity
refit starts afresh with that application's settings. Settings are saved in
provenance and per-crossing audits.

After changing fit settings, rerun the notebook from the top; changing only the
active crossing selects a previously computed fit. Rerun SOLO exports before
LATTE when its SOLO seeds change. Existing notebook outputs are historical until
rerun. Synthetic checks (no voyage CSV or external ESP backend required):

```bash
MPLBACKEND=Agg python -m unittest discover \
  -s seacofs_eddy_tilt_analysis/observations -p 'test_radius_tracking.py' -v
```
