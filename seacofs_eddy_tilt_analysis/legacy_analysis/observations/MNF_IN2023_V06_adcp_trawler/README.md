# CE23 — IN2023_V06 ADCP SOLO and LATTE

All times are UTC and selection endpoints are inclusive. Data remain at the
scratch path in `transects.json`; no raw CSV is copied into the repository.

1. Run `CE23_solo.ipynb`, controlled by `CE23_solo.json`. C1 covers
   **17 October 2023 13:30 → 18 October 2023 13:30 UTC**, with initial origin
   **154.3°E, 37.9°S**. Its export supplies LATTE's initial centre.
2. Run `CE23_latte.ipynb`, controlled by `CE23_latte.json`. It selects all usable
   observations during **17 October 2023 10:27:31 → 21 October 2023
   22:27:30.333 UTC**, with no geographic filter. Extra CSV data on either side
   are excluded. `solo_transect` selects which SOLO centre seeds LATTE.

To add a SOLO crossing, add a unique `id`, `eddy: "CE23"`, UTC `start`/`end`,
`lon` and `lat` to `transects.json`, then copy a complete `fit_by_transect`
block under that ID in `CE23_solo.json` and tune it. All configured crossings
are fitted independently; `active_transect` selects the detailed plots.
The SOLO study-window control bounds the overview and permitted crossings.

Fit controls match the 2015 notebooks: independent settings, optional adaptive
inner/outer radii, acceptance checks, and along-track-bin sensitivity. Settings
are copied from Murphy as starting values and have not been optimised for CE23.
No background subtraction, motion correction or vertical interpolation is added.
Projected binning across the full survey can combine separate ship passes.

Exports are `ce23_solo_depth_fits.csv`, `ce23_solo_provenance.json`,
`ce23_latte_depth_fits.csv` and `ce23_latte_provenance.json` in the input CSV's
`outputs/` directory. After changing SOLO settings or crossing definitions,
rerun SOLO with `save_results: true` before LATTE. Changing LATTE's time window
or fit settings requires rerunning LATTE. Satellite overlays in SOLO and LATTE are optional
and may download daily maps; they do not affect the fit or export.

The notebooks start without saved outputs. “Surface” means the shallowest
observed bin. Cyclonic relative vorticity is negative in this workflow.

Select `adcp_instrument` in each JSON: `"OS150"` (default) or `"WH75"`.
Overview plots and fits use only that source, with independent depth tracking.
Set both configs to the same instrument and rerun SOLO with exports, then LATTE.
Selection audits record the chosen instrument and source filename. Existing saved
notebook outputs and scratch exports need refreshing after changing this setting.
