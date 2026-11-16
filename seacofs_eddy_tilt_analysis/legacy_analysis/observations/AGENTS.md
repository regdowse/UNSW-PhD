# Observational ESP work

Use this with `../AGENTS.md`. This directory applies SOLO and LATTE to voyage
observations; the SEACOFS model workflows have different inputs and conventions.

- Start with the relevant voyage's README or `info.md` and notebook controls.
  For IN2015_V03, the voyage README maps the notebooks and configs:
  `transects.json` defines crossing windows and origins, and each notebook's JSON
  defines its settings. Its `info.md` is study context, not a config.
- Shared methods are in `esp_depth_tools.py`; the Katana ESP backend is
  `/home/z5297792/ESP_zonodo/functions.py`. Raw ADCP, satellite caches, and
  exports live under `/srv/scratch/z5297792/tilt_observation_data/`. Verify current
  configs and provenance before interpreting saved notebook outputs.
- Use UTC and actual measured depth bins. Keep velocities in m/s and check the
  documented units of ESP positions and radii. The shallowest ADCP bin is not the
  sea surface; the deepest bin does not establish the eddy bottom.
- Preserve independent crossings, failed-fit gaps, and the distinction between
  accepted inner centres and accepted Gaussian reconstructions. Here CE fitted
  relative vorticity is negative and AE is positive.
- Check each notebook's selection: IN2015_V03 regional LATTE uses its geographic
  footprint without a time filter, whereas buffered LATTE uses all usable ADCP
  samples in each expanded crossing window. Other voyages use their own windows.
  If LATTE uses a SOLO seed, rerun the SOLO export after changing its crossings or
  fit settings, then rerun LATTE.
- Keep notebooks short; put shared calculations in `esp_depth_tools.py`. Preserve
  user edits and outputs. For code changes, validate configs and notebook syntax,
  inspect the diff, and run a representative Katana data check. Report what ran
  and which saved outputs need refreshing.
