# IN2015_V03: Murphy and Freddy

Apply SOLO to each supplied crossing independently, with one notebook per eddy:
`murphy_solo.ipynb` and `freddy_solo.ipynb`. The common implementation lives in
`../esp_depth_tools.py`. LATTE is deferred.

Data path: `/srv/scratch/z5297792/tilt_observation_data/MNF_IN2015_V03_adcp_trawler/MNF_IN2015_V03_adcp_trawler.csv`

## Data and time coverage

- Requested export: 2015-06-04 00:00 UTC to 2015-06-12 00:00 UTC.
- CSV inspected on 2026-10-13: 83,332 rows; actual timestamps
  2015-06-04 00:01:22 UTC to 2015-06-11 21:21:43 UTC.
- Instrument/file: `in2015_v03_os75nb.nc`; 50 depth bins, 30.4-814.4 m,
  at 16 m spacing. Coverage varies by time and depth. The shallowest bin is
  not the sea surface and the deepest observation does not establish the eddy bottom.
- Columns used: `Start (UTC)`, `Mean Lon`, `Mean Lat`, `Depth`, `U (m/sec)`,
  `V (m/sec)`. Preserve SI units and UTC. The audit reports missing `Avg QC`
  and zero bin attendance; no quality status is inferred from blank fields.
- Eddy-level envelopes supplied by the user: Murphy, 5 June 02:00-17:30 UTC;
  Freddy, 8 June 22:45 to 10 June 03:15 UTC. These envelopes are not single
  straight transects and must not be pooled into one SOLO fit.

## Transects (single source for both notebooks)

Windows retain the user's supplied times, checked against the midtimes and
sampling durations in Roughan et al. (2017), Table 4 (PDF page 14, journal p.2507).
Coordinates below are the published MOCA centre estimates and serve only as
editable initial guesses / coordinate origins, not as constraints or validation
truth. Each crossing starts independently at its own origin.

| ID | Eddy | Start UTC | End UTC | Centre lon | Centre lat |
| --- | --- | --- | --- | --- | --- |
| M1 | Murphy | 2015-06-05T02:00:00Z | 2015-06-05T05:00:00Z | 154.8351 | -28.6168 |
| M2 | Murphy | 2015-06-05T14:30:00Z | 2015-06-05T17:30:00Z | 154.6091 | -28.6726 |
| F1 | Freddy | 2015-06-08T22:45:00Z | 2015-06-09T01:45:00Z | 153.2057 | -32.5763 |
| F2 | Freddy | 2015-06-09T03:00:00Z | 2015-06-09T04:30:00Z | 153.1891 | -32.5825 |
| F3 | Freddy | 2015-06-09T05:50:00Z | 2015-06-09T07:20:00Z | 153.2188 | -32.5691 |
| F4 | Freddy | 2015-06-09T13:45:00Z | 2015-06-09T14:45:00Z | 153.2487 | -32.6099 |
| F5 | Freddy | 2015-06-09T21:00:00Z | 2015-06-09T22:30:00Z | 153.2435 | -32.6550 |
| F6 | Freddy | 2015-06-10T02:15:00Z | 2015-06-10T03:15:00Z | 153.2461 | -32.6425 |

Table 4 also lists Freddy on 12 June at 16:15 UTC (1.5 h window, hence
15:30-17:00 UTC), plus two crossings on 16 June. These are outside this export
and are intentionally excluded. The first Freddy window crosses midnight.

## Method and interpretation

Reference: Roughan et al. (2017), *A tale of two eddies: The biophysical
characteristics of two contrasting cyclonic eddies in the East Australian
Current System*, JGR Oceans, [doi:10.1002/2016JC012241](https://doi.org/10.1002/2016JC012241).
Relevant locations: section 3.2 / Figure 3 (MOCA), Table 4 (times and centres),
section 5.2 / Figure 10 (vertical centre trajectories).

- Both eddies are cyclonic: use **negative relative vorticity** in this
  Southern Hemisphere coordinate convention.
- MOCA uses the normal-velocity zero crossing and locally uniform rotation to
  estimate centres. SOLO uses the ESP polynomial core fits and assumes circular
  geometry; the Gaussian outer stage is an additional model fit. These notebooks
  are a SOLO application, not an exact reproduction of the paper's MOCA estimates.
- The paper writes tangential speed as `V = Omega*r`. Its rotation-rate quantity
  must not be equated directly with the ESP relative vorticity: for SOLO,
  `zeta = 2*Omega`. The Gaussian `Rc` is also distinct from a reported eddy radius
  or diameter; its normalised velocity peak is at `Rc/sqrt(2)`.
- Time-window isolation is essential: the paper reports substantial translation
  between crossings. Profiles are not averaged across transects. Per-depth
  fitting still assumes approximately steady flow within each crossing.
- The 2015 notebooks explicitly select `solo_initialisation='local'`: the
  scaled, bounded crossing initialisation previously used in
  `UNSW-PhD/esp/local_estimators.py`. This retains the SOLO equations but avoids
  unconstrained root extrapolation. Set `'backend'` to use the original backend.
  The default for existing 2021 notebooks remains `'backend'`.
- Initial settings are editable, not tuned estimates of eddy size: Murphy uses
  a 25 km along-track core half-window; Freddy uses 10 km. Both require at least
  8 core/outer observations. The outer stage retains the original sign filter,
  but reports failures instead of silently returning the initial guess.
- Radius support is assessed separately from centre support. By default,
  Gaussian reconstructions require observations on both sides of `Rc/sqrt(2)`.
  Short near-core crossings can therefore yield centre trajectories without
  accepted full-field reconstructions. A finite fit or high R-squared is not
  proof of radius identifiability or accuracy.
- Centre offsets point from each crossing's shallowest accepted centre toward
  deeper centres. This differs from the model's saved deep-to-shallow delta tilt.
  Figure 10 should be read with its plotted vertical coordinate before comparing
  a sign or bearing; its centres move east/north with increasing depth.
- The paper's approximately 3.3 km/100 m (Murphy) and 1.5 km/100 m (Freddy)
  describe its inverse-variance-weighted MOCA tilt fits. They are contextual
  reference values, not targets used to select parameters here. These notebooks
  do not reproduce those confidence intervals or combine crossings into that fit.

## Notebook use

Run from this study directory. Edit controls at the top; both notebooks read
this table through `esp.load_study`. Each run fits all listed crossings for its
eddy and presents a comparison table/figure. `ACTIVE_TRANSECT` selects which
crossing is shown in detailed diagnostics, maps, 3D contours and fixed-plane
zonal/meridional sections. Change it and rerun the selection/plot cells to inspect
another crossing; there is no need to refit all crossings.

Saved outputs identify their execution environment and input fingerprints.
Optional CSV exports retain transect IDs, actual time bounds, rejected levels,
coordinates and settings. Do not treat numerical acceptance as independent
scientific validation. The shared-helper relocation preserves existing 2021
settings and saved outputs.
