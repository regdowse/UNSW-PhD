# IN2015_V03: Murphy and Freddy

Apply SOLO to each supplied crossing independently in `murphy_solo.ipynb` and
`freddy_solo.ipynb`. `murphy_latte.ipynb` and `freddy_latte.ipynb` each start
from a selected SOLO crossing's shallowest accepted inner centre, then fit all
ADCP observations in the corresponding regional subset at each depth; LATTE is
not restricted to the SOLO crossing windows. Freddy LATTE separates the early
(F1-F7) and later (F8-F9) geographic footprints, but applies no time filter
within the chosen footprint. Its early footprint contains both OS75 and OS150
depth grids from different dates, so a centreline can mix temporal and vertical
changes.
Both LATTE notebooks also have a separate buffered-transect analysis after the
regional export. It takes each crossing time from this table, adds editable hours
before and after, and fits all ADCP samples in that expanded window without a
straight-track or geographic mask. Each crossing starts from its own shallowest
accepted SOLO inner centre. The buffered fits have separate scratch exports;
rerun the SOLO export after editing this table so those centres stay aligned.
The common implementation lives in `../esp_depth_tools.py`. Notebook exports
go to an `outputs/` directory beside the CSV on scratch.

Data path: `/srv/scratch/z5297792/tilt_observation_data/MNF_IN2015_V03_adcp_trawler/MNF_IN2015_V03_adcp_trawler.csv`

## Data and time coverage

- Expanded export inspected on 2026-10-13: 112,412 rows; timestamps
  2015-06-04 00:01:22 UTC to 2015-06-17 00:57:38 UTC.
- OS75 (`in2015_v03_os75nb.nc`): 83,332 rows; 50 depth bins from 30.4
  to 814.4 m at 16 m spacing. OS150 (`in2015_v03_os150nb.nc`): 29,080
  rows; 40 depth bins from 18.5 to 330.5 m at 8 m spacing. The available
  bins within an individual crossing can span less than the instrument range.
  The shallowest bin is not the sea surface and the deepest observation does
  not establish the eddy bottom.
- Columns used: `Start (UTC)`, `Mean Lon`, `Mean Lat`, `Depth`, `U (m/sec)`,
  `V (m/sec)`. Preserve SI units and UTC. The audit reports missing `Avg QC`
  and zero bin attendance; no quality status is inferred from blank fields.
- The original eddy envelopes supplied by the user are Murphy, 5 June
  02:00-17:30 UTC, and Freddy visit 1, 8 June 22:45 to 10 June 03:15 UTC.
  The expanded export now also covers Freddy on 12 and 16 June. Each window
  below is fitted separately; the broader envelopes are not straight transects.

## Transects (single source for both notebooks)

The original windows retain the user's supplied times. F7-F9 follow the
midtimes and sampling durations in Roughan et al. (2017), Table 4 (PDF page 14,
journal p.2507).
Coordinates below are the published MOCA centre estimates and serve only as
editable initial guesses / coordinate origins, not as constraints or validation
truth. Each crossing starts independently at its own origin.

| ID | Eddy | Start UTC | End UTC | Centre lon | Centre lat |
| --- | --- | --- | --- | --- | --- |
| M1 | Murphy | 2015-06-04T22:00:00Z | 2015-06-05T07:00:00Z | 154.8351 | -28.6168 |
| M2 | Murphy | 2015-06-05T14:30:00Z | 2015-06-05T17:30:00Z | 154.6091 | -28.6726 |
| F1 | Freddy | 2015-06-08T22:15:00Z | 2015-06-09T01:55:00Z | 153.2057 | -32.5763 |
| F2 | Freddy | 2015-06-09T03:00:00Z | 2015-06-09T04:30:00Z | 153.1891 | -32.5825 |
| F3 | Freddy | 2015-06-09T05:50:00Z | 2015-06-09T07:20:00Z | 153.2188 | -32.5691 |
| F4 | Freddy | 2015-06-09T13:45:00Z | 2015-06-09T14:45:00Z | 153.2487 | -32.6099 |
| F5 | Freddy | 2015-06-09T21:00:00Z | 2015-06-09T22:30:00Z | 153.2435 | -32.6550 |
| F6 | Freddy | 2015-06-10T02:15:00Z | 2015-06-10T03:15:00Z | 153.2461 | -32.6425 |
| F7 | Freddy | 2015-06-12T15:30:00Z | 2015-06-12T17:00:00Z | 153.4448 | -32.8224 |
| F8 | Freddy | 2015-06-16T15:30:00Z | 2015-06-16T16:30:00Z | 153.9855 | -33.2769 |
| F9 | Freddy | 2015-06-16T23:00:00Z | 2015-06-17T00:00:00Z | 153.8311 | -33.2702 |

<!-- Original -->
<!-- | ID | Eddy | Start UTC | End UTC | Centre lon | Centre lat |
| --- | --- | --- | --- | --- | --- |
| M1 | Murphy | 2015-06-05T02:00:00Z | 2015-06-05T05:00:00Z | 154.8351 | -28.6168 |
| M2 | Murphy | 2015-06-05T14:30:00Z | 2015-06-05T17:30:00Z | 154.6091 | -28.6726 |
| F1 | Freddy | 2015-06-08T22:45:00Z | 2015-06-09T01:45:00Z | 153.2057 | -32.5763 |
| F2 | Freddy | 2015-06-09T03:00:00Z | 2015-06-09T04:30:00Z | 153.1891 | -32.5825 |
| F3 | Freddy | 2015-06-09T05:50:00Z | 2015-06-09T07:20:00Z | 153.2188 | -32.5691 |
| F4 | Freddy | 2015-06-09T13:45:00Z | 2015-06-09T14:45:00Z | 153.2487 | -32.6099 |
| F5 | Freddy | 2015-06-09T21:00:00Z | 2015-06-09T22:30:00Z | 153.2435 | -32.6550 |
| F6 | Freddy | 2015-06-10T02:15:00Z | 2015-06-10T03:15:00Z | 153.2461 | -32.6425 |
| F7 | Freddy | 2015-06-12T15:30:00Z | 2015-06-12T17:00:00Z | 153.4448 | -32.8224 |
| F8 | Freddy | 2015-06-16T15:30:00Z | 2015-06-16T16:30:00Z | 153.9855 | -33.2769 |
| F9 | Freddy | 2015-06-16T23:00:00Z | 2015-06-17T00:00:00Z | 153.8311 | -33.2702 | -->

F7 is the 12 June continuation of Freddy visit 1; F8 and F9 are visit 2
on 16 June. All three are now present in the export and use the OS150
instrument. The first and last Freddy windows cross midnight. Table 4
midtimes/durations define these windows; the new CSV has 18, 12 and 12
observation times within F7, F8 and F9 respectively. Within the windows,
available depth levels span 34.5-250.5 m (F7), 18.5-290.5 m (F8), and
18.5-290.5 m (F9). Keep the visits distinct when interpreting centre
movement or vertical structure.

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



