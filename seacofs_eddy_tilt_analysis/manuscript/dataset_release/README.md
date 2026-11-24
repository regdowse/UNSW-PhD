# SEACOFS mesoscale eddy dataset

Dataset version 1.0.0.

**Creator:** Reg W. Dowse (sole dataset creator).  
**Licence:** [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).  
**Dataset DOI:** [10.5281/zenodo.24486070](https://doi.org/10.5281/zenodo.24486070).

This dataset contains velocity-defined mesoscale eddy detections, tracked and
processed surface fits, fitted vertical centre profiles, and Delta tilt estimates
for the East Australian Current region. It was extracted from daily SEACOFS model
velocity fields covering **2 January 1994–28 February 2019** (9,189 calendar days).
The final surface dataset contains **3,097 eddies**: 1,497 anticyclonic eddies (AEs)
and 1,600 cyclonic eddies (CEs), contributing 133,038 processed eddy-days.

The release describes modelled eddies, not independent ocean observations. Each
eddy-day is a daily tracked eddy record; depth-resolved tables contain multiple
rows per eddy-day. Intermediate stages are included to document extraction.
Model velocity fields, manuscript plotting code and manuscript-specific PV-gradient
products are separate resources and are not included here.

## ESP/DOPPIO and Delta: methods and references

**Elliptical Streamfunction Parametrisation (ESP)** describes eddy velocity
structure using an elliptical streamfunction model. Its **DOPPIO** method fits
two orthogonal velocity transects to estimate the eddy centre, inner-core
rotation and elliptical shape. In this pipeline, the ESP `doppio` function is
used at the surface and successively at depth; the separate
`out_core_param_fit` function fits the outer Gaussian parameters.
The [ESP methods paper](https://doi.org/10.5194/egusphere-2026-3599)
describes the framework and DOPPIO. The
[ESP v1.1 software DOI](https://doi.org/10.5281/zenodo.24182671)
identifies the implementation used in the recorded generation run.

**Delta** is a separate tilt estimator applied to the resulting fitted centre
profiles. It combines adjacent-depth centre displacements across neighbouring
days, weights their temporal consistency, and fits a representative axis to
obtain one tilt distance and direction per supported eddy-day. Delta is described
in the accompanying eddy-tilt manuscript: **manuscript DOI: XXXX**
(placeholder; DOI not yet assigned). The extraction and tilt settings needed to
interpret this dataset are also documented below. Its implementation is archived
in [Ocean Eddy Pipeline v0.1.1](https://doi.org/10.5281/zenodo.24457378),
in `src/seacofs_eddy_dataset/core/tilt.py`.

The ESP methods-paper DOI, ESP software DOI, pipeline DOI and eventual Delta
manuscript DOI identify different resources. Replace `XXXX` with the manuscript
DOI when available; it is not an active DOI or a link to the ESP paper.

## Files and starting points

| File or directory | Contents | Rows |
|---|---|---:|
| `processed/eddy_dataset_processed.parquet` | Final processed surface eddy-days; start here | 133,038 |
| `vertical_profiles_confirmed/profiles.parquet` | Profiles passing the recorded-level-count QC | 2,994,569 |
| `tilt/tilt_dataset.parquet` | Delta evaluation days, including missing estimates | 120,348 |
| `analysis/tilt_summary.parquet` | One-row summary of valid Delta estimates | 1 |
| `analysis/profile_tilt_summaries.parquet` | Raw surface-relative profile geometry; **not Delta** | 126,681 |
| `detections/` | Initial surface vector-geometry candidates | 414,949 |
| `surface_eddies/` | Surface-fit attempts, including failed fits | 414,949 |
| `tracked/eddy_tracks.parquet` | Cleaned and tracked observations before final processing | 260,761 |
| `vertical_profiles/` | Fitted profile records before count-based QC | 3,003,294 |
| `z_r.npy` | Static physical depths for native model levels, in analysis-array order | (272, 317, 30) |
| `static_grid.nc` | Static geometry: model coordinates, longitude/latitude, bathymetry, mask, angles and target depths | See below |
| `variables.csv` | Field-by-field definitions, types, units and qualifications | — |
| `generation_config.yaml` | Configuration distributed with pipeline v0.1.1 | — |
| `CITATION.cff` | Machine-readable dataset creator, version and licence | — |
| `LICENSE.txt` | Reuse licence and source-model attribution | — |

Each of `detections/`, `surface_eddies/` and `vertical_profiles/` contains 307
`fnumber=NNNNN.parquet` partitions, one per source model file. Partition labels
identify model files, not eddies or single days. The last source file contains
nine daily records; preceding files contain 30.

## Identifiers, dates and joins

- Within the **final** products, `(Eddy, Day)` identifies an eddy-day. Join the
  processed, confirmed-profile, tilt and profile-summary tables using both fields.
  `(Eddy, Day, Depth)` is unique in the profile tables. A profile-to-surface join
  is many-to-one; do not join on `Day` alone.
- Final `Eddy` IDs are renumbered after filtering. They **do not preserve the
  tracked-stage ID namespace**. Detection/surface-fit tables precede final IDs.
  `eddy_idx` is a within-day tracking identifier, not a final track ID.
- `Day` is an integer number of days since 1990-01-01; `Date` is the corresponding
  timestamp. The processed range is Day 1462–10650, inclusive.
- `Age` is the **total number of processed daily records in the track**, repeated
  on each row. It is not elapsed age on that day. Elapsed age can be calculated
  from `Day - min(Day)` within a final `Eddy`.
- `fname` retains the original machine-local model-file path. Use its basename
  and the source-model archive to locate a file on another installation.

For example, after extracting the release:

```python
import pandas as pd
surface = pd.read_parquet('processed/eddy_dataset_processed.parquet')
tilt = pd.read_parquet('tilt/tilt_dataset.parquet')
daily = surface.merge(tilt, on=['Eddy', 'Day'], how='left', validate='one_to_one')
```

## Coordinates and units

Horizontal fit coordinates and radii are in **kilometres**. `lon` and `lat` use
signed geographic degrees east and north, respectively; negative latitude means
south. The analysis grid transposes source horizontal arrays into `(x, y)` order:
`ic` indexes x and `jc` indexes y, both zero-based. These are original-grid indices,
including `nic/njc` obtained by mapping initial detection centres to that grid.

`static_grid.nc` preserves the pipeline's `x` and `y` coordinate vectors and
transposed `lon_rho`, `lat_rho`, `h`, `mask_rho` and `angle` arrays. Cartesian
coordinate vectors are cumulative spherical distances along the middle grid row
and column (Earth radius 6,357 km); they are not a standard map projection.
`pipeline_reference_angle` is the single angle used by extraction, approximately
-20 degrees. Full model `angle` is also supplied for context. With this angle a,
the pipeline rotates input velocities as `u_model = u*cos(a) + v*sin(a)` and
`v_model = v*cos(a) - u*sin(a)`.

`z_r.npy` stores signed, negative-below-sea-level depths in metres, ordered from
the uppermost native level downwards. It aligns with `(x, y, native_level)`; the
first native level is below the surface, rather than zero. The profile **target**
depths are `abs(z_r[150,150,1:])`, common physical reference levels across the
horizontal domain. Profiles add a fitted surface at `Depth=0`. Its input velocities come from the
uppermost native model level; this is a surface label, not an extrapolation of
velocity to zero depth. `Depth` is positive
metres downward; `z` is a reference-level index, not the local native sigma index.
All 30 possible profile labels (surface plus 29 target depths) are supplied as
`profile_depth` in `static_grid.nc`. Available profiles extend to about 4,749 m.
The 1,000 m cap below applies to Delta, not to the released full profiles.

Stored `Omega0`, `Omega` and `w` are in s⁻¹ after the pipeline's factor-10⁻³
conversion from fits using kilometre coordinates and m/s velocities. `w` denotes
**relative vorticity**, not vertical water velocity. In this Southern Hemisphere
application, AE relative vorticity is positive and CE relative vorticity negative.
The Gaussian streamfunction parameter `psi0` retains the implementation's mixed
units **m km s⁻¹**; multiply it by 1,000 for m²/s. `R` is the mean available
four-direction peak-speed radius; `Rc` is the Gaussian scale. They are distinct.
Shape-matrix components and aspect ratio are dimensionless.

## Construction workflow

### 1. Surface detection

The vector-geometry method of Nencioli et al. identifies surface velocity
reversals, a local speed minimum and consistent circulation. The pipeline rotates
velocities into model coordinates, interpolates them onto a 1 km regular grid,
and uses windows `a=4`, `b=3`. Candidate centres are mapped to the nearest original
grid cells before fitting. See [Nencioli et al. (2010)](https://doi.org/10.1175/2009JTECHO725.1).

### 2. Surface ESP fitting and acceptance

DOPPIO fits two orthogonal transects extending 30 km each side of the initial
centre. It estimates centre position, inner-core angular velocity `Omega0`,
relative vorticity and symmetric shape matrix Q. For r relative to the fitted
centre, the elliptical coordinate is `rho² = rᵀ Q r`; the aspect ratio is
`sqrt(lambda_max/lambda_min)`. The smaller eigenvalue defines the major axis.

Directional radii are measured along four model-grid directions to the first
local decrease following the tangential-speed maximum; their available mean is
R. Outer Gaussian fitting uses `rho <= max(min(1.75*R, 200), 30)` km and a
rectangular extraction half-width three times this selection radius. This is an
observation-domain limit, not an optimisation bound on Rc. The outer fit returns
Rc, psi0 and Omega. Surface continuation requires compatible DOPPIO centres
through the first reference depth at or below 50 m (about 60 m here), with the
same circulation sign and successive jumps no larger than 100 km.

Surface fitting writes a row for each attempted candidate, retaining missing
parameters on failure. Before tracking, candidates require consistent detection
and DOPPIO cyclonicity, an interior fitted centre, finite nonzero
`abs(Omega0) <= 5e-5 s⁻¹`, and acceptable centre error e. The error ceiling is
`expm1(median(log1p(e)) + 4*MAD(log1p(e))/0.6745)`; zero MAD gives no upper ceiling.
An unavailable outer fit does not automatically reject an otherwise usable core
centre. See [the ESP method](https://doi.org/10.5194/egusphere-2026-3599).

### 3. Tracking and surface processing

Tracking searches for matches of the same cyclonicity within unit Euclidean
distance in `(x/50, y/50, Omega0/1e-5)`. It searches the closest preceding calendar
day first, up to four days back, then orders candidates by distance. An earlier
track ID can be assigned only once on the current day. Unmatched detections start
new tracks; this is a sequential assignment, not global trajectory optimisation.

Before gap filling, Q components receive a centred three-**observation** mean
with shorter end windows. Tracks must span at least 21 elapsed days. Aspect ratio
above 5 invalidates fitted parameters. Outer fits exceeding Rc=300 km,
Rc/R=1.75, abs(Omega)=5e-5 s⁻¹ or abs(psi0)=300 stored units instead fall back to
Omega0 and set Rc/psi0 missing. Missing calendar days are inserted. Parameters are
linearly interpolated, with endpoint filling and no configured maximum
outer-parameter gap. Tracks with an entirely missing required parameter or mean
R <=15 km are removed. IDs are renumbered, geographic coordinates/grid indices
assigned, and aspect ratio recalculated. Processed surface vorticity becomes
`w = Omega*(q11+q22)`.

These surface records can therefore contain smoothed/interpolated parameters.
There is no per-variable flag separating observations from interpolated values.
All fields in the released processed table are populated; that does not mean each
value was independently fitted on that day.

### 4. Vertical reconstruction and count-based QC

Three-dimensional velocity is transposed and reversed into surface-first order,
then linearly interpolated to the fixed reference depths without extrapolation.
Each processed surface solution seeds its profile. The preceding accepted centre
sets the next transect location. At each level, DOPPIO supplies the new core
centre/shape; continuity requires adequate transects, unchanged surface-relative
vorticity sign and jumps <=100 km. At least one directional radius is required.
The same outer-fit domain rule is used; both the rectangular and elliptical
selection domains must contain at least ten geometric points before ESP's own
finite-observation and rotation filtering.

If a continuity/fitting check returns no usable centre or raises an exception,
deeper following stops for that day. Missing outer parameters returned without
an exception can coexist with a recorded core centre. Confirmed-profile QC retains
eddy-days with at least **three recorded levels**, including surface. It is not a
complete finite-parameter, fit-quality or tilt-eligibility screen. The confirmed
release has 126,681 profiles; 6,357 processed days lack a confirmed profile.
62,054 confirmed depth records have missing Omega/Rc. Use the parameters required
by your analysis to define additional eligibility, rather than replacing missing
fits with zeros.

### 5. Delta tilt and pipeline summaries

The Delta estimator described in the accompanying manuscript (DOI: **XXXX**,
pending assignment) uses confirmed-profile centres within 0–1,000 m, leaving a deepest
available fitted level of 858.92 m. Centres are interpolated onto 10 m levels
within each day's retained range. Adjacent-depth increments are averaged over a
centred five-calendar-day window with Gaussian sigma 1 day; available-day weights
are renormalised separately at each interval. Only intervals supported by the
reference day's finite x/y centres enter averaging, variance estimation and
cumulative reconstruction. Neighbouring days cannot extend its depth range.

Depth weights are inverse sums of unweighted sample variances of x/y increments,
with regularisation `epsilon=1e-10 km²`; at least two contributing observations
are needed for a variance. Weighted SVD fits the reconstructed centreline using
the implementation's horizontal-km/depth-m numerical coordinates. It requires
at least five points spanning at least 200 m in upper-interval depth labels. For
example, a profile ending at 200 m supplies labels 0–190 m and fails that span.

`TiltDis` is fitted horizontal endpoint separation in km. `TiltDir` is the
**deep-to-shallow** bearing clockwise from geographic north, including the
configured +20-degree model-bearing offset. This differs from a raw surface-to-depth
centre offset and from maximum profile displacement. Evaluation windows must fit
within each eddy's confirmed-profile time bounds, excluding the first/last two
calendar days of that range. The YAML's `min_track_length_days` entry is not read
by this tilt stage; it is not an additional applied filter.

The tilt table contains all 120,348 temporal-window evaluation rows: 106,865 have
finite estimates and 13,483 have missing estimates. Missing reference profiles,
insufficient depth/temporal support or fitting failures can yield missing values;
these do **not** indicate physically untilted eddies. The other 12,690 processed
days have no evaluation row. Preserve the distinction by inspecting keys before
a left join, which otherwise represents both cases as missing values.

`analysis/tilt_summary.parquet` summarises finite Delta estimates.
`analysis/profile_tilt_summaries.parquet` instead uses each day's raw centre offsets
from its shallowest profile level: `max_tilt_distance` is the largest radial offset,
`max_depth` the deepest contributing level, and `tilt_slope` a linear distance/depth
slope. These are not alternative names for Delta tilt.

## Software, configuration and provenance

- [Ocean Eddy Pipeline v0.1.1](https://doi.org/10.5281/zenodo.24457378), archived
  commit `143b9a31ec0586206cf440df3da129fb886bc4b8`.
- The recorded generation run imported [ESP v1.1](https://doi.org/10.5281/zenodo.24182671).
  ESP v1.2 has a byte-identical computational `functions.py`; this does not change
  the historical generation version.
- Source model resources: [SEA-POD v3](https://doi.org/10.5281/zenodo.12703063)
  and [SEA-COFS](https://doi.org/10.5281/zenodo.8294716). The static grid is copied
  from reference file `outer_avg_01461.nc`, transposed consistently with extraction.
- `generation_config.yaml` is byte-identical to the configuration in the pipeline
  archive. Its absolute paths document the original installation; change them
  before running elsewhere. `parallel.skip_existing=false` allows replacement
  of outputs and should not be used on an irreplaceable downloaded copy.
- The recorded PBS run reported exit status 0 on 17 November 2026 and lists outputs
  through tilt analysis. It also recorded memory-limit events; the release audit
  separately verified expected partition coverage and primary-table consistency.
  The outputs do not embed a full historical environment/configuration manifest.
- `static_grid.nc` is a release-time extraction of reference geometry, not a new
  fitting product. Its provenance is recorded in `release_manifest.json`.
  Checksums identify the exact files selected for this release; they do not
  reconstruct undocumented historical changes.

## Access, citation and limitations

Parquet tables can be read with pandas/pyarrow; `z_r.npy` with NumPy;
`static_grid.nc` with netCDF4 or xarray. The geometry file contains no model
velocities or manuscript-specific analysis products. Model velocity fields are
needed to rerun extraction, but not to read the released eddy tables.

Cite this dataset as:

> Dowse, R. W. SEACOFS mesoscale eddy dataset (Version 1.0.0) [Data set].
> Zenodo. https://doi.org/10.5281/zenodo.24486070

Also cite the pipeline and ESP generation versions above when describing
construction. Reg W. Dowse is the sole creator
of this derived eddy dataset. Anyone may use, adapt and redistribute it, including
commercially, under CC BY 4.0 with appropriate credit, a licence link and an
indication of changes. See `LICENSE.txt` for source-model attribution; the static
grid and supplied depth array originate from the CC BY 4.0 model resource by
Junde Li, Colette Kerry and Moninya Roughan. Software licences are separate.
Replace the Delta manuscript DOI
placeholder `XXXX` and add its full citation when available. Tables give
consistently processed model fits, not uncertainty estimates
for individual centres; interpretation remains sensitive to model resolution,
identification, gap filling, profile depth support and the applied screening rules.
