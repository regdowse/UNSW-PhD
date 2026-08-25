# SEACOFS Eddy Dataset Modular Pipeline

This folder is a fresh modular version of the SEACOFS eddy dataset workflow.

The goal is to replace the current notebook-driven workflow with a reproducible, parallel pipeline that can run each stage in order, cache intermediate outputs, and resume safely after failures.

## Pipeline Stages

1. `detect_nencioli`
   - Reads `outer_avg_*.nc` model files.
   - Detects surface eddy candidates with the Nencioli method.
   - Writes partitioned detection tables by `fnumber`.

2. `fit_doppio_surface`
   - Reads Nencioli detections and the matching model file.
   - Fits DOPPIO surface eddy geometry and parameters.
   - Writes partitioned surface eddy tables by `fnumber`.

3. `track_eddies`
   - Reads all surface eddy partitions.
   - Applies QC and links eddies through time.
   - Writes the tracked eddy dataset.

4. `process_tracked_dataset`
   - Cleans and normalises the tracked dataset for downstream work.
   - Writes the processed eddy dataset.

5. `compute_vertical_profiles`
   - Reads the processed eddy dataset and 3D velocity fields.
   - Computes vertical profiles for each eddy-day.
   - Writes partitioned vertical profile outputs.

6. `qc_vertical_profiles`
   - Checks and confirms usable vertical profiles.
   - Writes the confirmed vertical profile dataset.

7. `compute_tilt`
   - Reads confirmed vertical profiles.
   - Computes tilt metrics and supporting diagnostics.
   - Writes the tilt dataset.

8. `analyse_tilt`
   - Reads the tilt dataset.
   - Produces figures, tables, and summary statistics.

## Design Principles

- Keep science kernels in importable Python modules, not notebooks.
- Keep notebooks for exploration, QC, and figures only.
- Write one output partition per independent work unit, usually `fnumber`, `day`, or `eddy_id`.
- Avoid repeatedly rewriting one large pickle during long jobs.
- Make every stage resumable by checking whether expected partition outputs already exist.
- Store paths, thresholds, and worker counts in config files.

## First Extraction Targets

The current notebooks map onto this folder roughly as follows:

- `nencioli_dataset.ipynb` -> `src/seacofs_eddy_dataset/stages/detection.py`
- `doppio_dataset.ipynb` -> `src/seacofs_eddy_dataset/stages/surface_fit.py`
- `eddy_tracking.ipynb` -> `src/seacofs_eddy_dataset/stages/tracking.py`
- `dataset_processing.ipynb` -> `src/seacofs_eddy_dataset/stages/processing.py`
- `vert_dataset.ipynb` -> `src/seacofs_eddy_dataset/stages/vertical_profiles.py`
- `tilt_dataset.ipynb` -> `src/seacofs_eddy_dataset/stages/tilt.py`

## Function Modules

Reusable functions live under `src/seacofs_eddy_dataset/core/`:

- `grid.py`: reference grid loading, Cartesian grid construction, filename helpers.
- `velocity.py`: fill-value cleaning, velocity rotation, interpolation.
- `nencioli.py`: Nencioli eddy detection kernel.
- `doppio.py`: transect, radius, local eddy geometry, and DOPPIO table helpers.
- `tracking.py`: surface eddy QC before temporal linking.
- `vertical.py`: 3D-to-depth interpolation and profile table helpers.
- `esp.py`: adapter for the external ESP_zonodo DOPPIO functions.
- `tilt.py`: tilt displacement, weighted tilt fitting, and summary metrics.

Stage modules under `src/seacofs_eddy_dataset/stages/` should orchestrate I/O,
parallel execution, and resume behaviour. They should call `core/` functions
rather than accumulating science logic directly.

## Example Usage

Install the package into the Python environment you will use for the run:

```bash
cd MRes/seacofs_eddy_dataset_modular
python -m pip install -e .
```

```bash
python -m seacofs_eddy_dataset.cli run-stage detect_nencioli --config config/example.yaml
python -m seacofs_eddy_dataset.cli run-stage fit_doppio_surface --config config/example.yaml
python -m seacofs_eddy_dataset.cli run-all --config config/example.yaml
```

## Output And Parallelisation

Set output locations and parallel worker counts in your config file. Start by
copying `config/example.yaml` to `config/local.yaml`, then edit the paths:

```yaml
paths:
  model_root: /srv/scratch/z3533156/26year_BRAN2020
  output_root: /srv/scratch/z5297792/SEACOFS_26yr_eddy_dataset_modular
  z_r: /srv/scratch/z5297792/z_r.npy
  esp_zonodo: /home/z5297792/ESP_zonodo

parallel:
  workers: 16
  backend: process
  skip_existing: true
```

`output_root` is the parent directory for every generated file. The pipeline
creates this layout underneath it:

```text
output_root/
  detections/
    fnumber=01461.parquet
  surface_eddies/
    fnumber=01461.parquet
  tracked/
    eddy_tracks.parquet
  processed/
    eddy_dataset_processed.parquet
  vertical_profiles/
    fnumber=01461.parquet
  vertical_profiles_confirmed/
    profiles.parquet
  tilt/
    tilt_dataset.parquet
  analysis/
    tilt_summary.parquet
    profile_tilt_summaries.parquet
```

`workers` controls how many model files are processed at once for the long
file-level stages: detection, DOPPIO surface fitting, and vertical profile
extraction. With `skip_existing: true`, reruns skip output partitions that
already exist, so interrupted runs can resume without starting from scratch.

There is also a concise notebook control panel at
`notebooks/run_pipeline.ipynb`. It loads `config/local.yaml`, runs selected
pipeline stages, and shows quick output summaries.

All stages are now implemented as runnable Python modules. The long-running
file-level stages (`detect_nencioli`, `fit_doppio_surface`, and
`compute_vertical_profiles`) write one parquet partition per `outer_avg_*.nc`
file using names like `fnumber=01461.parquet`. Tracking, processing, QC, tilt,
and analysis stages write consolidated parquet tables under `output_root`.

The DOPPIO surface and vertical stages still depend on the external
ESP_zonodo `functions.py` used by the original notebooks. Set
`paths.esp_zonodo` in your local config to that directory before running those
stages.
