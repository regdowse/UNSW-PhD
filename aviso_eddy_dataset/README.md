# AVISO Eddy Dataset

This package applies the SEACOFS surface eddy framework to daily AVISO
altimetry. It uses absolute geostrophic velocities (`ugos`, `vgos`) and runs:

1. Nencioli candidate detection.
2. DOPPIO surface fitting.
3. Continuous temporal tracking across all annual files.
4. Quality control and final dataset processing.

There are no subsurface, vertical-profile, or tilt stages.

## Grid handling

AVISO longitude and latitude are converted to increasing Cartesian kilometre
axes. The unrotated velocities are treated internally in `(longitude,
latitude)` order.

Only Nencioli uses the configurable 1 km interpolated grid. Every detected
centre is mapped back to the nearest native AVISO grid cell. DOPPIO then uses
native-resolution `ugos`, `vgos`, and native Cartesian coordinates.

## Reused SEACOFS functions

The package imports the established science kernels from the sibling
`seacofs_eddy_dataset_modular` package: Nencioli detection, DOPPIO geometry
helpers, the ESP function loader, tracking, processing helpers, and atomic
Parquet I/O. All AVISO-specific loading and orchestration remains in this
folder; the SEACOFS source is not modified.

## Installation on Katana

From the `UNSW-PhD` repository root:

```bash
python -m pip install -e ./seacofs_eddy_dataset_modular
python -m pip install -e ./aviso_eddy_dataset
```

Copy and edit the example configuration if a machine-specific version is
needed:

```bash
cp aviso_eddy_dataset/config/example.yaml aviso_eddy_dataset/config/local.yaml
```

The committed example already points to:

```text
input:  /srv/scratch/z5502183/AVISO_0.125
output: /srv/scratch/z5297792/aviso_eddy_dataset
ESP:    /home/z5297792/ESP_zonodo
```

## Running

List stages:

```bash
aviso-eddy list-stages
```

Run one stage:

```bash
aviso-eddy run-stage detect_nencioli --config aviso_eddy_dataset/config/example.yaml
aviso-eddy run-stage fit_doppio_surface --config aviso_eddy_dataset/config/example.yaml
aviso-eddy run-stage track_eddies --config aviso_eddy_dataset/config/example.yaml
aviso-eddy run-stage process_tracked_dataset --config aviso_eddy_dataset/config/example.yaml
```

Run the complete pipeline:

```bash
aviso-eddy run-all --config aviso_eddy_dataset/config/example.yaml
```

Detection and surface fitting write one resumable Parquet partition per annual
source file. Tracking reads every partition together, sorts by the continuous
day coordinate, and therefore tracks across year boundaries.

## Outputs

```text
/srv/scratch/z5297792/aviso_eddy_dataset/
  detections/source=AVISO_0.125_EAC_1994.parquet
  surface_eddies/source=AVISO_0.125_EAC_1994.parquet
  tracked/eddy_tracks.parquet
  processed/eddy_dataset_processed.parquet
```

Set `parallel.skip_existing: true` to resume the expensive file-level stages.
Set it to `false` when parameters change and partitions must be regenerated.
