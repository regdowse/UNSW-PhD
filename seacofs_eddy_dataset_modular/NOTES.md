# Modularisation Notes

## Immediate Plan

1. Extract shared grid setup from the notebooks into a reusable module.
2. Move `nencioli_dataset.ipynb` computational code into `stages/detection.py`.
3. Make detection write one parquet file per `fnumber`.
4. Move `doppio_dataset.ipynb` computational code into `stages/surface_fit.py`.
5. Make surface fitting read detection partitions and write one parquet file per `fnumber`.
6. Move tracking QC and linking into `stages/tracking.py`.
7. Move vertical profile extraction into `stages/vertical_profiles.py`.
8. Replace manual `part0` to `part9` notebooks with a parameterised parallel runner.

## Output Layout

The configured `output_root` should use a structure like:

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

## Current Notebook Mapping

- Detection: `MRes/SEACOFS_dataset/detect_track_tilt_dataset/nencioli_dataset.ipynb`
- Surface fitting: `MRes/SEACOFS_dataset/detect_track_tilt_dataset/doppio_dataset.ipynb`
- Tracking: `MRes/SEACOFS_dataset/detect_track_tilt_dataset/eddy_tracking.ipynb`
- Processing: `MRes/SEACOFS_dataset/detect_track_tilt_dataset/dataset_processing.ipynb`
- Vertical profiles: `MRes/SEACOFS_dataset/detect_track_tilt_dataset/vert_dataset.ipynb`
- Tilt dataset: `MRes/SEACOFS_dataset/detect_track_tilt_dataset/tilt_dataset.ipynb`
