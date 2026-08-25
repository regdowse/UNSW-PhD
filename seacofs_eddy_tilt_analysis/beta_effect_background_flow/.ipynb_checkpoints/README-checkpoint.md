# Background-relative beta-effect analysis

This folder separates the expensive model-data processing from the statistical notebook.

Files:

- `01_build_background_cache.ipynb`: brief Jupyter launcher with an editable worker count.
- `background_flow_tools.py`: one-pass, restartable, file-parallel background-flow cache builder.
- `build_background_cache.py`: command-line entry point for a Katana batch job.
- `background_relative_beta_effect.ipynb`: residual propagation, vorticity-budget, and tilt tests.

The primary background estimate is the instantaneous surface velocity in an annulus around each eddy. The workflow also calculates a thickness-weighted 0–500 m annulus velocity and monthly climatological surface and 0–500 m velocities.

`z_r` is treated as a three-dimensional field of physical depths at sigma-level centres. It is not treated as a fixed z-level vector. Layer weights are reconstructed locally from adjacent sigma-centre depths, the surface, and bathymetry, then clipped at 500 m. The level-read optimization finds the deepest surface-counted sigma index needed anywhere in the sampled annuli. That number can be all sigma levels where sampled water columns are shallower than 500 m; it never assigns one fixed physical depth to a sigma index.

Run the cache builder on a compute node before executing the notebook:

```bash
cd ~/UNSW-MRes/MRes/seacofs_eddy_tilt_analysis/beta_effect_background_flow
python build_background_cache.py --workers 4
```

Start with four workers on the shared scratch filesystem. Increasing the worker count can reduce performance if NetCDF reads become I/O-bound. Completed file partitions are skipped on rerun.

Each archive file is opened once. During that pass it contributes to the full-archive monthly climatology, and—only where an eddy is present—to cropped annulus estimates. This replaces separate archive passes and repeated full-grid annulus masks.
