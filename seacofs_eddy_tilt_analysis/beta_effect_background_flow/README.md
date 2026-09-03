# Background-relative beta-effect analysis

This folder separates the expensive model-data processing from the statistical notebook.

Files:

- `01_build_background_cache.ipynb`: brief Jupyter launcher with an editable worker count.
- `background_flow_tools.py`: one-pass, restartable, file-parallel background-flow cache builder.
- `build_background_cache.py`: command-line entry point for a Katana batch job.
- `background_relative_beta_effect.ipynb`: residual propagation, vorticity-budget, and tilt tests.
- `vorticity_loss_during_beta_drift.ipynb`: off-shelf beta-aligned drift, daily weakening, cumulative `Delta w = -Delta f`, and whole-eddy bootstrap tests.
- `planetary_eddy_vertical_flow_shear.ipynb`: tests whether environmental vertical shear and directly tracked depth-specific centre propagation are consistent with the observed tilt of planetary-dominated eddies.
- `vertical_propagation_tools.py`: interpolates the vertical-profile centre positions to fixed physical depths and calculates geographic propagation velocities.

The annulus method has been removed. The primary background estimate is now a centred 91-day moving seasonal climatology at the eddy's grid cell. It is constructed from the 26-year monthly climatological fields by weighting each month according to the number of days it contributes to the 45-day-before/45-day-after window. This removes fixed-month discontinuities while retaining a tractable cache; it is a monthly-resolution approximation to a daily moving climatology. The full-archive mean remains a sensitivity test. Both estimates contain surface, thickness-weighted 0–200 m, and thickness-weighted 0–500 m velocities.

`z_r` is treated as a three-dimensional field of physical depths at sigma-level centres. It is not treated as a fixed z-level vector. Layer weights are reconstructed locally from adjacent sigma-centre depths, the surface, and bathymetry, then clipped at 500 m. The level-read optimization finds the deepest surface-counted sigma index needed anywhere in the sampled annuli. That number can be all sigma levels where sampled water columns are shallower than 500 m; it never assigns one fixed physical depth to a sigma index.

Run the cache builder on a compute node before executing the notebook:

```bash
cd ~/UNSW-PhD/seacofs_eddy_tilt_analysis/beta_effect_background_flow
python build_background_cache.py --workers 4
```

Start with four workers on the shared scratch filesystem. Increasing the worker count can reduce performance if NetCDF reads become I/O-bound. Completed file partitions are skipped on rerun.

Each archive file is opened once and contributes sums and counts to the monthly and full-archive climatological fields. Completed partitions are reused after interruption.

The revised workflow uses the `background_flow_cache_all_eddies_v4` directory. Run `01_build_background_cache.ipynb` once to build it; older annulus-cache partitions are not mixed into the new result.
