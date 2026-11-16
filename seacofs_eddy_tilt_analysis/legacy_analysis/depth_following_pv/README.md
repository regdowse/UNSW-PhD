# Depth-following PV-gradient cache

`build_depth_following_pv_cache.ipynb` follows every eddy's fitted vertical
centreline from the surface to 1000 m. At each fitted depth it averages the
completed shallow-water PV-gradient components inside that level's own ellipse.
It saves:

- a depth-resolved table with one row per eddy/day/depth; and
- a thickness-weighted snapshot table with one row per eddy/day.

Directions and magnitudes in the snapshot table are reconstructed after the
east/north vector components are averaged. Run this expensive processing
notebook once on the HPC system, then load its Parquet outputs in subsequent
analyses.

The corresponding API is:

```python
snapshot_df, depth_df = tilt.add_pv_gradient_terms(
    eddies,
    grid,
    core_mean=True,
    depth_following=True,
    vertical=tilt.load_vert(paths, dic_form=False),
    max_depth_m=1000,
)
```

Calls without `depth_following=True` remain unchanged and return one DataFrame.

After the cache has been built, all three analysis inputs can be selected from
the same public function:

```python
# Calculate the original surface-centred result.
original_df = tilt.add_pv_gradient_terms(
    eddies, grid, core_mean=True, source="original"
)

# Load the saved one-row-per-eddy/day depth-following result.
snapshot_df = tilt.add_pv_gradient_terms(source="depth_snapshot")

# Load the saved one-row-per-eddy/day/depth result.
depth_df = tilt.add_pv_gradient_terms(source="depth")
```

Pass `cache_root=...` to load an equivalent cache from another directory.
