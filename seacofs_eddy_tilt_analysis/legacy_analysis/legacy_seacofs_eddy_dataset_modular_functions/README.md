# Frozen SEACOFS helpers for legacy analysis

This local package replaces the dependency on the removed sibling
`seacofs_eddy_dataset_modular/`. The maintained pipeline lives at
https://github.com/regdowse/ocean-eddy-pipeline.

Source: `regdowse/UNSW-PhD` commit `c9a008f523572d6663083e9a769f7c25cfa4645d`, folder
`seacofs_eddy_dataset_modular/src/seacofs_eddy_dataset/`. Extraction date:
26 November 2026. This is a frozen reference subset, not a second pipeline.

Only the functions currently used by legacy analyses/tests, plus their required
helpers, are retained:

| Local module | Original module | Retained definitions |
| --- | --- | --- |
| `geometry.py` | `core/doppio.py` | `find_directional_radii`, `tangential_velocity` |
| `tilt.py` | `core/tilt.py` | `compute_weighted_tilt`, `bearing`, `summarise_tilt`, `add_tilt_displacement` |
| `tilt_stage.py` | `stages/tilt.py` | `run`, `run_analysis`, `_profiles_path` (used by the existing synthetic tests) |
| `io.py` | `io.py` | `write_partition`, `ensure_dir` |
| `config.py` | `config.py` | `PipelineConfig`, retaining only `output_root` and `skip_existing` properties |

Function bodies and their scientific conventions are unchanged. Imports now
refer to this package; no import of `seacofs_eddy_dataset`, old dataset source
path, or installed copy of the production pipeline is needed. The beta-flow
notebooks use their own `beta_effect_background_flow/analysis_tools.py`, so the
unrelated old dataset plotting helpers are not copied.

Add the `legacy_analysis/` directory to `sys.path`, then import from
`legacy_seacofs_eddy_dataset_modular_functions`. The affected notebook setup
cells and test suites already do this. Delta-method notebooks discover the
legacy tree from their current directory or its ancestors, including launches
from the repository or analysis root. Restart old notebook kernels to discard
previously imported pipeline modules.

Python 3.10+, NumPy and pandas are required; actual Parquet I/O also needs a
Parquet engine. This removes the old code-folder dependency, not the existing
external data, ESP backend, or analysis-specific scientific libraries. Scratch
paths, existing notebook outputs, and scientific results are unchanged. No
multi-year caches or datasets need regenerating for this code relocation.

The unused old dataset notebooks are not copied here. Their original versions
remain recoverable from the source commit above.

From the UNSW-PhD root, check the affected workflows with:

```bash
MPLBACKEND=Agg python -m unittest discover -s seacofs_eddy_tilt_analysis/legacy_analysis/delta_tilt_method -p 'test_*.py' -v
MPLBACKEND=Agg python -m unittest discover -s seacofs_eddy_tilt_analysis/legacy_analysis/esp_population_composites -p 'test_*.py' -v
```
