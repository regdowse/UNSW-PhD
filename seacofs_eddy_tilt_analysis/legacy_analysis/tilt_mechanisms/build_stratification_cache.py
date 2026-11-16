"""Katana CLI for the restartable, file-parallel eddy-core N2 cache."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ANALYSIS_ROOT = HERE.parent
if str(ANALYSIS_ROOT) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_ROOT))
if str(ANALYSIS_ROOT / "beta_effect_background_flow") not in sys.path:
    sys.path.insert(0, str(ANALYSIS_ROOT / "beta_effect_background_flow"))

import seacofs_tilt_tools as tilt
from mechanism_tools import N2CacheConfig, build_n2_cache_xroms


def main():
    defaults = N2CacheConfig()
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--point-batch-size", type=int, default=128)
    parser.add_argument("--model-root", type=Path, default=defaults.model_root)
    parser.add_argument("--output", type=Path, default=defaults.output_path)
    parser.add_argument("--grid", type=Path, default=tilt.Paths().grid)
    parser.add_argument("--z-r", type=Path, default=tilt.Paths().z_r)
    parser.add_argument("--overwrite-partitions", action="store_true")
    args = parser.parse_args()

    eddies, _ = tilt.load_tilt_tables(tilt.Paths())
    result = build_n2_cache_xroms(
        eddies,
        model_root=args.model_root,
        output_path=args.output,
        grid_path=args.grid,
        z_r_path=args.z_r,
        workers=args.workers,
        point_batch_size=args.point_batch_size,
        skip_existing=not args.overwrite_partitions,
    )
    print(args.output)
    print(f"Rows: {len(result):,}; eddies: {result['Eddy'].nunique():,}")


if __name__ == "__main__":
    main()
