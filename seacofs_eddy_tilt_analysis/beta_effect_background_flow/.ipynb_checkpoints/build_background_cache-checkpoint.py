"""Katana command-line entry point for the background-flow cache."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parentelected = eddies.copy()
ANALYSIS_ROOT = HERE.parent
if str(ANALYSIS_ROOT) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_ROOT))

import seacofs_tilt_tools as tilt
from background_flow_tools import BackgroundConfig, build_background_cache


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--depth-max", type=float, default=500.0)
    parser.add_argument("--annulus-inner", type=float, default=1.5)
    parser.add_argument("--annulus-outer", type=float, default=3.0)
    args = parser.parse_args()

    paths = tilt.Paths()
    grid = tilt.load_grid(paths.grid, paths.z_r)
    eddies, _ = tilt.load_tilt_tables(paths)
    eddies = tilt.add_pv_gradient_terms(eddies, grid)
    # eddies = eddies[eddies["topo_plan_ratio"] < 0].copy()

    config = BackgroundConfig(
        depth_max_m=args.depth_max,
        annulus_inner_rc=args.annulus_inner,
        annulus_outer_rc=args.annulus_outer,
    )
    result = build_background_cache(eddies, grid, config=config, workers=args.workers)
    print(config.background_table_path)
    print(f"Rows: {len(result):,}; eddies: {result['Eddy'].nunique():,}")


if __name__ == "__main__":
    main()


