from __future__ import annotations

import numpy as np
import pandas as pd

from seacofs_eddy_dataset.config import PipelineConfig
from seacofs_eddy_dataset.core.tilt import compute_weighted_tilt, summarise_tilt
from seacofs_eddy_dataset.io import write_partition


def _profiles_path(config: PipelineConfig):
    return config.output_root / "vertical_profiles_confirmed" / "profiles.parquet"


def run(config: PipelineConfig) -> None:
    out_path = config.output_root / "tilt" / "tilt_dataset.parquet"
    if config.skip_existing and out_path.exists():
        print(out_path)
        return

    source = _profiles_path(config)
    if not source.exists():
        raise FileNotFoundError(f"Missing confirmed vertical profiles: {source}")
    profiles = pd.read_parquet(source)
    if profiles.empty:
        write_partition(pd.DataFrame(columns=["Eddy", "Day", "TiltDis", "TiltDir"]), out_path)
        print(out_path)
        return

    settings = config.raw.get("tilt", {})
    frames = [
        compute_weighted_tilt(
            profiles,
            int(eddy),
            num=int(settings.get("smoothing_days", 6)),
            depth_int=int(settings.get("depth_interval_m", 10)),
            max_depth=int(settings.get("max_depth_m", 1000)),
            min_depth_range=int(settings.get("min_depth_range_m", 200)),
            min_points=int(settings.get("min_points", 5)),
            bearing_offset=float(settings.get("bearing_offset_deg", 20.0)),
        )
        for eddy in sorted(profiles["Eddy"].dropna().unique())
    ]
    tilt = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    write_partition(tilt, out_path)
    print(out_path)


def run_analysis(config: PipelineConfig) -> None:
    tilt_path = config.output_root / "tilt" / "tilt_dataset.parquet"
    profiles_path = _profiles_path(config)
    if not tilt_path.exists():
        raise FileNotFoundError(f"Missing tilt dataset: {tilt_path}")

    tilt = pd.read_parquet(tilt_path)
    profiles = pd.read_parquet(profiles_path) if profiles_path.exists() else pd.DataFrame()
    out_dir = config.output_root / "analysis"

    summary = {
        "n_tilt_rows": len(tilt),
        "n_valid_tilts": int(tilt["TiltDis"].notna().sum()) if "TiltDis" in tilt else 0,
        "median_tilt_distance_km": float(np.nanmedian(tilt["TiltDis"])) if len(tilt) else np.nan,
        "mean_tilt_distance_km": float(np.nanmean(tilt["TiltDis"])) if len(tilt) else np.nan,
    }
    write_partition(pd.DataFrame([summary]), out_dir / "tilt_summary.parquet")

    if not profiles.empty:
        rows = []
        for (eddy, day), profile in profiles.groupby(["Eddy", "Day"]):
            item = {"Eddy": eddy, "Day": day}
            item.update(summarise_tilt(profile))
            rows.append(item)
        write_partition(pd.DataFrame(rows), out_dir / "profile_tilt_summaries.parquet")

    print(out_dir)
