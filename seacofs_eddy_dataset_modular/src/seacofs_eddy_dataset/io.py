from __future__ import annotations

from pathlib import Path

import pandas as pd


def ensure_dir(path: str | Path) -> Path:
    directory = Path(path)
    directory.mkdir(parents=True, exist_ok=True)
    return directory


def partition_path(root: str | Path, stage: str, key: str, suffix: str = ".parquet") -> Path:
    return Path(root) / stage / f"{key}{suffix}"


def write_partition(df: pd.DataFrame, path: str | Path) -> Path:
    out_path = Path(path)
    ensure_dir(out_path.parent)
    tmp_path = out_path.with_suffix(out_path.suffix + ".tmp")
    df.to_parquet(tmp_path, index=False)
    tmp_path.replace(out_path)
    return out_path


def read_partitions(paths: list[str | Path]) -> pd.DataFrame:
    frames = [pd.read_parquet(path) for path in paths]
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)
