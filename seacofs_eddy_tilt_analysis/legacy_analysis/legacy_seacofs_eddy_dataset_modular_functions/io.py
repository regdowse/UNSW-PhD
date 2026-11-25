from __future__ import annotations

from pathlib import Path

import pandas as pd


def ensure_dir(path: str | Path) -> Path:
    directory = Path(path)
    directory.mkdir(parents=True, exist_ok=True)
    return directory


def write_partition(df: pd.DataFrame, path: str | Path) -> Path:
    out_path = Path(path)
    ensure_dir(out_path.parent)
    tmp_path = out_path.with_suffix(out_path.suffix + ".tmp")
    df.to_parquet(tmp_path, index=False)
    tmp_path.replace(out_path)
    return out_path
