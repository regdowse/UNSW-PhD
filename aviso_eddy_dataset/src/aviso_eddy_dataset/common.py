from __future__ import annotations

from pathlib import Path

import pandas as pd
import xarray as xr

from .config import PipelineConfig


def find_data_files(config: PipelineConfig) -> list[Path]:
    pattern = config.raw.get("files", {}).get("pattern", "AVISO_0.125_EAC_*.nc")
    files = sorted(config.data_root.glob(pattern))
    if not files:
        raise FileNotFoundError(f"No AVISO files matching {pattern!r} in {config.data_root}")
    return files


def file_key(path: Path) -> str:
    return path.stem.replace(" ", "_")


def partition_path(config: PipelineConfig, stage: str, path: Path) -> Path:
    return config.output_root / stage / f"source={file_key(path)}.parquet"


def open_aviso(path: Path, config: PipelineConfig) -> xr.Dataset:
    variables = config.raw.get("variables", {})
    required = {
        variables.get("longitude", "longitude"),
        variables.get("latitude", "latitude"),
        variables.get("time", "time"),
        variables.get("u", "ugos"),
        variables.get("v", "vgos"),
    }
    dataset = xr.open_dataset(path, decode_times=True)
    missing = required.difference(dataset.variables)
    if missing:
        dataset.close()
        raise KeyError(f"{path.name} is missing variables: {sorted(missing)}")
    return dataset


def day_number(value, config: PipelineConfig) -> int:
    origin = pd.Timestamp(config.raw.get("processing", {}).get("date_origin", "1990-01-01"))
    timestamp = pd.Timestamp(value)
    delta = timestamp.normalize() - origin.normalize()
    if delta.components.hours or delta.components.minutes or delta.components.seconds:
        raise ValueError(f"AVISO timestamp is not on an integral day relative to {origin}: {timestamp}")
    return int(delta.days)


def source_file_by_day(config: PipelineConfig) -> pd.Series:
    """Map every source time coordinate to its file, checking for overlaps."""

    time_name = config.raw.get("variables", {}).get("time", "time")
    records: list[tuple[int, str]] = []
    for path in find_data_files(config):
        with open_aviso(path, config) as dataset:
            records.extend((day_number(value, config), str(path)) for value in dataset[time_name].values)
    mapping = pd.DataFrame(records, columns=["Day", "source_file"])
    duplicates = mapping.loc[mapping.duplicated("Day", keep=False)]
    if not duplicates.empty:
        days = sorted(duplicates.Day.unique().tolist())
        raise ValueError(f"Overlapping AVISO source times for day numbers: {days[:10]}")
    return mapping.set_index("Day")["source_file"].sort_index()
