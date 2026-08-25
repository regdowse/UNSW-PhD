from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml


@dataclass(frozen=True)
class PipelineConfig:
    raw: dict[str, Any]
    config_path: Path

    @property
    def model_root(self) -> Path:
        return Path(self.raw["paths"]["model_root"])

    @property
    def output_root(self) -> Path:
        return Path(self.raw["paths"]["output_root"])

    @property
    def workers(self) -> int:
        return int(self.raw.get("parallel", {}).get("workers", 1))

    @property
    def skip_existing(self) -> bool:
        return bool(self.raw.get("parallel", {}).get("skip_existing", True))


def load_config(path: str | Path) -> PipelineConfig:
    config_path = Path(path).expanduser().resolve()
    with config_path.open("r", encoding="utf-8") as handle:
        raw = yaml.safe_load(handle)
    if not isinstance(raw, dict):
        raise ValueError(f"Config must be a YAML mapping: {config_path}")
    return PipelineConfig(raw=raw, config_path=config_path)
