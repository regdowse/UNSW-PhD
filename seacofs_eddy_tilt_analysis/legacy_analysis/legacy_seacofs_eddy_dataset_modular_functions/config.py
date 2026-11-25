from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class PipelineConfig:
    raw: dict[str, Any]
    config_path: Path

    @property
    def output_root(self) -> Path:
        return Path(self.raw["paths"]["output_root"])

    @property
    def skip_existing(self) -> bool:
        return bool(self.raw.get("parallel", {}).get("skip_existing", True))
