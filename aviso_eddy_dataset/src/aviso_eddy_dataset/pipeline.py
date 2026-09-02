from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

from . import detection, processing, surface_fit, tracking
from .config import PipelineConfig


@dataclass(frozen=True)
class Stage:
    name: str
    description: str
    run: Callable[[PipelineConfig], None]


STAGES = [
    Stage("detect_nencioli", "Detect daily candidates from 1 km interpolated ugos/vgos.", detection.run),
    Stage("fit_doppio_surface", "Fit DOPPIO on native-resolution ugos/vgos.", surface_fit.run),
    Stage("track_eddies", "Track surface eddies continuously across all source files.", tracking.run),
    Stage("process_tracked_dataset", "Apply QC and create the processed AVISO dataset.", processing.run),
]
STAGE_BY_NAME = {stage.name: stage for stage in STAGES}


def run_stage(name: str, config: PipelineConfig) -> None:
    STAGE_BY_NAME[name].run(config)


def run_all(config: PipelineConfig) -> None:
    for stage in STAGES:
        print(f"==> {stage.name}: {stage.description}")
        stage.run(config)
