from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

from seacofs_eddy_dataset.config import PipelineConfig
from seacofs_eddy_dataset.stages import detection, processing, surface_fit, tilt, tracking, vertical_profiles

StageRunner = Callable[[PipelineConfig], None]


@dataclass(frozen=True)
class Stage:
    name: str
    runner: StageRunner
    description: str


STAGES: tuple[Stage, ...] = (
    Stage("detect_nencioli", detection.run, "Detect surface eddy candidates with Nencioli."),
    Stage("fit_doppio_surface", surface_fit.run, "Fit DOPPIO surface eddy parameters."),
    Stage("track_eddies", tracking.run, "Link surface eddies into time-continuous tracks."),
    Stage("process_tracked_dataset", processing.run, "Clean tracked eddy dataset."),
    Stage("compute_vertical_profiles", vertical_profiles.run, "Compute vertical profiles by eddy-day."),
    Stage("qc_vertical_profiles", vertical_profiles.run_qc, "Confirm usable vertical profiles."),
    Stage("compute_tilt", tilt.run, "Compute tilt metrics from confirmed vertical profiles."),
    Stage("analyse_tilt", tilt.run_analysis, "Produce tilt summaries, tables, and figures."),
)

STAGE_BY_NAME = {stage.name: stage for stage in STAGES}


def run_stage(name: str, config: PipelineConfig) -> None:
    try:
        stage = STAGE_BY_NAME[name]
    except KeyError as exc:
        valid = ", ".join(STAGE_BY_NAME)
        raise ValueError(f"Unknown stage {name!r}. Valid stages: {valid}") from exc
    stage.runner(config)


def run_all(config: PipelineConfig) -> None:
    for stage in STAGES:
        stage.runner(config)
