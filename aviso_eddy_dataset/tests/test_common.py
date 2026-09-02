from pathlib import Path

import pandas as pd

from aviso_eddy_dataset.common import day_number, file_key, partition_path
from aviso_eddy_dataset.config import PipelineConfig


def config(tmp_path):
    return PipelineConfig(
        raw={
            "paths": {
                "data_root": str(tmp_path),
                "output_root": str(tmp_path / "output"),
                "esp_zonodo": str(tmp_path),
            },
            "processing": {"date_origin": "1990-01-01"},
        },
        config_path=tmp_path / "config.yaml",
    )


def test_day_number_is_continuous_across_year_boundary(tmp_path):
    cfg = config(tmp_path)
    assert day_number("1994-12-31", cfg) + 1 == day_number("1995-01-01", cfg)


def test_partition_names_are_source_specific(tmp_path):
    cfg = config(tmp_path)
    source = Path("/data/AVISO_0.125_EAC_1994.nc")
    assert file_key(source) == "AVISO_0.125_EAC_1994"
    assert partition_path(cfg, "detections", source) == (
        tmp_path / "output/detections/source=AVISO_0.125_EAC_1994.parquet"
    )
