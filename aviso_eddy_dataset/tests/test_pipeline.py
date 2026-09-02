from aviso_eddy_dataset.pipeline import STAGES


def test_pipeline_has_surface_only_stages_in_order():
    assert [stage.name for stage in STAGES] == [
        "detect_nencioli",
        "fit_doppio_surface",
        "track_eddies",
        "process_tracked_dataset",
    ]
