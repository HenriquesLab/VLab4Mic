from vlab4mic import sweep_generator


def test_analysis_dict_is_per_instance():
    """Regression: `analysis` used to be a class attribute, so every
    sweep_generator instance shared one dict and results bled across runs.
    It must now be initialised per-instance in __init__."""
    a = sweep_generator.sweep_generator()
    b = sweep_generator.sweep_generator()

    assert a.analysis is not b.analysis
    assert a.analysis["unsorted"] is not b.analysis["unsorted"]

    a.analysis["unsorted"]["only_in_a"] = 123
    a.analysis["dataframes"] = "df_a"

    assert "only_in_a" not in b.analysis["unsorted"]
    assert b.analysis["dataframes"] is None
