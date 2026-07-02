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


def test_set_sweep_parameters_reports_applied(capsys):
    g = sweep_generator.sweep_generator()
    capsys.readouterr()
    g.set_sweep_parameters(probe_DoL=[0.3, 0.6], pixelsize_nm=[80, 100])
    out = capsys.readouterr().out
    assert "Parameters set for sweep:" in out
    assert "probe_DoL" in out and "pixelsize_nm" in out
    assert "WARNING" not in out


def test_set_sweep_parameters_warns_on_ignored(capsys):
    g = sweep_generator.sweep_generator()
    capsys.readouterr()
    g.set_sweep_parameters(
        exp_time=[10, 20],
        peptide_motif={"motif": 1},
        minimal_distance=5,
        bogus_param=7,
    )
    out = capsys.readouterr().out
    assert "Parameters set for sweep: exp_time" in out
    assert "WARNING" in out
    for ignored in ("peptide_motif", "minimal_distance", "bogus_param"):
        assert ignored in out
    assert "exp_time" not in out.split("WARNING", 1)[1]


def test_set_sweep_parameters_silent_when_empty(capsys):
    g = sweep_generator.sweep_generator()
    capsys.readouterr()
    g.set_sweep_parameters()
    out = capsys.readouterr().out
    assert "Parameters set for sweep:" not in out
    assert "WARNING" not in out
