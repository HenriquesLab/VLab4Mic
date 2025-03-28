from supramolsim import experiments 

def test_experiment_parameterised():
    selected_mods = dict(
        STED=None,
        Confocal=None,
    )
    structure_and_labels = dict(
        structure_id="7R5K",
        structure_label="NPC_Nup96_Cterminal_direct",
        fluorophore_id="AF647",
    )
    defects_eps_d = dict(eps1=300, eps2=600)
    sweep_pars = dict(
        labelling_efficiency=dict(start=0.2, end=1, nintervals=3),
        defects=dict(start=0, end=0.5, nintervals=3),
    )
    savging = dict(experiment_id="SupraMolSim_experiment", output_directory="")
    experiment_generator = experiments.create_experiment_parametrisation(
        structure_and_labels=structure_and_labels,
        modalities_acquisition=selected_mods,
        defects_params=defects_eps_d,
        savging=savging,
        params2sweep=sweep_pars,
    )


def test_experiment_runsimulation():
    selected_mods = dict(
        STED=None,
        Confocal=None,
    )
    structure_and_labels = dict(
        structure_id="2RCJ",
        structure_label="NHS_ester",
        fluorophore_id="AF647",
    )
    savging = dict(experiment_id="SupraMolSim_experiment", output_directory="")
    experiment_generator = experiments.create_experiment_parametrisation(
        structure_and_labels=structure_and_labels,
        modalities_acquisition=selected_mods,
        savging=savging,
    )
    experiment_generator.run_simulation()


def test_gen_virtual_sample():
    vsample, testexperiment = experiments.generate_virtual_sample()
    print(vsample)
    assert len(vsample["reference_point"]) == 3


def test_image_sample():
    images, experiment = experiments.image_vsample()

def test_build_virtual_microscope():
    vmicroscope, experiment = experiments.build_virtual_microscope()
