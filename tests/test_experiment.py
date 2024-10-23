from supramolsim.workflows import param_sweep_generator
from supramolsim.experiments import create_experiment_parametrisation


def test_experiment_parameterised():
    selected_mods = dict(
        STED_demo=None,
        Confocal_demo=None,
    )
    structure_and_labels = dict(
        structure_id="7R5K",
        structure_label="7R5K_Nup96_Cterminal_direct",
        fluorophore_id="AF647",
    )
    defects_eps_d = dict(eps1=300, eps2=600)
    sweep_pars = dict(
        labelling_efficiency=dict(start=0.2, end=1, nintervals=3),
        defects=dict(start=0, end=0.5, nintervals=3),
    )
    savging = dict(experiment_id="SupraMolSim_experiment", output_directory="")
    experiment_generator = create_experiment_parametrisation(
        structure_and_labels, selected_mods, defects_eps_d, sweep_pars, savging
    )
    param_sweep_generator(
        linspaces_dict=experiment_generator.sweep_linspaces,
        structure=experiment_generator.structure,
        imager=experiment_generator.imager,
        structure_label=experiment_generator.structure_label,
        fluorophore_id=experiment_generator.fluorophore_id,
        configuration_path=experiment_generator.configuration_path,
        defects_eps=experiment_generator.defect_eps,
        exp_name=experiment_generator.experiment_id,
        output_dir=experiment_generator.output_directory,
        write=False,
    )
