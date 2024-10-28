from tqdm import tqdm
from supramolsim.experiments import ExperimentParametrisation
import itertools
import numpy as np
from IPython.utils import io
from ..utils.transform.datatype import truncate
import matplotlib.pyplot as plt
import os
import pandas as pd
from .metrics import img_compare


def parameter_sweep_reps(
    Experiment: ExperimentParametrisation,
    sweep_parameters,
    write=False,
    repetitions=1,
    **kwargs,
):
    sweep_outputs = []
    # setup sweep parameters
    for parametername, pars in sweep_parameters.items():
        Experiment.sweep_pars[parametername] = pars
    Experiment._param_linspaces()
    reference = Experiment.gen_reference()
    out_dir = Experiment.output_directory
    # prepare combination of parameters
    linspaces_dict = Experiment.sweep_linspaces
    linspace_list = [linspaces_dict["labelling_efficiency"], linspaces_dict["defects"]]
    total_lengths = [len(i) for i in linspace_list]
    total_len = np.prod(np.array(total_lengths))
    # iterate over parameter combination
    print("Parameter sweep")
    iteration_params = []
    for combination in tqdm(itertools.product(*linspace_list), total=total_len):
        with io.capture_output() as captured:
            labeff = combination[0]
            defect = combination[1]
            iteration_output_reps = []
            for rep in range(repetitions):
                iteration_name = (
                    Experiment.experiment_id
                    + "LEff_"
                    + str(truncate(labeff, 3))
                    + "_Defect_"
                    + str(truncate(defect, 3))
                )
                # start of replicate
                _particle = Experiment._build_particle(
                    lab_eff=labeff, defect=defect, keep=True
                )
                if write:
                    _particle.show_instance(with_sources=True)
                    fig_name = iteration_name + "rep" + str(rep) + ".png"
                    name_path = os.path.join(out_dir, fig_name)
                    plt.savefig(name_path)
                    plt.close()
                exported_field = Experiment._build_coordinate_field(
                    use_self_particle=True
                )
                Experiment.imager.import_field(**exported_field)
                iteration_output = Experiment.run_simulation(
                    name=iteration_name, save=write
                )
                iteration_output_reps.append(iteration_output)
            sweep_outputs.append(iteration_output_reps)
            # end of replicate
            iteration_params.append(dict(labelling_effiency=labeff, defect=defect))
    return sweep_outputs, iteration_params, reference


def _reformat_img_stack(img, subregion=False, **kwargs):
    single_img = img[0]
    if subregion:
        single_img = single_img[
            subregion[0] : subregion[1], subregion[0] : subregion[1]
        ]
    return single_img
