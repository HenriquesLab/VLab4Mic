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


def analyse_sweep_reps_vectors(
    img_outputs, img_params, reference, analysis_case_params, **kwargs
):
    conditions = list(reference.keys())
    references = dict()
    queries = dict()
    measurement_vectors = list()
    for case in conditions:
        item_vector = []
        queries[case] = []
        data_pivot = []
        measurement_per_combination = []
        sd_measurement_per_combination = []
        reference_img = _reformat_img_stack(
            reference[case], **analysis_case_params[case]
        )
        references[case] = reference_img
        # print(reference_img.shape)
        for i in range(len(img_params)):
            sweep_case_measurements = []
            case_iteration_replicas = []
            r = 0
            for simu_replica in img_outputs[i]:
                item_vector = []
                # print(simu_replica.keys(), i, case)
                case_iteration = _reformat_img_stack(
                    simu_replica[case], **analysis_case_params[case]
                )
                rep_measurement = img_compare(
                    reference_img, case_iteration, **analysis_case_params[case]
                )
                item_vector.append(case)
                item_vector.append(img_params[i]["labelling_effiency"])
                item_vector.append(img_params[i]["defect"])
                item_vector.append(rep_measurement)
                item_vector.append(r)
                # print(item_vector)
                measurement_vectors.append(item_vector)
                sweep_case_measurements.append(rep_measurement)
                case_iteration_replicas.append(case_iteration)
                r = r + 1
            queries[case].append(dict(im=case_iteration_replicas, params=img_params[i]))
            measurement_per_combination.append(np.mean(sweep_case_measurements))
            sd_measurement_per_combination.append(np.std(sweep_case_measurements))
        measurement_array = np.array(measurement_vectors)
        print(measurement_array)
        data_frame = pd.DataFrame(
            data={
                "Labelling efficiency": np.array(
                    measurement_array[:, 1], dtype=np.float32
                ),
                "Fracitonal defect": np.array(
                    measurement_array[:, 2], dtype=np.float32
                ),
                "Condition": measurement_array[:, 0],
                "Metric": np.array(measurement_array[:, 3], dtype=np.float32),
                "replica": measurement_array[:, 4],
            }
        )
    return data_frame, references, queries
