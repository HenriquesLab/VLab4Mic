from .experiments import ExperimentParametrisation
from .analysis import sweep
import matplotlib.pyplot as plt


class sweep_generator:
    experiment = ExperimentParametrisation()
    structures = [
        "1XI5",
    ]
    probes = [
        "NHS_ester",
    ]
    modalities = ["Widefield", "Confocal", "STED", "SMLM"]
    sweep_repetitions = 3
    probe_parameters = None
    defect_parameters = None
    vsample_parameters = None
    modality_parameters = None
    acquisition_parameters = None
    # reference parameters
    reference_structure = "1XI5"
    reference_probe = "NHS_ester"
    reference_probe_parameters = {"labelling_efficiency": 1.0}
    # outputs
    reference_virtual_sample = None
    reference_virtual_sample_params = None
    reference_image = None
    virtual_samples = None
    virtual_samples_parameters = None
    acquisition_outputs = None
    acquisition_outputs_parameters = None

    def generate_virtual_samples(self):
        self.experiment, self.virtual_samples, self.virtual_samples_parameters = (
            sweep.sweep_vasmples(
                structures=self.structures,
                probes=self.probes,
                probe_parameters=self.probe_parameters,
                virtual_samples=self.vsample_parameters,
                repetitions=self.sweep_repetitions,
            )
        )

    def image_virtual_samples(self):
        if self.virtual_samples is None:
            self.generate_virtual_samples()
        (
            self.experiment,
            self.acquisition_outputs,
            self.acquisition_outputs_parameters,
            mod_pixelsizes,
        ) = sweep.sweep_modalities_updatemod(
            experiment=self.experiment,
            vsample_outputs=self.virtual_samples,
            vsampl_pars=self.virtual_samples_parameters,
            modalities=self.modalities,
            modality_acq_prams=self.acquisition_parameters,
        )

    def generate_reference_sample(self):
        self.reference_virtual_sample, self.reference_virtual_sample_params = (
            sweep.generate_global_reference_sample(
                structure=self.reference_structure,
                probe=self.reference_probe,
                probe_parameters=self.reference_probe_parameters,
            )
        )

    def generate_reference_image(self):
        if self.reference_image is None:
            self.generate_reference_sample()
        self.reference_image, self.reference_image_parameters = (
            sweep.generate_global_reference_modality(
                reference_vsample=self.reference_virtual_sample,
                reference_vsample_params=self.reference_virtual_sample_params,
            )
        )

    def preview_acquisition_output(self, return_image=False):
        param_id = list(self.acquisition_outputs_parameters.keys())[0]
        repetition = 0
        frame = 0
        if return_image:
            return self.acquisition_outputs[param_id][repetition][frame]
        else:
            plt.imshow(self.acquisition_outputs[param_id][repetition][frame])
            print(self.acquisition_outputs_parameters[param_id])

    def preview_reference_image(self, return_image=False):
        if return_image:
            return self.reference_image
        else:
            plt.imshow(self.reference_image[0])
            print(self.reference_image_parameters)
