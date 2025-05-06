from .experiments import ExperimentParametrisation
from .analysis import sweep

class sweep_generator:
    experiment = ExperimentParametrisation()
    structures = ["1XI5", ]
    probes = ["NHS_ester", ]
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
    virtual_samples = None
    virtual_samples_parameters = None


    def generate_virtual_samples(self):
        self.experiment, self.virtual_samples, self.virtual_samples_parameters = sweep.sweep_vasmples(
            structures=self.structures,
            probes=self.probes, 
            probe_parameters=self.probe_parameters, 
            virtual_samples=self.vsample_parameters,
            repetitions=self.sweep_repetitions)
