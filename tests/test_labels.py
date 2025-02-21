from supramolsim.generate import labels
from supramolsim.workflows import probe_model
from supramolsim.utils.io.yaml_functions import load_yaml
import os


def test_createlabel_from_config(configuration_directory):
    labelcongf_filename = "NPC_Nup96_Cterminal_direct.yaml"
    local_path_congfig = os.path.join(
        configuration_directory, "probes", labelcongf_filename
    )
    fluorophoreID = "AF647"
    labelling_efficiency = 0.9
    label_object, label_params = labels.construct_label(
        local_path_congfig, fluorophoreID, labelling_efficiency
    )


def test_probe_model(configuration_directory):
    label_file = "Mock_primary_antibody" + ".yaml"
    local_path_congfig = os.path.join(configuration_directory, "probes", label_file)
    probe_params = load_yaml(local_path_congfig)
    probe, probe_emitters, anchor, ab_ref, probe_epitope = probe_model(
        **probe_params, config_dir=configuration_directory
    )
    assert anchor.shape == (3,)
    assert ab_ref.shape == (3,)
