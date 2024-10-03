from supramolsim.generate import labels
import os


def test_createlabel_from_config(configuration_directory):
    labelcongf_filename = "7R5K_Nup96_Cterminal_direct.yaml"
    local_path_congfig = os.path.join(
        configuration_directory, "labels", labelcongf_filename
    )
    fluorophoreID = "AF647"
    labelling_efficiency = 0.9
    label_object, label_params = labels.construct_label(
        local_path_congfig, fluorophoreID, labelling_efficiency
    )
    assert label_params["id"] == labelcongf_filename.split(".")[0]
