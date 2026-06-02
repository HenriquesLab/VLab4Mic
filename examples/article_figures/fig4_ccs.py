from vlab4mic import experiments
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import os
import tempfile
import time
from pathlib import Path
from urllib.parse import unquote, urlparse

import requests


DOWNLOAD_HEADERS = {
    "User-Agent": "VLab4Mic/0.1 fig4_ccs downloader",
    "Accept": "application/octet-stream,*/*",
}


def download_cif_file(
    url: str,
    timeout: tuple[int, int] = (10, 120),
    retries: int = 4,
    chunk_size: int = 1024 * 1024,
) -> str:
    """
    Download a CIF file from a URL and save it to a temporary directory.
    
    Parameters
    ----------
    url : str
        The URL of the CIF file to download.
    timeout : tuple[int, int], optional
        Connect and read timeout in seconds (default: (10, 120)).
    retries : int, optional
        Number of download attempts before failing (default: 4).
    chunk_size : int, optional
        Bytes read at a time while streaming the response (default: 1 MiB).
    
    Returns
    -------
    str
        Full path to the downloaded file.
    
    Raises
    ------
    Exception
        If the download fails (connection error, HTTP error, etc.).
    """
    cache_dir = Path(tempfile.gettempdir()) / "vlab4mic_fig4_ccs"
    cache_dir.mkdir(parents=True, exist_ok=True)

    filename = Path(unquote(urlparse(url).path)).name or "structure.cif"
    filepath = cache_dir / filename
    partial_path = filepath.with_suffix(filepath.suffix + ".part")

    expected_size = None
    try:
        response = requests.head(
            url,
            headers=DOWNLOAD_HEADERS,
            allow_redirects=True,
            timeout=timeout,
        )
        response.raise_for_status()
        expected_size = int(response.headers.get("content-length", 0)) or None
    except requests.RequestException as exc:
        print(f"Could not check remote file size for {url}: {exc}")

    if filepath.exists() and filepath.stat().st_size > 0:
        if expected_size is None or filepath.stat().st_size == expected_size:
            print(f"Using cached CIF file: {filepath}")
            return str(filepath)
        print(
            f"Cached file is incomplete ({filepath.stat().st_size} of "
            f"{expected_size} bytes). Downloading again..."
        )

    last_error = None
    for attempt in range(1, retries + 1):
        try:
            print(f"Downloading {url} (attempt {attempt}/{retries})...")
            with requests.get(
                url,
                headers=DOWNLOAD_HEADERS,
                stream=True,
                timeout=timeout,
            ) as response:
                if response.status_code == 429:
                    retry_after = response.headers.get("Retry-After")
                    raise RuntimeError(
                        "Zenodo rate limit reached"
                        + (f"; retry after {retry_after} seconds" if retry_after else "")
                    )
                if response.status_code == 403:
                    raise RuntimeError("Zenodo returned HTTP 403; this may be an IP block")

                response.raise_for_status()
                total_size = int(response.headers.get("content-length", 0)) or expected_size
                bytes_written = 0

                with open(partial_path, "wb") as out_file:
                    for chunk in response.iter_content(chunk_size=chunk_size):
                        if not chunk:
                            continue
                        out_file.write(chunk)
                        bytes_written += len(chunk)

                if total_size is not None and bytes_written != total_size:
                    raise RuntimeError(
                        f"incomplete download: got {bytes_written} of {total_size} bytes"
                    )

                os.replace(partial_path, filepath)
                print(f"Downloaded to: {filepath}")
                return str(filepath)

        except (requests.RequestException, RuntimeError) as exc:
            last_error = exc
            if partial_path.exists():
                partial_path.unlink()
            if attempt < retries:
                sleep_seconds = 2 ** (attempt - 1)
                print(f"Download failed: {exc}. Retrying in {sleep_seconds}s...")
                time.sleep(sleep_seconds)

    raise RuntimeError(f"Failed to download CIF file from {url}: {last_error}")


########################## 
random_seed = 24
dome_url = "https://zenodo.org/records/20377070/files/Dome_model.cif"
flat_url = "https://zenodo.org/records/20377070/files/Flat_lattice_model.cif"

# Download the CIF files
dome_model = download_cif_file(dome_url)
flat_model = download_cif_file(flat_url)
primary = dict(
    probe_template = "Antibody",
    probe_name="custom",
    probe_target_type = "Sequence",
    probe_target_value = "NNRIA",
    probe_distante_to_epitope = 0,
    probe_DoL=4,  
)
############################ DOME #######################################
images, noiseless, my_experiment = experiments.image_vsample(
    structure=dome_model,
    structure_is_path=True,
    probe_list=[primary,],
    clear_experiment=True,
    multimodal=["Widefield", "Confocal", "AiryScan", "STED", "SMLM"],
    run_simulation=False,
    sample_dimensions=[500,500,10],
    random_seed=random_seed
)
my_experiment.update_modality(modality_name="STED", depth_of_field_nm=1000)
my_experiment.update_modality(modality_name="SMLM", depth_of_field_nm=1000)
images, noiseless = my_experiment.run_simulation()
fig = plt.figure(figsize=[10,20])
ax = fig.add_subplot(1, 3, 1, projection="3d")
total = my_experiment.structure.num_assembly_atoms
fraction = 50000/total
my_experiment.structure.show_target_labels(
    with_assembly_atoms=True, 
    axis_object=ax, 
    show_axis=False, 
    assembly_fraction = fraction,
    view_init=[90,0,0],
    target_size = 0,
    atoms_size = 1,
    atoms_alpha = 0.01,
    reference_point = False,
    axesoff=True,
    )
ax = fig.add_subplot(1, 3, 2, projection="3d")
my_experiment.structure.show_target_labels(
    with_assembly_atoms=False, 
    axis_object=ax, 
    show_axis=False, 
    view_init=[90,0,0],
    target_size = 10,
    reference_point = False,
    )
ax = fig.add_subplot(1, 3, 3, projection="3d")
my_experiment.particle.gen_axis_plot(
 axis_object=ax,
 view_init=[90,0,0],
 emitter_plotsize=5,
 xlim=[-1250,1250], ylim=[-1250,1250], zlim=[0,1000]
)

filename = my_experiment.date_as_string + 'vlab4mic_ccs_DOME_labelled_structure.png'
filename2 = os.path.join(my_experiment.output_directory, filename)
fig.savefig(filename2,dpi=300, bbox_inches='tight')
plt.close()

fig = plt.figure(figsize=[20,30])
for i, mod_name in enumerate(images.keys()):
    ax = fig.add_subplot(1,5,i+1)
    ax.imshow(images[mod_name]["ch0"][0], cmap="grey")
    ax.set_title(mod_name)
    ax.set_axis_off()
length_nm = 200
nm = 1e-9
pixelsize = (my_experiment.imaging_modalities["SMLM"]["detector"]["scale"] / nm) * my_experiment.imaging_modalities["SMLM"]["detector"]["pixelsize"]
pixelsize = np.ceil(pixelsize)
length_px = length_nm / pixelsize
hight_px = 5 
scalebar = AnchoredSizeBar(ax.transData,
                           length_px, 
                           "",
                           "lower right",
                           borderpad=0.5,
                           color='white',
                           frameon=False,
                           size_vertical=hight_px)
ax.add_artist(scalebar)

filename = my_experiment.date_as_string + 'vlab4mic_ccs_DOME_image_simulation.png'
filename2 = os.path.join(my_experiment.output_directory, filename)
fig.savefig(filename2,dpi=300, bbox_inches='tight')
plt.close()

########################## FLAT MODEL #########################

my_experiment.select_structure(
    structure_path=flat_model,
    build=False
)
my_experiment.set_virtualsample_params(
    random_rotations=True,
    rotation_angles=[85],
)
my_experiment.build()
images, noiseless = my_experiment.run_simulation()

fig = plt.figure(figsize=[10,20])
ax = fig.add_subplot(1, 3, 1, projection="3d")
total = my_experiment.structure.num_assembly_atoms
fraction = 50000/total
my_experiment.structure.show_target_labels(
    with_assembly_atoms=True, 
    axis_object=ax, 
    show_axis=False, 
    assembly_fraction = fraction,
    view_init=[90,0,0],
    target_size = 0,
    atoms_size = 1,
    atoms_alpha = 0.01,
    reference_point = False,
    )
ax = fig.add_subplot(1, 3, 2, projection="3d")
my_experiment.structure.show_target_labels(
    with_assembly_atoms=False, 
    axis_object=ax, 
    show_axis=False, 
    view_init=[90,0,0],
    target_size = 10,
    reference_point = False,
    )
ax = fig.add_subplot(1, 3, 3, projection="3d")
my_experiment.particle.gen_axis_plot(
 axis_object=ax,
 view_init=[90,0,0],
 emitter_plotsize=5,
 xlim=[-1250,1250], ylim=[-1250,1250], zlim=[0,1000],
)

filename = my_experiment.date_as_string + 'vlab4mic_ccs_FLAT_labelled_structure.png'
filename2 = os.path.join(my_experiment.output_directory, filename)
fig.savefig(filename2,dpi=300, bbox_inches='tight')
plt.close()

fig = plt.figure(figsize=[20,30])
for i, mod_name in enumerate(images.keys()):
    ax = fig.add_subplot(1,5,i+1)
    ax.imshow(images[mod_name]["ch0"][0], cmap="grey")
    ax.set_title(mod_name)
    ax.set_axis_off()
length_nm = 200
nm = 1e-9
pixelsize = (my_experiment.imaging_modalities["SMLM"]["detector"]["scale"] / nm) * my_experiment.imaging_modalities["SMLM"]["detector"]["pixelsize"]
pixelsize = np.ceil(pixelsize)
length_px = length_nm / pixelsize
hight_px = 5 
scalebar = AnchoredSizeBar(ax.transData,
                           length_px, 
                           "",
                           "lower right",
                           borderpad=0.5,
                           color='white',
                           frameon=False,
                           size_vertical=hight_px)
ax.add_artist(scalebar)

filename = my_experiment.date_as_string + 'vlab4mic_ccs_FLAT_image_simulation.png'
filename2 = os.path.join(my_experiment.output_directory, filename)
fig.savefig(filename2,dpi=300, bbox_inches='tight')
plt.close()
