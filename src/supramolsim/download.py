import os
from pathlib import Path

import requests
import yaml
from tqdm import tqdm

headers = {
    "User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/102.0.0.0 Safari/537.36",
    "Accept": "application/json, text/plain, */*",
    "Accept-Language": "en-GB",
    # 'Accept-Encoding': 'gzip, deflate, br',
    "Referer": "https://henriqueslab.org/",
    "Origin": "https://henriqueslab.org",
    "Connection": "keep-alive",
    "Sec-Fetch-Dest": "empty",
    "Sec-Fetch-Mode": "cors",
    "Sec-Fetch-Site": "same-site",
    "Pragma": "no-cache",
    "Cache-Control": "no-cache",
}


def download_suggested_structures(data_path: str = "data") -> None:
    """
    Downloads PDB *.cif.gz files for the suggested structures in the given path.

    Args:
        data_path: The path to the data directory. Defaults to "data".

    Returns:
        None
    """

    structures_path = Path(data_path) / "structures"

    for filename in structures_path.glob("*.yaml"):
        if filename.stem == "_template":
            continue

        cifs_filename = filename.with_suffix(".cif.gz")

        with open(filename, "r") as f:
            data = yaml.load(f, Loader=yaml.FullLoader)

        id = data["id"]
        title = data["title"]

        if cifs_filename.exists():
            print(f"{title} - [{cifs_filename}] already exists. Skipping...")
            continue

        print(f"Downloading {title} - [{cifs_filename}]...")
        url = f"https://files.rcsb.org/download/{id}.cif.gz"
        download_file(url, filename.with_suffix(".cif.gz"))


def download_file(url: str, fname: str, chunk_size: int = 1024) -> None:
    """Downloads a file from the given URL and saves it to the specified file name.

    Args:
        url: The URL of the file to download.
        fname: The file name to save the downloaded file as.
        chunk_size: The size of the chunks to download the file in. Defaults to 1024.

    Returns:
        None
    """
    resp = requests.get(url, stream=True, headers=headers)
    total = int(resp.headers.get("content-length", 0))

    desc = os.path.basename(fname)
    # Truncate description if too long.
    if len(desc) > 50:
        desc = desc[:30] + "..."

    with (
        open(fname, "wb") as file,
        tqdm(
            desc=desc,
            total=total,
            unit="iB",
            unit_scale=True,
            unit_divisor=1024,
        ) as bar,
    ):
        for data in resp.iter_content(chunk_size=chunk_size):
            size = file.write(data)
            bar.update(size)
