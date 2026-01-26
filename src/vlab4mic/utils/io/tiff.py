import tifffile as tif
from datetime import datetime
import numpy as np
import os


def write_tif(image_stack: np.ndarray, dirname: str, notes: str):
    now = datetime.now()  # dd/mm/YY H:M:S
    dt_string = now.strftime("%Y%m%d") + "_"
    # adjust intensities to match pixel depth
    if image_stack.max() < (2**16 - 1):
        image_stack = image_stack.astype(np.uint16)
    else:
        image_stack = image_stack.astype(np.uint32)
    filename = dirname + dt_string + notes + ".tif"
    current_files = os.listdir(dirname)
    if filename.split("/")[-1] in current_files:
        count = 1
        new_filename = (
            dirname
            + dt_string
            + notes
            + "_"
            + str(count)
            + ".tif"
        )
        while new_filename.split("/")[-1] in current_files:
            count += 1
            new_filename = (
                dirname
                + dt_string
                + notes
                + "_"
                + str(count)
                + ".tif"
            )
        filename = new_filename
    print(f"image saved in {filename}")
    tif.imwrite(filename, image_stack, bigtiff=True)
