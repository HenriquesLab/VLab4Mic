from datetime import datetime
import numpy as np
import os

def write_txt(list_of_lines: list, dirname: str, notes: str):
    now = datetime.now()  # dd/mm/YY H:M:S
    sim_dt_string = now.strftime("%Y%m%d") + "_"
    filename = dirname + sim_dt_string + notes + ".csv"
    current_files = os.listdir(dirname)
    if filename.split("/")[-1] in current_files:
        count = 1
        new_filename = (
            dirname
            + sim_dt_string
            + notes
            + "_"
            + str(count)
            + ".csv"
        )
        while new_filename.split("/")[-1] in current_files:
            count += 1
            new_filename = (
                dirname
                + sim_dt_string
                + notes
                + "_"
                + str(count)
                + ".csv"
            )
        filename = new_filename
    text_file = open(filename, "w")
    for li in range(len(list_of_lines)):
        text_file.write(list_of_lines[li])
        text_file.write("\n")
    text_file.close()
    print(f"csv written: {filename}")
