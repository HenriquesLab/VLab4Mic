from datetime import datetime
import numpy as np


def write_txt(list_of_lines: list, dirname: str, notes: str):
    now = datetime.now()  # dd/mm/YY H:M:S
    sim_dt_string = now.strftime("%Y%m%d") + "_"
    filename = dirname + sim_dt_string + notes + ".csv"
    text_file = open(filename, "w")
    for li in range(len(list_of_lines)):
        text_file.write(list_of_lines[li])
        text_file.write("\n")
    text_file.close()
    print(f"csv written: {filename}")
