"""Read and write emitter and localisation tables as CSV files.

The format follows the ThunderSTORM localisation table (comma-separated,
quoted headers with units), which ThunderSTORM, SMAP, Picasso (via its CSV
import), SuReSim and most localisation analysis tools can read.
"""

import csv

import numpy as np

COLUMNS = [
    "id",
    "frame",
    "x [nm]",
    "y [nm]",
    "z [nm]",
    "intensity [photon]",
    "uncertainty [nm]",
]


def write_localisation_table(
    path, coordinates, photons=None, uncertainty=None, frame=None
):
    """
    Write positions to a ThunderSTORM-style CSV table.

    Parameters
    ----------
    path : str
        Output file path.
    coordinates : numpy.ndarray
        Nx3 array of x, y, z positions in nm.
    photons : numpy.ndarray or float, optional
        Photons per position. Empty column if None.
    uncertainty : numpy.ndarray or float, optional
        Lateral localisation uncertainty (standard deviation) in nm.
        Empty column if None (e.g. for true emitter positions).
    frame : numpy.ndarray or int, optional
        Frame of each position (1-based). Default 1.

    Returns
    -------
    str
        The path written.
    """
    coordinates = np.asarray(coordinates, dtype=float).reshape(-1, 3)
    n = coordinates.shape[0]

    def column(values, default):
        if values is None:
            return [default] * n
        values = np.asarray(values)
        if values.ndim == 0:
            return [values.item()] * n
        return list(values)

    frames = column(frame, 1)
    intensities = column(photons, "")
    uncertainties = column(uncertainty, "")
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, quoting=csv.QUOTE_NONNUMERIC)
        writer.writerow(COLUMNS)
        for i in range(n):
            writer.writerow(
                [
                    i + 1,
                    int(frames[i]),
                    round(float(coordinates[i, 0]), 3),
                    round(float(coordinates[i, 1]), 3),
                    round(float(coordinates[i, 2]), 3),
                    intensities[i] if intensities[i] == "" else float(intensities[i]),
                    uncertainties[i] if uncertainties[i] == "" else float(uncertainties[i]),
                ]
            )
    return path


def read_localisation_table(path):
    """
    Read a ThunderSTORM-style CSV table.

    Parameters
    ----------
    path : str
        Table path.

    Returns
    -------
    dict
        Column name -> numpy array. Empty cells become NaN.
    """
    with open(path, newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        rows = [row for row in reader]
    table = {}
    for j, name in enumerate(header):
        table[name] = np.array(
            [float(row[j]) if row[j] != "" else np.nan for row in rows]
        )
    return table
