import math
import numpy as np


def truncate(number, digits) -> float:
    # from Erwin Mayer at
    # https://stackoverflow.com/questions/8595973/truncate-to-three-decimals-in-python
    nbDecimals = len(str(number).split(".")[1])
    if nbDecimals <= digits:
        return number
    stepper = 10.0**digits
    return math.trunc(stepper * number) / stepper


def inCube(X, corners):
    """
    Check if a point `X` is inside of cube with `corners` (two points)
    # Based from: Michael Silverstein's answer
    at
    stackoverflow.com/questions/29720910/
    fastest-way-to-search-if-a-coordinate-is-inside-a-cube
    """
    # Where is X > corners?
    greater = X > corners
    # If X is greater than both corners of any axis, it is outside
    inside = ~np.any(np.equal(*greater))
    return inside


def generate_coordinate_textline(x, y, z, separator=",", decimals=3, **kwargs):
    line = (
        str(truncate(x, decimals))
        + separator
        + str(truncate(y, decimals))
        + separator
        + str(truncate(z, decimals))
    )
    return line


def dictionary2string(mydict):
    str_dict0 = str(mydict)
    str_dict1 = str_dict0.replace(" ", "")
    str_dict2 = str_dict1.replace(",", "_")
    str_dict3 = str_dict2.replace("'", "")
    str_dict4 = str_dict3.replace(":", "")
    str_dict5 = str_dict4.replace(":", "")
    return str_dict5[1:-1] + "_"
