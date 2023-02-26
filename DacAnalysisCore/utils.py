from .exceptions import AnalyticalException
import numpy as np
from .ddn import DictDotNotation as ddn


def change_time_unit(data, orig_unit, chng_unit):
    """
    time data unit converter
    
    args:
    - data : 1 dimentional list, must be time infomation.
    - orig_unit : unit of data provided. "s" or "m" or "h"
    - chng_unit : unit of data returns. "s" or "m" or "h"
    
    returns:
    - data : converted 1 dimentional list

    raises:
    - AnalyticalException : When orig_unit or chng_unit is not in "s" or "m" or "h"
    """
    if orig_unit == chng_unit:
        return data
    try:
        on = {"s": 0, "m": 1, "h":2}[orig_unit]
        cn = {"s": 0, "m": 1, "h":2}[chng_unit]
        conv = [1, 60, 3600, 1/3600, 1/60][on-cn]
    except IndexError:
        raise AnalyticalException("unexpected argument")
    else:
        return [i*conv for i in data]


def change_conc_unit(data, orig_unit, chng_unit):
    """
    conc data unit converter
    
    args:
    - data : 1 dimentional list, must be concentration infomation.
    - orig_unit : unit of data provided. "%" or "ppm"
    - chng_unit : unit of data returns. "%" or "ppm"
    
    returns:
    - data : converted 1 dimentional list

    raises:
    - AnalyticalException : When orig_unit or chng_unit is not in "%" or "ppm"
    """
    if orig_unit == chng_unit:
        return data
    elif [orig_unit, chng_unit] == ["ppm", "%"]:
        return [i if i is None else i/10000 for i in data]
    elif [orig_unit, chng_unit] == ["%", "ppm"]:
        return [i if i is None else i*10000 for i in data]
    else:
        raise AnalyticalException("unexpected arguments")


def get_range(values):
    """
    calculate dynamic range of data
    
    args:
    - values : 1 dimentional list
    
    returns:
    - ddn : DictDotNotation that contains min, max, range infomation
    """
    values = [i for i in values if i is not None]
    mx = values[np.argmax(values)]
    mn = values[np.argmin(values)]
    return ddn({
        "max": mx,
        "min": mn,
        "diff": mx-mn
    })