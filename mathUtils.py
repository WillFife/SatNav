"""
Math utilities for efficiency and functionality

@author: William Fife
"""

# import dependencies
import astropy.units as u
import numpy as np
import matlab.engine


def npfloat2float(element):
    return element.item()


def ndarray2matlab(arr):
    #eng1 = matlab.engine.start_matlab()
    n_rows = arr.shape[0]
    n_cols = arr.shape[1]
    n_elem = arr.size
    
    flat      = arr.flatten('F')
    flat_list = [float(i) for i in flat]

    mat = matlab.double(flat_list, size=(n_rows, n_cols))
    return mat


class MathUtils:
    """
    Math utilities class.
    """

    def __init__(self):
        pass
    
    def check_units(self, quantity, unit):
        """
        Check if quantity has unit specified. If not,
        convert quantity to unit.

        Inputs:
            Astropy Quantity quantity: Quantity to have units checked.
            Astropy Unit: Unit to attach to quantity.
        Returns:
            Astropy Quantity: Quantity with correct unit
        """
        try:
            if quantity.unit.is_equivalent(unit) == False:
                corrected_unit = quantity.to(unit)
                return corrected_unit
            return quantity
        except Exception as e:
            raise