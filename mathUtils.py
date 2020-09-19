"""
Math utilities for efficiency and functionality

@author: William Fife
"""

# import dependencies
import astropy.units as u
import numpy as np


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