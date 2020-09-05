"""
Math utilities for efficiency and functionality

@author: William Fife
"""

# import dependencies
import astropy.units as u
import numpy as np


class mathUtils:
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
            if not quantity.unit.is_equivalent(unit):
                corrected_unit = quantity.to(unit)
                return corrected_unit
        except Exception as e:
            raise