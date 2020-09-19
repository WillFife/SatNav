"""
Data file class to load and process data for SatNav.

@author: William Fife
"""


# import dependencies
import numpy as np
import pandas as pd
import astropy.units as u
from mathUtils import MathUtils


class DataFile():
    """
    Data file manipulation class.
    """
    def __init__(self, filename, data_dir='data/'):
        self.file_df = pd.read_csv(data_dir + filename)
        self.mutils  = MathUtils()
        self.time_s  = np.arange(0., len(self.file_df), 1.0) * u.second

        self.col_units = {'latitude {deg}' : u.degree,
                          'longitude {deg}': u.degree,
                          'h_ellipsiod {m}': u.meter}


    def attach_units(self, element, unit):
        return element * unit
    

    def setup(self):
        """
        * Attach units to specified columns in IMU files.
        * Create time series.
        * modify data types in other columns.
        """
        # Attach units
        for col in self.col_units.keys():
            self.file_df[col].astype(np.float64, copy=False)
            self.file_df[col].apply(self.attach_units, args=(self.col_units[col],))