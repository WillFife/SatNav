"""
Constants/formulae given by the WGS84 reference system.

@author: William Fife (wnf222)
"""

# import dependencies
import astropy.units as u
import astropy.constants as const
from numpy import sin

# Constants
E_SMA     = 6378137.0 * u.meter
E_FLAT    = 298.257223563 
E_ROT_VEL = 7292115e-11 * u.rad / u.second
E_MU      = 398600441800000.0 * const.GM_earth.unit
E_ECC_SQ  = 2. * E_FLAT - E_FLAT**2


def R_N(lat_gd):
	pass