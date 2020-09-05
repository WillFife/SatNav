"""
Constants/formulae given by the WGS84 reference system.

@author: William Fife (wnf222)
"""


# import dependencies
import astropy.units as u
import astropy.constants as const
import numpy as np
from numpy import sin, cos, sqrt

# import SatNav modules
from mathUtils import mathUtils


# Constants
E_SMA      = 6378137.0 * u.meter                     # SMA for earth shape
E_FLAT_INV = 298.257223563
E_FLAT     = 1./E_FLAT_INV
E_ROT_VEL  = 7292115e-11 * u.rad / u.second
E_MU       = 398600441800000.0 * const.GM_earth.unit 
E_ECC_SQ   = 2. * E_FLAT - E_FLAT**2


class WGS84:
	self.mathUtils    = mathUtils()
	self.unit_checker = self.mathUtils.check_units

	def R_N(self, lat_gd):
		"""
		Compute and return meridian radius of curvature.

		Inputs:
			Astropy Quantity lat_gd : Geodetic latitude
		Returns:
			Astropy Quantity meridian_roc : Meridian radius of curvature
		"""
		# check units
		lgd          = self.unit_checker(lat_gd, u.radian)
		meridian_roc = E_SMA/sqrt(1. - E_ECC_SQ*sin(lgd))
		return meridian_roc


	def geodetic_point_in_ecef(self, lat_gd, lon, h):
		"""
		Compute ECEF position of a geodetic position.

		Inputs:
			Astropy Quantity lat_gd : Geodetic latitude
			Astropy Quantity lon    : Latitude
			Astropy Quantity h      : Height
		Returns:
			Astropy Quantity r_ecef (meter) : 1x3 array for ECEF position
		"""

		# check units
		lat_gd = self.unit_checker(lat_gd, u.radian)
		lon    = self.unit_checker(lon, u.radian)
		h      = self.unit_checker(h, u.meter)

		r_n = R_N(lat_gd)
		x   = (r_n + h) * cos(lat_gd) * cos(lon)
		y   = (r_n + h) * cos(lat_gd) * sin(lon)
		z   = (r_n*(1 - E_ECC_SQ) + h)* sin(lat_gd)

		r_ecef = np.array([x, y, z]) * u.meter

		return r_ecef


	def ecef_to_geodetic(self, r_ecef):
		pass







