"""
Constants/formulae given by the WGS84 reference system.

@author: William Fife (wnf222)
"""

# import dependencies
import astropy.units as u
import astropy.constants as const
import numpy as np
from numpy import sin, cos, sqrt

# Constants
E_SMA      = 6378137.0 * u.meter                     # Not SMA for earth point mass about Sun
E_FLAT_INV = 298.257223563
E_FLAT     = 1./E_FLAT_INV
E_ROT_VEL  = 7292115e-11 * u.rad / u.second
E_MU       = 398600441800000.0 * const.GM_earth.unit 
E_ECC_SQ   = 2. * E_FLAT - E_FLAT**2


def R_N(lat_gd):
	if lat_gd.unit.is_equivalent(u.degree):
		return E_SMA/sqrt(1. - E_ECC_SQ*sin(lat_gd.to(u.degree)))

	return E_SMA/sqrt(1. - E_ECC_SQ*sin(lat_gd))


def geodetic_point_in_ecef(lat_gd, lon, h):
	"""
	Compute ECEF position of a geodetic position.

	Inputs:
		Astropy Quantity lat_gd : Geodetic latitude
		Astropy Quantity lon    : Latitude
		Astropy Quantity h      : Height

	Returns:
		Astropy Quantity r_ecef (meter) : 1x3 array for ECEF position
	"""

	# convert units
	if lat_gd.unit.is_equivalent(u.radian) not True:
		lat_gd = lat_gd.to(u.radian)
	if lon.unit.is_equivalent(u.radian) not True:
		lon = lon.to(u.radian)
	if h.unit.is_equivalent(u.meter) not True:
		h = h.to(u.meter)

	R_N = R_N(lat_gd)
	x   = (R_N + h) * cos(lat_gd) * cos(lon)
	y   = (R_N + h) * cos(lat_gd) * sin(lon)
	z   = (R_N*(1 - E_ECC_SQ) + h)* sin(lat_gd)

	r_ecef = np.array([x, y, z]) * u.meter

	return r_ecef


def ecef_to_geodetic(r_ecef):
	pass







