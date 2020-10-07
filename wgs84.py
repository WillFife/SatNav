"""
Constants/formulae given by the WGS84 reference system.

@author: William Fife (wnf222)
"""


# import dependencies
import astropy.units as u
import astropy.constants as const
import numpy as np
from numpy import sin, cos, sqrt
import pandas as pd

# import SatNav modules
from mathUtils import MathUtils


class wgs84Constants():
	def __init__(self):
		# Constants
		self.E_SMA      = 6378137.0 * u.meter                     # SMA for earth shape
		self.E_FLAT_INV = 298.257223563
		self.E_FLAT     = 1./self.E_FLAT_INV
		self.E_ROT_VEL  = 7292115e-11 * u.rad / u.second
		self.E_MU       = 398600441800000.0 * const.GM_earth.unit 
		self.E_ECC_SQ   = 2. * self.E_FLAT - self.E_FLAT**2


class WGS84():
	def __init__(self):
		self.MathUtils    = MathUtils()
		self.unit_checker = self.MathUtils.check_units
		self.consts       = wgs84Constants()


	def R_N(self, lat_gd):
		"""
		Compute and return meridian radius of curvature.

		Inputs:
			Astropy Quantity lat_gd : Geodetic latitude
		Returns:
			Astropy Quantity meridian_roc : Meridian radius of curvature
		"""
		meridian_roc = self.consts.E_SMA/sqrt(1. - self.consts.E_ECC_SQ*sin(lat_gd))
		return meridian_roc


	def geodetic_point_in_ecef(self, lat_gd, lon, h):
		"""
		Compute ECEF position of a geodetic position.

		Inputs:
			Astropy Quantity lat_gd : Geodetic latitude
			Astropy Quantity lon    : Latitude
			Astropy Quantity h      : Height above reference ellipsoid
		Returns:
			Astropy Quantity r_ecef (m) : 3x1 array for ECEF position
		"""
		r_n = self.R_N(lat_gd)
		x   = (r_n + h) * cos(lat_gd) * cos(lon)
		y   = (r_n + h) * cos(lat_gd) * sin(lon)
		z   = (r_n*(1 - self.consts.E_ECC_SQ) + h)* sin(lat_gd)

		r_ecef = u.Quantity([x, y, z], u.meter)

		return r_ecef.reshape((3,1))

	
	def geodetic_point_in_ecef_arr(self, lat_gd, lons, hs):
		"""
		Compute ECEF positions from arrays of latitude, longitude, and heights.
		"""
		r_ecef = np.zeros((len(lat_gd), 3))

		for ii in range(len(lat_gd)):
			lat = lat_gd[ii]
			lon = lons[ii]
			h   = hs[ii]
			r   = self.geodetic_point_in_ecef(lat, lon, h)
			r_ecef[ii,:] = r.flatten()

		return r_ecef


	def ecef_to_geodetic(self, r_ecef):
		"""
		Compute geodetic position from ECEF cartesian vector

		Inputs:
			Astropy Quantity r_ecef (m): 3x1 array ECEF position
		Returns:
			Astropy Quantity lat_gd : Geodetic latitude
			Astropy Quantity lon    : Latitude
			Astropy Quantity h      : Height above reference ellipsoid
		"""

		# check units
		r_ecef = r_ecef.to(u.meter)
		
		# unpack and compute intermediate values
		x     = r_ecef[0,0]
		y     = r_ecef[1,0]
		z     = r_ecef[2,0]
		r_mag = np.linalg.norm(r_ecef)
		rho   = np.sqrt(x**2 + y**2)
		
		# compute longitude
		lon = np.arctan2(y, x)

		# initial latitude and height guess
		lat = np.arcsin(z/r_mag)
		h   = rho/np.cos(lat) - self.R_N(lat)

		# loop through to find latitude
		max_iter = 100
		i        = 0
		tol      = 1e-7
		while i < max_iter:
			R_N     = self.R_N(lat)
			numer   = z + R_N*self.consts.E_ECC_SQ*np.sin(lat)
			new_lat = np.arctan2(numer, rho)
			if abs(new_lat - lat) <= tol:
				lat = new_lat
				h   = R_N
				return lat, lon, h
			
			lat = new_lat
			i += 1
		
		return lat, lon, h


	def ecef_to_enu(self, lat_gd, lon):
		"""
		Generate the rotation matrix used to express a vector in ECEF
		as a vector in ENU at the position defined by geodetic
		coordinates.

		Inputs:
			Astropy Quantity lat_gd: Geodetic latitude
			Astropy Quantity lon: Longitude
		Returns:
			Astropy Quantity T_ecef_enu: 3x3 Rotation matrix from ECEF to ENU
		"""
		# check units
		lat_gd = self.unit_checker(lat_gd, u.radian)
		lon    = self.unit_checker(lon, u.radian)

		T_ecef_enu      = np.eye(3)
		T_ecef_enu[0,0] = -np.sin(lon)
		T_ecef_enu[0,1] = np.cos(lon)
		T_ecef_enu[0,2] = 0.
		T_ecef_enu[1,0] = -np.sin(lat_gd)*np.cos(lon)
		T_ecef_enu[1,1] = -np.sin(lat_gd)*np.sin(lon)
		T_ecef_enu[1,2] = np.cos(lat_gd)
		T_ecef_enu[2,0] = np.cos(lon)*np.cos(lat_gd)
		T_ecef_enu[2,1] = np.sin(lon)*np.cos(lat_gd)
		T_ecef_enu[2,2] = np.sin(lat_gd)

		return T_ecef_enu


class OrbitalMath():
	def __init__(self):
		self.var = 'placeholder'


	def E_newton(self, M, e, tol=10e-10, max_iter=20):
		"""
		Use newton-raphson method to solve for E
		"""
		if M < np.pi:
			E_prev = M + e/2
		else:
			E_prev = M - e/2

		for i in range(max_iter):
			E_next = E_prev - (E_prev - e*np.sin(E_prev) - M)/(1 - e*np.cos(E_prev))
			if E_next - E_prev < tol:
				return E_next
			
			E_prev = E_next
			if i==max_iter-1:
				return E_next


	def true_anomaly_from_E(self, E, e):
		"""
		Compute true anomaly from eccentric anomaly
		"""
		num = np.tan(E/2.) * np.sqrt(1 + e)
		den = np.sqrt(1 - e)

		return 2. * np.arctan2(num, den)