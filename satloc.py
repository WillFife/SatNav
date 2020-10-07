"""
Retrieve the CSV file for navigation data for the
32 possible satellites at a time given by GPS week and GPS seconds.

Also create a CSV of the ECEF position and velocities from the given navigation data.

@author: William Fife
"""


# import dependencies
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import astropy.units as u
import sys

# SatNav imports
from mathUtils import MathUtils
from datafile import DataFile
from wgs84 import wgs84Constants, WGS84, OrbitalMath
from time_conv import *


DIR = 'data/lab2/'


def createNavDataCsv(gpsWeek, gpsSec, navdir=DIR[:-1]):
    print('--- Creating navData csv for GPS week {}, GPS Sec {} s ---'.format(gpsWeek, gpsSec))

    # import matlab runner for extracting nav data
    sys.path.insert(1, '/Applications/MATLAB_R2020a.app')
    import matlab.engine

    # Compute UTC time to give to matlab function 'RetrieveNavigationData.m'
    utc = gps2utc(gpsWeek, gpsSec)
    utc_str = '{}'.format(str(utc))
    print(utc_str)
    
    # Access matlab function and create csv
    try:
        eng = matlab.engine.start_matlab()
        eng.retrieveNavigationData(gpsWeek, gpsSec, utc_str, 0, navdir)
        eng.quit()
    except ValueError:
        print("\nFile created. Exception from Matlab return type handled.")


def satloc(gpsWeek, gpsSec, svID, navdir=DIR, usefile=None, printSat=False):
    """
    Return satellite position and velocity expressed in and relative
    to the ECEF reference frame.
    
    Inputs:
        gpsWeek - Week of true time at which satellite ECEF state is desired

        gpsSec  - Seconds of week

        svID    - The satellite PRN number

    Outputs:
        r_sv_ecef - Position of satellite in ECEF [m]
        v_sv_ecef - Velocity of satellite in ECEF [m/s]
    """

    # Formulate filename and retrieve nav data
    navFile='filler'
    if type(gpsSec) == float:
        navFile = 'gpsW_SOW_{}_{:.4f}_navData.csv'.format(int(gpsWeek), gpsSec)
    else:
        navFile = 'gpsW_SOW_{}_{}_navData.csv'.format(int(gpsWeek), gpsSec)
    navDF   = pd.read_csv(navdir + navFile)
    if usefile != None:
        navDF = pd.read_csv(navdir+usefile)
    sat     = navDF[navDF['SVID']==svID]

    # load in variables needed
    GM      = 3.986005e14            # Earth gravitational parameter m3/s2
    OmegaE  = 7.2921151467e-5        # Earth mean rotation rate rad/s
    t_oe    = float(sat['te'])       # ephemeris epoch
    t       = float(sat['tc'])       # clock time
    A       = float(sat['sqrta']**2)
    dn      = float(sat['dn'])
    M0      = float(sat['M0'])
    a0      = float(sat['af0'])
    a1      = float(sat['af1'])
    a2      = float(sat['af2'])
    e       = float(sat['e'])
    i0      = float(sat['i0'])
    L0      = float(sat['L0'])
    i_dot   = float(sat['didt'])
    lan_dot = float(sat['dOdt'])
    omega0  = float(sat['omega0'])
    C_uc    = float(sat['Cuc'])
    C_us    = float(sat['Cus'])
    C_rc    = float(sat['Crc'])
    C_rs    = float(sat['Crs'])
    C_ic    = float(sat['Cic'])
    C_is    = float(sat['Cis'])

    if printSat:
        print(sat)

    # initial computations for perifocal state
    dtc  = a0 + a1*(gpsSec - t) + a2*((gpsSec - t)**2)
    tc   = t_oe - dtc
    tk   = gpsSec - tc
    n0   = np.sqrt(GM / A**3)
    n    = n0 + dn
    M    = M0 + n*tk
    E    = M

    # Kepler's equation Newton's method for Eccentric Anomaly
    orbmath = OrbitalMath()
    # Only 20 iterations used
    E = orbmath.E_newton(M, e, max_iter=20)

    v_k    = orbmath.true_anomaly_from_E(E, e) # true anomaly
    argl_k = v_k + omega0
    
    # correction terms
    dargl = C_us*np.sin(2*argl_k) + C_uc*np.cos(2*argl_k)
    dr    = C_rs*np.sin(2*argl_k) + C_rc*np.cos(2*argl_k)
    dinc  = C_is*np.sin(2*argl_k) + C_ic*np.cos(2*argl_k)

    # corrected terms
    argl = argl_k + dargl
    r_k  = A*(1 - e*np.cos(E)) + dr
    i_k  = i0 + dinc + i_dot*tk

    # position in perifocal frame
    p_x = r_k*np.cos(argl)
    p_y = r_k*np.sin(argl)

    # corrected longitude of ascending node
    lan_k = L0 + (lan_dot - OmegaE)*tk - OmegaE * t_oe

    # ECEF position
    x_ecef = p_x*np.cos(lan_k) - p_y*np.cos(i_k)*np.sin(lan_k)
    y_ecef = p_x*np.sin(lan_k) + p_y*np.cos(i_k)*np.cos(lan_k)
    z_ecef = p_y*np.sin(i_k)

    # intermediate terms for satellite velocity
    E_k_dot   = n / (1 - e*np.cos(E))
    v_k_dot   = E_k_dot*np.sqrt(1-e**2) / (1 - e*np.cos(E))
    i_k_dot   = i_dot + 2*v_k_dot*(C_is*np.cos(2*argl_k) - C_ic*np.sin(2*argl_k))
    argl_dot  = v_k_dot + 2*v_k_dot*(C_us*np.cos(2*argl_k) - C_uc*np.sin(2*argl_k))
    r_k_dot   = e*A*E_k_dot*np.sin(E) + 2*v_k_dot*(C_rs*np.cos(2*argl_k) - C_rc*np.sin(2*argl_k))
    lan_k_dot = lan_dot - OmegaE

    # perifocal velocity
    p_vx = r_k_dot*np.cos(argl_k) - r_k*argl_dot*np.sin(argl_k)
    p_vy = r_k_dot*np.sin(argl_k) + r_k*argl_dot*np.cos(argl_k)

    print('p_x = ', p_x)
    print('p_vx = ', p_vx)
    print('lan_k = ', lan_k)

    # ECEF velocity
    Vx_ecef = -p_x*lan_k_dot*np.sin(lan_k) + \
              p_vx*np.cos(lan_k) - \
              p_vy*np.sin(lan_k)*np.cos(i_k) - \
              p_y*( lan_k_dot*np.cos(lan_k)*np.cos(i_k) - i_k_dot*np.sin(lan_k)*np.sin(i_k) )
    
    Vy_ecef = p_x*lan_k_dot*np.cos(lan_k) + \
              p_vx*np.sin(lan_k) + \
              p_vy*np.cos(lan_k)*np.cos(i_k) - \
              p_y*( lan_k_dot*np.sin(lan_k)*np.cos(i_k) + i_k_dot*np.cos(lan_k)*np.sin(i_k) )

    Vz_ecef = p_vy*np.sin(i_k) + p_y*i_k_dot*np.cos(i_k)

    # place into arrays and return
    R_ecef = np.array([x_ecef, y_ecef, z_ecef])
    V_ecef = np.array([Vx_ecef, Vy_ecef, Vz_ecef])

    return R_ecef, V_ecef


def main():
    gpsWeek_t, gpsSec_t = (1653.0, 570957.338101566)
    print('---TESTING SATLOC AGAINST GPSweek {}, GPSsec {}'.format(gpsWeek_t, gpsSec_t))

    pos_t, vel_t = satloc(gpsWeek_t, gpsSec_t, 1, usefile='gpsW_SOW_1653_570957.3381_navData.csv')

    print('\n Test position {}, \nTest velocity = {}'.format(pos_t.reshape((3,1)), vel_t.reshape((3,1))))
    
    datetimeformat = "%Y-%m-%d %H:%M:%S"
    utc = datetime.strptime("2013-09-01 12:00:00",datetimeformat)
    gpsWeek, gpsSec = utc2gps(utc)
    print('\n---TESTING SATLOC AGAINST GPSweek {}, GPSsec {}'.format(gpsWeek, gpsSec))

    pos_2, vel_2 = satloc(gpsWeek, gpsSec, 2, usefile='gpsW_SOW_1756_43165_navData.csv')
    pos_5, vel_5 = satloc(gpsWeek, gpsSec, 5, usefile='gpsW_SOW_1756_43165_navData.csv')

    print('\n --PRN 2-- Test position {}, \nTest velocity = {}'.format(pos_2.reshape((3,1)), vel_2.reshape((3,1))))
    print('\n --PRN 5-- Test position {}, \nTest velocity = {}'.format(pos_5.reshape((3,1)), vel_5.reshape((3,1))))

    utc = datetime.strptime("2013-09-01 12:00:01",datetimeformat)
    gpsWeek, gpsSec = utc2gps(utc)
    print('\n---TESTING SATLOC AGAINST GPSweek {}, GPSsec {}'.format(gpsWeek, gpsSec))

    pos_2t, vel_2t = satloc(gpsWeek, gpsSec, 2, usefile='gpsW_SOW_1756_43166_navData.csv')
    pos_5t, vel_5t = satloc(gpsWeek, gpsSec, 5, usefile='gpsW_SOW_1756_43166_navData.csv')

    print('\n --PRN 2-- Test position {}, \nTest velocity = {}'.format(pos_2t.reshape((3,1)), vel_2t.reshape((3,1))))
    print('\n --PRN 5-- Test position {}, \nTest velocity = {}'.format(pos_5t.reshape((3,1)), vel_5t.reshape((3,1))))


if __name__ == "__main__":
    main()