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
import os

# SatNav imports
from mathUtils import *
from datafile import DataFile
from wgs84 import wgs84Constants, WGS84, OrbitalMath
from time_conv import *


DIR = 'data/lab3/'


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
        eng.retrieveNavigationData(float(gpsWeek), float(gpsSec), utc_str, 0, navdir, nargout=0)
        eng.quit()
    except Exception as e:
        print("\nSomething went wrong with MATLAB... here is the exception")
        print(e)


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
    if type(gpsSec) == float and usefile == None:
        navFile = 'gpsW_SOW_{}_{:.4f}_navData.csv'.format(int(gpsWeek), gpsSec)
    elif type(gpsSec) != float and usefile == None:
        navFile = 'gpsW_SOW_{}_{}_navData.csv'.format(int(gpsWeek), gpsSec)
    elif usefile != None:
        navFile = usefile
    
    navDF   = pd.read_csv(navdir + navFile)
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
    R_ecef = np.array([x_ecef, y_ecef, z_ecef]).reshape((3,1)) * u.meter
    V_ecef = np.array([Vx_ecef, Vy_ecef, Vz_ecef]).reshape((3,1)) * u.meter/u.second

    return R_ecef, V_ecef


def satelaz(r_sv_ecef, r_rx_ecef):
    """
    Compute satellite azimuth and elevation angle with respect to the
    receiver in ECEF.

    Inputs:
        r_sv_ecef (m) : 3x1 position vector of satellite in ECEF
        r_rx_ecef (m) : 3x1 position vector of receiver in ECEF
    Outputs:
        az (rad) : azimuth angle 
        el (rad) : elevation angle
    """

    # Compute satellite wrt receiver in ECEF
    r_sv_rx_ecef = r_sv_ecef - r_rx_ecef

    # transform relative vector to ENU
    wgs84       = WGS84()
    lat, lon, h = wgs84.ecef_to_geodetic(r_rx_ecef)
    T_ecef_enu  = wgs84.ecef_to_enu(lat, lon)
    r_sv_rx_enu = np.matmul(T_ecef_enu, r_sv_rx_ecef)

    # Compute azimuth and elevation
    east  = r_sv_rx_enu[0,0]
    north = r_sv_rx_enu[1,0]
    up    = r_sv_rx_enu[2,0]

    az = np.arctan2(east, north)
    el = (0.5*np.pi * u.radian) - np.arccos( up / np.linalg.norm(r_sv_rx_enu) )

    return az, el


def satmap(navFile, r_rx_ecef, el_mask_deg, gpsWeek, gpsSecVec, navdir=DIR, plot_flag=False):
    """
    Generate plotting data for SV's above a particular
    receiver position over a span of GPS seconds.

    Inputs:
        navFile          : csv navigation file
        r_rx_ecef  (m)   : 3x1 receiver position in ECEF
        el_mas_deg (deg) : elevation cutoff
        gpsWeek          : Gps week number
        gpsSecVec (s)    : Gps seconds array
        plot_flag        : flag to create sky plot (plot if True)
    Outputs:
        svIds   : Unique SV ID numbers to be plotted
        sv_data : Nt*Nsv by 4 array in the form
                [svId, gpsSec, az, el
                 svId, gpsSec, az, el
                 .
                 .
                 svId, gpsSec, az, el]
    """
    # 32 SVIDs in each navfile
    SVIDs = range(1,33)
    ephem = pd.read_csv(navdir + navFile)
    elrad = np.deg2rad(el_mask_deg) * u.radian

    # sats_in_view will hold all svids in view at each GpsSec
    svIds        = []
    sv_data      = np.zeros(4)
    for gpsSec in gpsSecVec:
        for sv in SVIDs:
            if np.any(ephem['SVID'] == sv):
                # grab sat ecef position 
                R_sat, V_sat = satloc(gpsWeek, gpsSec, sv, usefile=navFile)
                
                # grab azimuth, elevation from receiver
                az, el = satelaz(R_sat, r_rx_ecef)

                # check if equal to or above elevation threshold
                if el >= elrad:
                    if sv not in svIds:
                        svIds.append(sv)
                    
                    data    = [sv, gpsSec, el.value, az.value]
                    sv_data = np.vstack((sv_data, data))
    
    # delete first row, it was just used as a initializer
    sv_data = np.delete(sv_data, 0, 0)

    if plot_flag:
        # import matlab runner
        sys.path.insert(1, '/Applications/MATLAB_R2020a.app')
        import matlab.engine

        try:
            eng     = matlab.engine.start_matlab()
            satdata = ndarray2matlab(sv_data)
            eng.plotsat(satdata, float(gpsWeek), float(elrad), nargout=0)
            eng.quit()

        except Exception as e:
            print("\nSomething went wrong with MATLAB... here is the exception\n")
            print(e)

    return svIds, sv_data


def channel2navsol(gpsWeek, gpsSec, svID, sec_n=None, rx_ecef=None, createf=False):
    # first, get the data from the matlab script
    if createf:
        createNavDataCsv(gpsWeek, gpsSec)

    # grab that file to get position and then delete later
    usesec = gpsSec
    if sec_n != None:
        usesec = sec_n

    navFile='filler'
    if type(gpsSec) == float:
        navFile = 'gpsW_SOW_{}_{:.4f}_navData.csv'.format(int(gpsWeek), usesec)
    else:
        navFile = 'gpsW_SOW_{}_{}_navData.csv'.format(int(gpsWeek), usesec)

    # get position of SV
    r_sv_ecef, v_ecef = satloc(gpsWeek, gpsSec, svID, usefile=navFile)

    if rx_ecef != None:
        az, el = satelaz(r_sv_ecef, rx_ecef)

        return r_sv_ecef, el
    
    return r_sv_ecef
    

def main():
    print('Using satloc functionality...')

if __name__ == "__main__":
    main()