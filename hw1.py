"""
Homework 1 code.

@author: William Fife
"""


# import dependencies
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import astropy.units as u
import sys
sys.path.append('/Users/williamfife/Projects/dev/NavPod')

# SatNav imports
from mathUtils import MathUtils
from datafile import DataFile
from wgs84 import wgs84Constants, WGS84

# NavPod imports
from plotter import Plotter

FIG_DIR = '/Users/williamfife/Documents/SchoolandWork/Academics/Fall2020/SatNav/Reports/HW1/plots/'
plotter = Plotter()
plotter.fig_dir = FIG_DIR
wgs84 = WGS84()

def run_files(filelist):
    fig_s  = plt.figure('pos_deltas')
    ax_pos = fig_s.add_subplot(111)
    fig_s.suptitle('delta_lat_lon')
    ax_pos.set(ylabel='dlat (deg)')
    ax_pos.set(xlabel='dlon (deg)')
    ax_pos.grid()
    ax_pos.axis('equal')
    ax_pos.set_ylim(-5e-5, 5e-5)
    ax_pos.set_xlim(-5e-5, 5e-5)

    fig_en  = plt.figure('enu_scatter')
    ax_en = fig_en.add_subplot(111)
    fig_en.suptitle('enu_scatter')
    ax_en.set(ylabel='dNorth (m)')
    ax_en.set(xlabel='dEast (m)')
    ax_en.grid()
    ax_en.axis('equal')
    ax_en.set_ylim(-10, 10)
    ax_en.set_xlim(-10, 10)

    fig_p  = plt.figure('delta_pos_enu')
    ax_e = fig_p.add_subplot(311)
    fig_p.suptitle('delta_pos_enu')
    ax_e.set(ylabel='dE (m)')
    ax_e.set(xlabel='Time (s)')
    ax_e.grid()
    ax_e.set_ylim(-10, 10)

    ax_n = fig_p.add_subplot(312)
    ax_n.set(ylabel='dN (m)')
    ax_n.set(xlabel='Time (s)')
    ax_n.grid()
    ax_n.set_ylim(-10, 10)

    ax_u = fig_p.add_subplot(313)
    ax_u.set(ylabel='dU (m)')
    ax_u.set(xlabel='Time (s)')
    ax_u.grid()
    ax_u.set_ylim(-10, 10)
    for f in filelist:
        print('---- WORKING '+f+' ----')
        fname = f.split('.')[0]
        df = DataFile(f, data_dir='data/hw1/')
        df.setup()
        lat_mean  = np.mean(df.file_df['latitude {deg}']) * u.degree
        lon_mean  = np.mean(df.file_df['longitude {deg}']) * u.degree
        h_mean    = np.mean(df.file_df['h_ellipsiod {m}']) * u.meter
        sats_mean = np.mean(df.file_df['num_sats'])

        dlats = df.file_df['latitude {deg}'] - lat_mean.value
        dlons = df.file_df['longitude {deg}'] - lon_mean.value
        dhs   = df.file_df['h_ellipsiod {m}'] - h_mean.value

        dhs_std = np.std(dhs)
        lat_std = np.std(df.file_df['latitude {deg}'])
        lon_std = np.std(df.file_df['longitude {deg}'])
        CEP_latlon = 0.62*lat_std + 0.56*lon_std

        res = ax_pos.scatter(dlons, dlats, s=10, label=fname)
        color  = res.get_facecolor()[0]
        circle = plt.Circle((0, 0), CEP_latlon, fill=False, color=color)
        ax_pos.add_artist(circle)

        print('-- Means -- lat = {}, lon = {}, h_mean = {}, num_sats = {}'.format(lat_mean, lon_mean, h_mean, sats_mean))

        # compute ECEF position vector using mean LAT, LON, H
        r_ecef_mean = wgs84.geodetic_point_in_ecef(lat_mean,
                                                   lon_mean,
                                                   h_mean)

        rs_ecef = wgs84.geodetic_point_in_ecef_arr(df.file_df['latitude {deg}'].to_numpy() * u.degree,
                                                   df.file_df['longitude {deg}'].to_numpy() * u.degree,
                                                   df.file_df['h_ellipsiod {m}'].to_numpy() * u.meter)
        
        x_ecef_std = np.std(rs_ecef[:,0])
        y_ecef_std = np.std(rs_ecef[:,1])
        z_ecef_std = np.std(rs_ecef[:,2])
        SEP = 0.51*(x_ecef_std + y_ecef_std + z_ecef_std)
        print('-- SEP -- {} {}'.format(SEP, u.meter.to_string()))
        print('-- dhs_std -- {} {}'.format(dhs_std, u.meter.to_string()))
        print('-- CEP_latlon -- {} {}'.format(CEP_latlon, u.degree.to_string()))

        print('ECEF_from_MEAN_GEO = ', r_ecef_mean)

        R_ecef_enu = wgs84.ecef_to_enu(lat_mean, lon_mean)
        rs_enu = np.zeros((len(rs_ecef), 3))

        for ii in range(len(rs_ecef)):
            r_ecef  = rs_ecef[ii, :]
            r_ecef_row = r_ecef_mean.flatten()
            dr_ecef = r_ecef - r_ecef_row.value
            dr_enu  = np.matmul(R_ecef_enu, dr_ecef)
            rs_enu[ii, :] = dr_enu

        ax_e.plot(df.time_s, rs_enu[:,0], label='SEP={:.2f} {}'.format(SEP, u.meter.to_string()))
        ax_n.plot(df.time_s, rs_enu[:,1], label=f)
        ax_u.plot(df.time_s, rs_enu[:,2], label=fname)

        ax_en.scatter(rs_enu[:,0], rs_enu[:,1], s=10, label=fname)


        # newline
        print("")
    
    ax_pos.legend()
    ax_e.legend()
    ax_u.legend()
    ax_en.legend()


if __name__ == "__main__":
    filelist = ['nmea20200828Jahaus.csv', 'nmea20200831Jahaus.csv', 'nmeaJahaus20200807.csv', 'nmeaJahaus20200909.csv']
    run_files(filelist)
    plt.show()