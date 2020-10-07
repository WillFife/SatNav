"""
Time functions and conversions

@author: William Fife
"""


# import dependencies
import numpy as np
from datetime import datetime, timedelta
from leapseconds import *


def utc2gps(n):
    """
    Convert UTC time to GPS time expressed in GPS week
    and GPS second of week
    
    Inputs:
        n -- The UTC date and time expressed as a datetime object
             with format '%Y-%m-%d %H:%M:%S'
        
    Outputs:
        gpsWeek -- The unambiguous GPS week (no rollover at 1024 weeks) [integer]

        gpsSec -- The GPS time of week expressed as GPS seconds
    
    References:
        For the leap seconds -- https://gist.github.com/zed/92df922103ac9deb1a05
    """

    # Create the UTC/GPS epoch in a datetime format
    datetimeformat = "%Y-%m-%d %H:%M:%S"
    epoch = datetime.strptime("1980-01-06 00:00:00",datetimeformat)

    # Calculate the difference between GPS and UTC time
    # accounting for leap seconds
    tdiff = n - epoch - dTAI_UTC_from_utc(n)

    # get gps week from the time difference in days
    gpsWeek = tdiff.days // 7

    # get gps week from the time difference in seconds
    gpsSec  = tdiff.seconds + 86400 * (tdiff.days - 7 * gpsWeek)

    return gpsWeek, gpsSec


def gps2utc(gpsWeek, gpsSec):
    """
    Convert GPS time expressed in GPS week and GPS second of week to UTC
    
    Inputs:
        gpsWeek - unambiguous GPS week [integer]

        gpsSec - GPS time of week expressed in GPS seconds [float]
    
    Outputs:
        n - UTC date and time expressed as a datetime with format
            '%Y-%m-%d %H:%M:%S'
    """

    # Create the epoch in a datetime format
    datetimeformat = "%Y-%m-%d %H:%M:%S"
    epoch          = datetime.strptime("1980-01-06 00:00:00",datetimeformat)

    # determine the leap seconds
    noleap = timedelta(days=(gpsWeek*7), seconds=(gpsSec))
    leap   = dTAI_UTC_from_utc(epoch + noleap)
    
    # compute UTC time
    delta = timedelta(days=(gpsWeek*7), seconds=(gpsSec + leap.seconds))
    n     = datetime.strftime(epoch + delta, datetimeformat)

    return n

    
def main():
    datetimeformat = "%Y-%m-%d %H:%M:%S"
    utc = datetime.strptime("2013-09-01 12:00:00",datetimeformat)
    gps = utc2gps(utc)
    print('UTC : ', utc)
    print('\n GPS Week = {}, GPS Sec = {}'.format(gps[0], gps[1]))
    print('\n gps2utc : ', gps2utc(gps[0], gps[1]))
    print('day of year: ', utc.timetuple().tm_yday)