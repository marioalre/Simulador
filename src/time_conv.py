# Library for time conversion

import numpy as np

def sid2rad(second):
    ''' Sidereal time to radians'''

    hour = second / 3600

    return hour * 15 * np.pi / 180

def rad2sid(radians):
    ''' Radians to sidereal time'''

    hour = radians * 180 / np.pi / 15

    return hour

def ConvTime(date, UTC, dUT1, dAT):
    '''The following formulas provide a means
    to determine the number of centuries elapsed
     from the epoch J2000.0 Ref: Vallado
    
    Parameters
    ----------
    date : list
        Date in the format [year, month, day]
    UTC : string
        UTC time in the format HH:MM:SS
    dUT1 : float
        Difference between UT1 and UTC in seconds
    dAT : float 
        Difference between TAI and UTC in seconds
    Returns
    -------
    '''

    UT1 = float(UTC[0:2]) + float(UTC[3:5]) / 60 + float(UTC[6:8]) / 3600 + dUT1 / 3600 # UT1 in hours
    TAI = float(UTC[0:2]) + float(UTC[3:5]) / 60 + float(UTC[6:8]) / 3600 + dAT / 3600  # TAI in hours
    GPS = TAI - 19 # GPS in hours
    TT = TAI + 32.184 # TT in hours

    pass

def J0(year, month, day):
    '''This function calculates the Julian date of a given date at 0h UT. 
    Between 1900 and 2100.
    Parameters
    ----------
    year : int
        Year
    month : int
        Month
    day : int
        Day
    Returns
    ----------
    J0 : float
        Julian date
    '''

    return 367 * year - np.floor(7 * (year + np.floor((month + 9) / 12)) / 4) + np.floor(275 * month / 9) + day + 1721013.5


def localSideralTime(year, month, day, ut, EL):
    '''This function calculates the local sideral time at a given date and time.
    Parameters
    ----------
    year : int
        Year
    month : int
        Month
    day : int
        Day
    ut : float
        Universal time in hours
    EL : float
        East longitude in degrees
    Returns
    ----------
    theta : float
        Local sideral time in degrees
    '''

    j0 = J0(year, month, day)
    j = (j0 - 2451545) / 36525

    g0 = 100.4606184 + 36000.77004 * j + 0.000387933 * j**2 - 2.583e-8 * j**3

    theta = g0 + 360.98564724 * ut + EL

    # Convert to the range [0, 360)
    g0 = zero2360(g0)

    lst = g0 + 360.98564724 * ut/24 + EL

    # Convert to the range [0, 360)
    lst = zero2360(lst)

    return lst

def zero2360(x):
    '''Converts an angle in degrees to the range [0, 360)'''

    if x >= 360:
        x = x - 360 * np.floor(x / 360)
    elif x < 0:
        x = x + 360 * np.ceil(np.abs(x) / 360)

    return x