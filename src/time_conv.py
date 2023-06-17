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

def ConvTime(date, UTC, dUT1, dAT, time_system='UT1'):
    '''The following formulas provide a means
    to determine the number of centuries elapsed
    from the epoch J2000.0 Ref: Vallado
    
    Parameters
    ----------
    date : list
        Date in the format [year, month, day]
    UTC : list
        UTC time in the format [hour, minute, second]
    dUT1 : float
        Difference between UT1 and UTC in seconds
    dAT : float 
        Difference between TAI and UTC in seconds
    time_system : string, optional
        Time system to be returned, by default 'UT1'
    Returns
    -------
    '''
    # UT1 = UTC + dUT1 (Universal Time 1)
    UT1 = UTC
    if UT1[2] + dUT1 >= 60:
        UT1[1] += 1
        UT1[2] = UT1[2] + dUT1 - 60
    else:
        UT1[2] += dUT1

    # TAI = UTC + dAT (International Atomic Time)
    TAI = UTC
    if TAI[2] + dAT >= 60:
        TAI[1] += 1
        TAI[2] = TAI[2] + dAT - 60
    else:
        TAI[2] += dAT

    # GPS = UTC + 19s  (Global Positioning System)
    GPS = UTC
    if GPS[2] + dAT - 19 >= 60:
        GPS[1] += 1
        GPS[2] = GPS[2] + dAT - 19 - 60
    else:
        GPS[2] += dAT - 19

    # TT = TAI + 32.184s (Terrestrial Time)
    TT = TAI
    if TT[2] + 32.184 >= 60:
        TT[1] += 1
        TT[2] = TT[2] + 32.184 - 60
    else:
        TT[2] += 32.184

    # Julian date TT
    JD_TT = J0(date[0], date[1], date[2], TT[0], TT[1], TT[2])
    T_TT = (JD_TT - 2451545) / 36525 # Julian centuries since J2000.0
    # Julian date TAI
    JD_TAI = J0(date[0], date[1], date[2], TAI[0], TAI[1], TAI[2])
    T_TAI = (JD_TAI - 2451545) / 36525 # Julian centuries since J2000.0

    if time_system == 'UT1':
        return UT1
    elif time_system == 'TAI':
        return TAI
    elif time_system == 'GPS':
        return GPS
    elif time_system == 'TT':
        return TT
    elif time_system == 'JD_TT':
        return JD_TT
    elif time_system == 'T_TT':
        return T_TT
    elif time_system == 'JD_TAI':
        return JD_TAI
    elif time_system == 'T_TAI':
        return T_TAI
    else:
        print('Time system not recognized')


def J0(year, month, day, h=0, m=0, s=0):
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
    h : float
        Hour
    m : float
        Minute
    s : float
        Second
    Returns
    ----------
    J0 : float
        Julian date
    '''

    return 367 * year - np.floor(7 * (year + np.floor((month + 9) / 12)) / 4) + np.floor(275 * month / 9) + day + 1721013.5 + (h + m / 60 + s / 3600) / 24

def MDJ(year, month, day, h, m, s):
    '''This function calculates the Modified Julian date of a given date at 0h UT'''
    return J0(year, month, day, h, m, s) - 2400000.5

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