import numpy as np

def rotation_matrix_1(angle):
    '''Rotation matrix around the x axis'''
    return np.array([[1, 0, 0], [0, np.cos(angle), -np.sin(angle)], [0, np.sin(angle), np.cos(angle)]])

def rotation_matrix_2(angle):
    '''Rotation matrix around the y axis'''
    return np.array([[np.cos(angle), 0, np.sin(angle)], [0, 1, 0], [-np.sin(angle), 0, np.cos(angle)]])

def rotation_matrix_3(angle):
    '''Rotation matrix around the z axis'''
    return np.array([[np.cos(angle), -np.sin(angle), 0], [np.sin(angle), np.cos(angle), 0], [0, 0, 1]])

def cart2efix(r, v, t, dt, t0):
    ''' Space fixed coordinate system to Earth fixed coordinate system
    Parameters
    ----------
    r : array_like
        Position vector in km
    v : array_like
        Velocity vector in km/s
    t : float
        Time in seconds after t0
    dt : float
        Time step in seconds
    t0 : float
        Time of reference in seconds
    Returns
    -------
    r : array_like
        Position vector in km
    v : array_like
        Velocity vector in km/s
    '''

    # Earth's rotation rate
    omega = 7.2921158553e-5

    hour = t0 / 3600
    sid = hour * 15 * np.pi / 180  

    t = np.arange(0, t+dt, dt)

    theta0 = omega * t + sid

    for i in range(len(theta0)):
        R = np.transpose(rotation_matrix_3(theta0[i]))
        r[i, :] = np.dot(R, r[i, :])
        v[i, :] = np.dot(R, v[i, :])

    return r, v


def efix2cart():
    ''' Earth fixed coordinate system to Space fixed coordinate system'''

def sid2rad(second):
    ''' Sidereal time to radians'''

    hour = second / 3600

    return hour * 15 * np.pi / 180

def rad2sid(radians):
    ''' Radians to sidereal time'''

    hour = radians * 180 / np.pi / 15

    return hour