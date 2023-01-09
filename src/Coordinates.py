import numpy as np
import folium
from IPython.display import display



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
    pass

def sid2rad(second):
    ''' Sidereal time to radians'''

    hour = second / 3600

    return hour * 15 * np.pi / 180

def rad2sid(radians):
    ''' Radians to sidereal time'''

    hour = radians * 180 / np.pi / 15

    return hour

def ecef2latlong(r):
    ''' ECEF to latitude and longitude
    Algorithm 12 from Vallado

    Parameters
    ----------
    r : array_like
        Position vector in km ECEF
    Returns
    -------
    lat : float
        Latitude in radians
    lon : float 
        Longitude in radians
    '''

    lon = np.arctan2(r[:,1], r[:, 0])
    lon  *= 180 / np.pi 
    lat = np.arctan2(r[:, 2], np.sqrt(r[:, 0]**2 + r[:, 1]**2))
    lat *= 180 / np.pi

    # Check quadrant
    # Between -180 and 180
    
    for i in range(len(lon)):
        if lon[i] < -180:
            lon[i] += 360
        elif lon[i] > 180:
            lon[i] -= 360

        if lat[i] < -90:
            lat[i] += 180
        elif lat[i] > 90:
            lat[i] -= 180

    return lat, lon

def plot_ground_track(lat = None, long = None):
    ''' Plot ground track
    Parameters
    ----------
    lat : array_like
        Latitude in degrees
    long : array_like
        Longitude in degrees
    Returns
    -------
    map : folium map
        Folium map object
    '''

    # Crea un mapa centrado en las coordenadas (25, 35)
    m = folium.Map(location=[0, 0], zoom_start=2, tiles='Stamen Terrain')

    # Dibuja la trayectoria del satélite en el mapa
    folium.PolyLine(
        locations=list(zip(lat, long)),
        color='red',
        weight=2,
        opacity=1
    ).add_to(m)

    # Guarda el mapa en un archivo HTML
    m.save('results/satellite_track.html')
    display(m)

    return m


    