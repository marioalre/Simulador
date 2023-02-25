import numpy as np
import folium
from src.time_conv import ConvTime

def rotation_matrix_1(angle):
    '''Rotation matrix around the x axis'''
    return np.array([
        [1, 0, 0],
        [0, np.cos(angle), np.sin(angle)],
        [0, -np.sin(angle), np.cos(angle)]
        ])

def rotation_matrix_2(angle):
    '''Rotation matrix around the y axis'''
    return np.array([
        [np.cos(angle), 0, -np.sin(angle)],
        [0, 1, 0],
        [np.sin(angle), 0, np.cos(angle)]
    ])

def rotation_matrix_3(angle):
    '''Rotation matrix around the z axis'''
    return np.array([
        [np.cos(angle), np.sin(angle), 0],
        [-np.sin(angle), np.cos(angle), 0],
        [0, 0, 1]
    ])

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

def latlongH2ecef(phi, lam, h):
    ''' Latitude, longitude and height to ECEF
    Parameters
    ----------
    phi : float
        Latitude in degrees
    lam : float
        Longitude in degrees
    h : float   
        Ellipsoidal height in m
    Returns
    -------
    r : array_like  
        Position vector in km ECEF
    '''
    # Earth Constants
    a = 6378137.0  # equatorial radius (meters)
    f = 1/298.257223563  # flattening of the ellipsoid
    b = a * (1 - f)  # polar radius (meters)
    e2 = 1 - (b/a)**2  # eccentricity squared

    # Convert to radians
    phi_rad = phi * np.pi / 180
    lam_rad = lam * np.pi / 180

    # Calculate the radius of curvature in the prime vertical N
    N = a / np.sqrt(1 - e2 * np.sin(phi_rad)**2)

    # Calculate the cartesian coordinates (X,Y,Z) in the ECEF system
    X = (N + h) * np.cos(phi_rad) * np.cos(lam_rad)
    Y = (N + h) * np.cos(phi_rad) * np.sin(lam_rad)
    Z = ((1 - e2) * N + h) * np.sin(phi_rad)

    # Vector of position in ECEF coordinates
    return [X, Y, Z]


def plot_ground_track(lat = None, long = None, map = None):
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
    if map is None:
        # Crea un mapa centrado en las coordenadas (25, 35)
        map = folium.Map(location=[0, 0], zoom_start=2.4, min_zoom=2.4 , tiles='Stamen Terrain')

    # Dibuja la trayectoria del satélite en el mapa
    '''
    folium.PolyLine(
        locations=list(zip(lat, long)),
        color='red',
        weight=2,
        opacity=1
    ).add_to(m)'''

    colors = ['red', 'blue', 'green', 'purple', 'orange', 'darkred', 'lightred', 'beige', 'darkblue', 'darkgreen', 'cadetblue', 'darkpurple', 'white', 'pink', 'lightblue', 'lightgreen', 'gray', 'black', 'lightgray']
    # Random colors
    color = np.random.choice(colors)
    for latlong in zip(lat, long):
        folium.CircleMarker(location=latlong, radius=0.6, color=color).add_to(map)

    return map

def fk5(self, r_gcrf, v_gcrf, date, UTC, dUT1, dAT, xp, yp):
        '''This function is used to convert the position and velocity vectors from GCRF to FK5
        celestial reference frame (GCRF) to terrestrial reference frame (ITRF).
        Parameters
        ----------
        r_gcrf : numpy array
            Position vector in GCRF
        v_gcrf : numpy array
            Velocity vector in GCRF
        date : list
            Date in the format [year, month, day]
        UTC : string
            UTC time in the format HH:MM:SS
        dUT1 : float
            Difference between UT1 and UTC
        dAT : float
            Difference between TAI and UTC
        xp : float
            Polar motion x-coordinate
        yp : float
            Polar motion y-coordinate
        Returns
        ------- 
        r_ITRF : numpy array
            Position vector in ITRF
        v_ITRF : numpy array
            Velocity vector in ITRF
        '''
        UT1 = ConvTime(date, UTC, dUT1, dAT, 'UT1')
        TAI = ConvTime(date, UTC, dUT1, dAT, 'TAI')
        TT = ConvTime(date, UTC, dUT1, dAT, 'TT')
        T_UT1 = ConvTime(date, UTC, dUT1, dAT, 'T_UT1')
        T_TT = ConvTime(date, UTC, dUT1, dAT, 'T_TT')

        # Find precession and nutation angles
        # Precession angles
        pass


def eci_to_ecef(r_ECI, v_ECI, t, theta):
    """
    Converts a position and velocity vector in ECI coordinates (r_ECI, v_ECI) 
    to a position and velocity vector in ECEF coordinates (r_ECEF, v_ECEF).
    
    Args:
        r_ECI (numpy.ndarray): Position vector in ECI coordinates [X, Y, Z] (meters).
        v_ECI (numpy.ndarray): Velocity vector in ECI coordinates [Vx, Vy, Vz] (meters/second).
        t (float): Time at which the transformation is desired (seconds).
        theta (float): Greenwich sidereal time corresponding to time t (radians).
    
    Returns:
        tuple: A tuple containing the position vector in ECEF coordinates [X, Y, Z] (meters) 
        and the velocity vector in ECEF coordinates [Vx, Vy, Vz] (meters/second).
    """
    # Define the Earth's rotation matrix around the z-axis (theta)
    Rz = np.array([
        [np.cos(theta), np.sin(theta), 0],
        [-np.sin(theta), np.cos(theta), 0],
        [0, 0, 1]
    ])
    
    # Calculate Greenwich Mean Time (GMT) from the sidereal time
    GMT = theta/np.pi/12*t
    
    # Define the Earth's rotation matrix around the x-axis (incl. nutation)
    delta_eps = 0.409105176667471 # mean obliquity of the ecliptic
    delta_psi = -0.000226965 * np.sin(np.radians(125.04 - 1934.136*GMT)) # nutation in longitude
    delta_eps = 0.00256 * np.cos(np.radians(125.04 - 1934.136*GMT)) # nutation in obliquity
    Rx = np.array([
        [1, 0, 0],
        [0, np.cos(delta_eps), np.sin(delta_eps)],
        [0, -np.sin(delta_eps), np.cos(delta_eps)]
    ]) @ np.array([
        [1, 0, 0],
        [0, np.cos(delta_psi), np.sin(delta_psi)],
        [0, -np.sin(delta_psi), np.cos(delta_psi)]
    ])
    
    # Total rotation matrix
    R = Rx @ Rz
    
    # Convert ECI position and velocity vectors to ECEF
    r_ECEF = R @ r_ECI
    v_ECEF = R @ v_ECI + np.cross(np.array([0, 0, 7.292115146706979e-5]), r_ECEF)
    
    return r_ECEF, v_ECEF



    