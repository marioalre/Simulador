import numpy as np
import folium
import math
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

def ecef2latlongh(r_ecef):
    '''
    ECEF to latitude, longitude and height
    Algorithm 12 from Vallado

    Verified with the example 3.3 from Vallado
    '''

    e = 0.081819221456
    Req = 6378.1363

    r_delta_sat = np.sqrt(r_ecef[0]**2 + r_ecef[1]**2)

    # Check quadrant

    alpha = np.arcsin(r_ecef[1] / r_delta_sat)

    lambda_ = alpha # first approximation

    delta = np.arctan(r_ecef[2] / r_delta_sat)

    # CI
    long_gd = delta # geodetic longitude
    r_delta = r_delta_sat
    r_k = r_ecef[2]

    long_old = 0
    if np.abs(long_gd - long_old) < 1e-10:
        long_old = -1


    while np.abs(long_gd - long_old) > 1e-10:

        long_old = long_gd
        C = Req / np.sqrt(1 - e**2 * np.sin(long_old)**2)
        S = Req* (1 - e**2) / np.sqrt(1 - e**2 * np.sin(long_old)**2)

        # Newton-Raphson method to solve the equation
        # long = np.arctan((r_k + C * e**2 * np.sin(long)) / r_delta)

        long_gd = long_old - (np.tan(long_old) - (C*np.sin(long_old)*e**2+r_k)/r_delta) \
            / (np.tan(long_old)**2 - (C*e**2*np.cos(long_old))/r_delta + 1)
        
    long_gc = np.arctan((1-e**2)*np.tan(long_gd))



    if -1e-10 < long_gd < 1e-10:
        long_gd = 1e-8
        long_gc = np.arctan((1-e**2)*np.tan(long_gd))

    if long_gd < np.pi/180:
        hell = r_k / np.sin(long_gd) - S
    else:
        hell = r_delta / np.cos(long_gd) - C

    print('Altitud: ', hell)
    print('Latitud: ', lambda_ * 180 / np.pi)
    print('Longitud geodésica: ', long_gd * 180 / np.pi)
    print('Longitud geocéntrica: ', long_gc * 180 / np.pi)

    return long_gd, long_gc, lambda_, hell

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

def fk4():
    # Choose one of the three methods to define the conversion matrix
    # Method 1: in book
    # 1950 - 2000
    fk4m = np.array([
        [0.9999256794956877, -0.0111814832204662, -0.0048590038153592],
        [0.0111814832391717,  0.9999374848933135, -0.0000271625947142],
        [0.0048590037723143, -0.0000271702937440,  0.9999881946043742]
    ])

    # Method 2: stk approach
    # New way is formed by multiplying the matrices on pages
    # 173 and 174 and adding in the correction to equinox given
    # on page 168 of the supplement to the astronomical almanac 
    # 1950 - 2000
    fk4m = np.array([
        [0.999925678612394, -0.011181874556714, -0.004858284812600],
        [0.011181874524964,  0.999937480517880, -0.000027169816135],
        [0.004858284884778, -0.000027156932874,  0.999988198095508]
    ])

    # Method 3: from Exp supp to Ast Almanac pg 185 6x6
    # 1950 - 2000
    '''
    fk4m = np.array([
        [0.9999256782, -0.0111820611, -0.0048579477,  0.00000242395018, -0.00000002710663, -0.00000001177656],
        [0.0111820610,  0.9999374784, -0.0000271765,  0.00000002710663,  0.00000242397878, -0.00000000006587],
        [0.0048579479, -0.0000271474,  0.9999881997,  0.00000001177656, -0.00000000006582,  0.00000242410173],
        [-0.000551     , -0.238565     ,  0.435739     ,  0.99994704     , -0.01118251     , -0.00485767     ],
        [ 0.238514     , -0.002667     , -0.008541     ,  0.01118251     ,  0.99995883     , -0.00002718     ],
        [-0.435623     ,  0.012254     ,  0.002117     ,  0.00485767     , -0.00002714     ,   1.00000956    ]
    ])'''

    return fk4m


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

# This class allows to convert ECI coordinates to ECEF coordinates following the algorithm described in the
# Fundamentals of Astrodynamics and Applications book by David A. Vallado.
class ECI2ECEF:
    ''' This class allows to convert ECI coordinates to ECEF coordinates following the algorithm described in the
        Fundamentals of Astrodynamics and Applications book by David A. Vallado.
    '''

    # Define the constructor that takes the input parameters
    def __init__(self, ttt, jdut1, lod, xp, yp, eqeterms, ddpsi, ddeps):
        ''' 
        This function initializes the class with the input parameters.
        Inputs:
            reci: ECI position vector (km)
            veci: ECI velocity vector (km/s)
            ttt: Julian centuries of terrestrial time
            jdut1: Julian date of UT1
            lod: Length of day (sec)
            xp: Polar motion coefficient (arcsec)
            yp: Polar motion coefficient (arcsec)
            eqeterms: Boolean to select if the extra terms for nutation are used (0 or 2)
            ddpsi: Correction to delta psi (arcsec)
            ddeps: Correction to delta eps (arcsec)
        '''
        self.ttt = ttt
        self.jdut1 = jdut1
        self.lod = lod
        self.xp = xp
        self.yp = yp
        self.eqeterms = eqeterms
        self.ddpsi = ddpsi
        self.ddeps = ddeps

    # Define the function that converts ECI coordinates to ECEF
    def eci2ecef(self, reci, veci, aeci):
        self.reci = reci
        self.veci = veci
        self.aeci = aeci

        # Call the auxiliary functions to get the precession, nutation and rotation parameters
        prec, psia, wa, ea, xa = self.precess(self.ttt)
        deltapsi, trueeps, meaneps, omega, nut = self.nutation(self.ttt, self.ddpsi, self.ddeps)
        st, stdot = self.sideral(self.jdut1, deltapsi, meaneps, omega, self.lod, self.eqeterms)
        pm = self.polarm(self.xp, self.yp, self.ttt)

        # Calculate the angular velocity of the Earth
        thetasa = 7.29211514670698e-05 * (1.0 - self.lod / 86400.0)
        omegaearth = np.array([0, 0, thetasa])

        # Convert ECI coordinates to ECEF using transformation matrices
        rpef = st.T @ nut.T @ prec.T @ self.reci
        recef = pm.T @ rpef

        vpef = st.T @ nut.T @ prec.T @ self.veci - np.cross(omegaearth, rpef)
        vecef = pm.T @ vpef

        temp = np.cross(omegaearth, rpef)

        # Two additional terms that are not necessary if the satellite is not on the surface of the Earth
        aecef = pm.T @ (st.T @ nut.T @ prec.T @ self.aeci) - np.cross(omegaearth,temp) - 2.0 * np.cross(omegaearth,vpef)

        # Return the results
        return recef, vecef, aecef
    
     # Define the function that converts ECI coordinates to ECEF
    def ecef2eci(self, recef, vecef, aecef):
        self.recef = recef
        self.vecef = vecef
        self.aecef = aecef

        # Call the auxiliary functions to get the precession, nutation and rotation parameters
        prec, psia, wa, ea, xa = self.precess(self.ttt)
        deltapsi, trueeps, meaneps, omega, nut = self.nutation(self.ttt, self.ddpsi, self.ddeps)
        st, stdot = self.sideral(self.jdut1, deltapsi, meaneps, omega, self.lod, self.eqeterms)
        pm = self.polarm(self.xp, self.yp, self.ttt)

        # Calculate the angular velocity of the Earth
        thetasa = 7.29211514670698e-05 * (1.0 - self.lod / 86400.0)
        omegaearth = np.array([0, 0, thetasa])

        # Convert ECI coordinates to ECEF using transformation matrices
        rpef = pm @ recef
        reci = prec @ nut @ st @ self.recef

        vpef = pm @ vecef
        veci = prec @ nut @ st @ (vpef + np.cross(omegaearth, rpef))

        temp = np.cross(omegaearth, rpef)

        # Two additional terms that are not necessary if the satellite is not on the surface of the Earth
        aeci = prec @ nut @ st @ (pm @ aecef) + np.cross(omegaearth,temp) + 2.0 * np.cross(omegaearth,vpef) 

        # Return the results
        return reci, veci, aeci

    # Define the auxiliary functions for precession, nutation and rotation

    def get_ttt(self, JD_tt):
        '''This function calculates the Julian centuries of terrestrial time
        Parameters
        ----------
        JD_tt : float
            Julian date of TT
        Returns
        ----------
        ttt : float
            Julian centuries of terrestrial time
        '''

        ttt = (JD_tt - 2451545.0) / 36525.0


    def precess(self, ttt):
        '''Based on Vallado's book, Fundamentals of Astrodynamics and Applications, 4th ed., Algorithm 23, p. 211
        and the Matlab code aviable at celestrak.com

        VALIDATED

        Parameters
        ----------
        ttt : float
            Julian centuries of terrestrial time
        opt : str
            '50' for fk4 b1950 precession angles
            '80' for iau 76 precession angles
            '00' for iau 06 precession angles
        Returns
        ----------
        prec : array
            Precession matrix
        '''
        # Code for precess function goes here
        convrt = np.pi / (180.0*3600.0)
        ttt2 = ttt * ttt
        ttt3 = ttt2 * ttt

        prec = np.zeros((3, 3))

        # iau 76 precession angles

        psia = 5038.7784 * ttt - 1.07259 * ttt2 - 0.001147 * ttt3
        wa = 84381.448 + 0.05127 * ttt2 - 0.007726 * ttt3
        ea = 84381.448 - 46.8150 * ttt - 0.00059 * ttt2 + 0.001813 * ttt3
        xa = 10.5526 * ttt - 2.38064 * ttt2 - 0.001125 * ttt3

        zeta = 2306.2181 * ttt + 0.30188 * ttt2 + 0.017998 * ttt3
        theta = 2004.3109 * ttt - 0.42665 * ttt2 - 0.041833 * ttt3
        z = 2306.2181 * ttt + 1.09468 * ttt2 + 0.018203 * ttt3

        # convert from arcseconds to radians
        psia *= convrt
        wa *= convrt
        ea *= convrt
        xa *= convrt
        zeta *= convrt
        theta *= convrt
        z *= convrt

        # calculate the precession matrix elements
        cz = np.cos(z)
        sz = np.sin(z)
        coszeta = math.cos(zeta)
        sinzeta = math.sin(zeta)
        costheta = math.cos(theta)
        sintheta = math.sin(theta)
        cosz = math.cos(z)
        sinz = math.sin(z)

        # Formar la matriz mod to j2000
        prec = np.zeros((3, 3))
        prec[0, 0] = coszeta * costheta * cosz - sinzeta * sinz
        prec[0, 1] = coszeta * costheta * sinz + sinzeta * cosz
        prec[0, 2] = coszeta * sintheta
        prec[1, 0] = -sinzeta * costheta * cosz - coszeta * sinz
        prec[1, 1] = -sinzeta * costheta * sinz + coszeta * cosz
        prec[1, 2] = -sinzeta * sintheta
        prec[2, 0] = -sintheta * cosz
        prec[2, 1] = -sintheta * sinz
        prec[2, 2] = costheta


        return prec, psia, wa, ea, xa


    def nutation(self, ttt, ddpsi, ddeps):
        '''Based on Vallado's book, Fundamentals of Astrodynamics and Applications, 4th ed., Algorithm 22, p. 205
        and the Matlab code aviable at celestrak.com

        VALIDATED
        '''

        deg2rad = math.pi / 180.0

        iar80, rar80 = self.iau80in()  # coeff in deg

        # ---- determine coefficients for iau 1980 nutation theory ----
        ttt2 = ttt * ttt
        ttt3 = ttt2 * ttt

        meaneps = -46.8150 * ttt - 0.00059 * ttt2 + 0.001813 * ttt3 + 84381.448
        meaneps = (meaneps / 3600.0) % 360.0
        meaneps = meaneps * deg2rad

        l, l1, f, d, omega, lonmer, lonven, lonear, lonmar, lonjup, lonsat, lonurn, lonnep, precrate = self.fundarg(ttt)

        deltapsi = 0.0
        deltaeps = 0.0
        for i in range(105, -1, -1):
            tempval = iar80[i - 1, 0] * l + iar80[i - 1, 1] * l1 + iar80[i - 1, 2] * f + iar80[i - 1, 3] * d + iar80[i - 1, 4] * omega
            deltapsi = deltapsi + (rar80[i - 1, 0] + rar80[i - 1, 1] * ttt) * math.sin(tempval)
            deltaeps = deltaeps + (rar80[i - 1, 2] + rar80[i - 1, 3] * ttt) * math.cos(tempval)

        # --------------- find nutation parameters --------------------
        a1 = deltapsi + ddpsi
        a2 = deltaeps + ddeps

        if a1 < 0.0:
            deltapsi = np.remainder(-a1,  2.0 * math.pi)
            deltapsi = -deltapsi
        else:
            deltapsi = np.remainder(a1,  2.0 * math.pi)

        if a2 < 0.0:
            deltaeps = np.remainder(-a2, 2.0 * math.pi)
            deltaeps = -deltaeps
        else:
            deltaeps = np.remainder(a2, 2.0 * math.pi)

        trueeps = meaneps + deltaeps

        cospsi = math.cos(deltapsi)
        sinpsi = math.sin(deltapsi)
        coseps = math.cos(meaneps)
        sineps = math.sin(meaneps)
        costrueeps = math.cos(trueeps)
        sintrueeps = math.sin(trueeps)

        nut = np.zeros((3, 3))
        nut[0, 0] = cospsi
        nut[0, 1] = costrueeps * sinpsi
        nut[0, 2] = sintrueeps * sinpsi
        nut[1, 0] = -coseps * sinpsi
        nut[1, 1] = costrueeps * coseps * cospsi + sintrueeps * sineps
        nut[1, 2] = sintrueeps * coseps * cospsi - sineps * costrueeps
        nut[2, 0] = -sineps * sinpsi
        nut[2, 1] = costrueeps * sineps * cospsi - sintrueeps * coseps
        nut[2, 2] = sintrueeps * sineps * cospsi + costrueeps * coseps

        return deltapsi, trueeps, meaneps, omega, nut
    
    def iau80in(self):
        # 0.0001" to rad
        convrt = 0.0001 * np.pi / (180 * 3600.0)

        nut80 = np.loadtxt("./data/nut80.dat")

        iar80 = nut80[:, 0:5]
        rar80 = nut80[:, 5:9]

        for i in range(106):
            for j in range(4):
                rar80[i, j] = rar80[i, j] * convrt

        return iar80, rar80
    

    def fundarg(self, ttt):
        deg2rad = np.pi / 180.0

        # ---- determine coefficients for iau 2000 nutation theory ----
        ttt2 = ttt * ttt
        ttt3 = ttt2 * ttt
        ttt4 = ttt2 * ttt2

        l = ((((0.064) * ttt + 31.310) * ttt + 1717915922.6330) * ttt) / 3600.0 + 134.96298139
        l1 = ((((-0.012) * ttt - 0.577) * ttt + 129596581.2240) * ttt) / 3600.0 + 357.52772333
        f = ((((0.011) * ttt - 13.257) * ttt + 1739527263.1370) * ttt) / 3600.0 + 93.27191028
        d = ((((0.019) * ttt - 6.891) * ttt + 1602961601.3280) * ttt) / 3600.0 + 297.85036306
        omega = ((((0.008) * ttt + 7.455) * ttt - 6962890.5390) * ttt) / 3600.0 + 125.04452222
        lonmer = 252.3 + 149472.0 * ttt
        lonven = 179.9 + 58517.8 * ttt
        lonear = 98.4 + 35999.4 * ttt
        lonmar = 353.3 + 19140.3 * ttt
        lonjup = 32.3 + 3034.9 * ttt
        lonsat = 48.0 + 1222.1 * ttt
        lonurn  =   0.0
        lonnep  =   0.0
        precrate=   0.0

        l = np.remainder(l, 360.0) * deg2rad
        l1 = np.remainder(l1, 360.0) * deg2rad
        f = np.remainder(f, 360.0) * deg2rad
        d = np.remainder(d, 360.0) * deg2rad
        omega = np.remainder(omega, 360.0) * deg2rad
        lonmer = np.remainder(lonmer, 360.0) * deg2rad
        lonven = np.remainder(lonven, 360.0) * deg2rad
        lonear = np.remainder(lonear, 360.0) * deg2rad
        lonmar = np.remainder(lonmar, 360.0) * deg2rad
        lonjup = np.remainder(lonjup, 360.0) * deg2rad
        lonsat = np.remainder(lonsat, 360.0) * deg2rad
        lonurn = np.remainder(lonurn, 360.0) * deg2rad
        lonnep = np.remainder(lonnep, 360.0) * deg2rad
        precrate = np.remainder(precrate, 360.0) * deg2rad

        return l, l1, f, d, omega, lonmer, lonven, lonear, lonmar, lonjup, lonsat, lonurn, lonnep, precrate


    def sideral(self, jdut1, deltapsi, meaneps, omega, lod, eqeterms):
        '''
        VALIDATED
        '''
        gmst = self.gstime(jdut1)

        if (jdut1 > 2450449.5) and (eqeterms > 0):
            ast = gmst + deltapsi * math.cos(meaneps) + 0.00264 * math.pi / (3600 * 180) * math.sin(omega) + 0.000063 * math.pi / (3600 * 180) * math.sin(2.0 * omega)
        else:
            ast = gmst + deltapsi * math.cos(meaneps)

        ast = ast % (2.0 * math.pi)
        thetasa = 7.29211514670698e-05 * (1.0 - lod / 86400.0)
        omegaearth = thetasa

        st = np.array([[math.cos(ast), -math.sin(ast), 0.0],
            [math.sin(ast), math.cos(ast), 0.0],
            [0.0, 0.0, 1.0]])

        stdot = np.array([[-omegaearth * math.sin(ast), -omegaearth * math.cos(ast), 0.0],
                [omegaearth * math.cos(ast), -omegaearth * math.sin(ast), 0.0],
                [0.0, 0.0, 0.0]])

        return st, stdot
    

    def gstime(self, jdut1):

        twopi = 2.0 * math.pi
        deg2rad = math.pi / 180.0

        tut1 = (jdut1 - 2451545.0) / 36525.0

        gst = -6.2e-6 * tut1**3 + 0.093104 * tut1**2 + (876600.0 * 3600.0 + 8640184.812866) * tut1 + 67310.54841

        gst = (gst * deg2rad / 240.0) % (twopi)

        if gst < 0.0:
            gst += twopi

        return gst
    

    def polarm(self, xp, yp, ttt):
        '''
        VALIDATED
        '''
        cosxp = math.cos(xp)
        sinxp = math.sin(xp)
        cosyp = math.cos(yp)
        sinyp = math.sin(yp)

        pm = np.zeros((3, 3))

        pm[0, 0] = cosxp
        pm[0, 1] = 0.0
        pm[0, 2] = -sinxp
        pm[1, 0] = sinxp * sinyp
        pm[1, 1] = cosyp
        pm[1, 2] = cosxp * sinyp
        pm[2, 0] = sinxp * cosyp
        pm[2, 1] = -sinyp
        pm[2, 2] = cosxp * cosyp
        
        return pm


if __name__ == "__main__":
    
    val = ecef2latlongh([6524.834, 6862.875, 6448.296])
    # Example usage
    # Create a j200 position vector in km
    rj200 = np.array([10000.,20000.,30000])

    # Convert to b195 position vector using method of choice
    rb195 = np.dot(fk4(), rj200)

    # Print the result
    print(rb195)

    # Ejemplo vallado 3.15

    # seguntos to arco segundos
    
    convrt = np.pi / (180.0 * 3600.0)

    r_ecef = [-1033.4793830,  7901.2952754,  6380.3565958]
    v_ecef = [-3.225636520,  -2.872451450,   5.531924446]
    a_ecef = [0.001,0.002,0.003]
    
    print(f"r_ecef: {r_ecef}")
    print(f"v_ecef: {v_ecef}")
    print(f"a_ecef: {a_ecef}")

    ttt = 0.0426236319
    jdut1 = 2453101.827406783
    lod = 0.001556300
    xp = -6.820455828585174e-07
    yp = 1.615927632369383e-06
    eqeterms = 2
    ddpsi = -2.530485008551223e-7
    ddeps = -1.878653014299452e-8
    conv = ECI2ECEF(ttt, jdut1, lod, xp, yp, eqeterms, ddpsi, ddeps)

    reci, veci, aeci = conv.ecef2eci(r_ecef, v_ecef, a_ecef)
    print(f"reci: {reci}")
    print(f"veci: {veci}")
    print(f"aeci: {aeci}")
    recef, vecef, aecef = conv.eci2ecef(reci, veci, aeci)
    print(f"recef: {recef}")
    print(f"vecef: {vecef}")
    print(f"aecef: {aecef}")

    