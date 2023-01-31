import numpy as np
from src.Orbit import Orbit
from src.orb2rv import Orb2rv
from src.Lambert import Lambert

class Interplanetary:
    '''Interplanetary trajectory'''

    def __init__(self, central_body):
        self.mu = central_body.mu
        self.radius = central_body.radius
        self.central_body = central_body

    def sv_planets(self, planet, date, time):
        '''Compute position and velocity vectors of a planet at a given date and time
        Parameters
        ----------
        planet : CelestialBody
            Planet name
        date: array
            Date in format [year, month, day]
        time: array
            Time in format [hour, minute, second]
        Returns
        -------
        r : ndarray
            Position vector in km
        v : ndarray
            Velocity vector in km/s
        '''
        deg = np.pi / 180

        # Julian date
        j0 = J0(date[0], date[1], date[2])

        ut = (time[0] + time[1] / 60 + time[2] / 3600) / 24

        jd = j0 + ut

        J2000_coe, cent_rates_coe = self.planetary_elements(planet)

        t0 = (jd - 2451545) / 36525

        elements = J2000_coe + cent_rates_coe * t0

        a = elements[0]
        e = elements[1]

        h = np.sqrt(self.mu * a * (1 - e**2))

        i = elements[2]
        
        Omega = zero2360(elements[3])
        w_hat = zero2360(elements[4])
        L = zero2360(elements[5]) # mean longitude

        omega = zero2360(w_hat - Omega) # argument of latitude
        M = zero2360(L - w_hat)

        orb = Orbit(self.central_body)
        E = orb.kepler_elliptic(M*deg, e)

        TA = zero2360(2*np.arctan(np.sqrt((1 + e)/(1 - e))*np.tan(E/2))/deg)

        coe = np.array([h, e, Omega, i, omega, TA, a, w_hat, L, M, E/deg])

        orb = Orb2rv(a = a, e = e, i = i, Omega = Omega, M=M, omega=omega , body=self.central_body)

        r = orb.r
        v = orb.v

        return r, v

    def planetary_elements(self, planet):
        '''Planetary elements from JPL Horizons database
        Ref Orbital Mechanics for Engineering Students, 3rd ed, Howard Curtis, pg 388
        
        Parameters
        ----------
        Planet : CelestialBody
            Planet name
        Returns
        -------
        J2000_elements : numpy array
            9x6 array with the following elements:
            a, e, i, Omega, omega, L
            (L is the mean longitude)
        cent_rates : numpy array
            9x6 array with the following elements:
            da/dt, de/dt, di/dt, dOmega/dt, domega/dt, dL/dt
        '''

        J2000_elements = np.array([[0.38709893, 0.20563069, 7.00487, 48.33167, 77.45645, 252.25084],
                                  [0.72333199, 0.00677323, 3.39471, 76.68069, 131.53298, 181.97973],
                                  [1.00000011, 0.01671022, 0.00005, -11.26064, 102.94719, 100.46435],
                                  [1.52366231, 0.09341233, 1.85061, 49.57854, 336.04084, 355.45332],
                                  [5.20336301, 0.04839266, 1.30530, 100.55615, 14.75385, 34.40438],
                                  [9.53707032, 0.05415060, 2.48446, 113.71504, 92.43194, 49.94432],
                                  [19.19126393, 0.04716771, 0.76986, 74.22988, 170.96424, 313.23218],
                                  [30.06896348, 0.00858587, 1.76917, 131.72169, 44.97135, 304.88003],
                                  [39.48168677, 0.24880766, 17.14175, 110.30347, 224.06676, 238.92881]])

        cent_rates = np.array([[0.00000066, 0.00002527, -23.51, -446.30, 573.57, 538101628.29],
                              [0.00000092, -0.00004938, -2.86, -996.89, -108.80, 210664136.06],
                              [-0.00000005, -0.00003804, -46.94, -18228.25, 1198.28, 129597740.63],
                              [-0.00007221, 0.00011902, -25.47, -1020.19, 1560.78, 68905103.78],
                              [0.00060737, -0.00012880, -4.15, 1217.17, 839.93, 10925078.35],
                              [-0.00301530, -0.00036762, 6.11, -1591.05, -1948.89, 4401052.95],
                              [0.00152025, -0.00019150, -2.09, -1681.4, 1312.56, 1542547.79],
                              [-0.00125196, 0.00002514, -3.64, -151.25, -844.43, 786449.21],
                              [-0.00076912, 0.00006465, 11.07, -37.33, -132.25, 522747.90]])

        J2000_coe = J2000_elements[planet.planet_id-1]
        cent_rates_coe = cent_rates[planet.planet_id-1]

        # Convert from Au to km
        J2000_coe[0] = J2000_coe[0] * 149597870.700
        cent_rates_coe[0] = cent_rates_coe[0] * 149597870.700

        # Convert from arcsec to fraction of a degree

        cent_rates_coe[2] = cent_rates_coe[2] / 3600
        cent_rates_coe[3] = cent_rates_coe[3] / 3600
        cent_rates_coe[4] = cent_rates_coe[4] / 3600
        cent_rates_coe[5] = cent_rates_coe[5] / 3600

        return J2000_coe, cent_rates_coe
        
    def planetary_mission(self, planet_departure, date_departure, planet_arr, date_arr):
        '''Interplanetary mission by patched conics method
        Parameters
        ----------
        planet_departure : CelestialBody
            Planet of departure
        date_departure : array
            Departure date
            Format: [year, month, day, hour, minute, second]
        planet_arr : CelestialBody
            Planet of arrival
        date_arr : array
            Arrival date
            Format: [year, month, day, hour, minute, second]
        Returns
        -------
        r : numpy array
            Position vector
        v : numpy array
            Velocity vector
        '''

        # Departure date
        r_departure, v_departure = self.sv_planets(planet_departure, date_departure[:2], date_departure[3:5])

        r_arr , v_arr = self.sv_planets(planet_arr, date_arr[:2], date_arr[3:5])

        # Calculate the delta v with the Lambert's problem
        deltaV = Lambert(r_departure, r_arr, self.central_body)

        # Calculate the total delta v
        # Optimized time of flight to minimize the total delta v with a genetic algorithm

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



if __name__ == '__main__':
    print('Testing the function in the module')
    
    # Create a planet object
    from src.CelestialBodies import CelestialBodies

    planet = CelestialBodies()
    Sun = CelestialBodies()
    Sun.sun()
    planet.earth()

    date = np.array([2003, 8, 27])
    time = np.array([12, 0, 0])

    print(type(Sun))

    # Calculate the position of the planet
    inter = Interplanetary(central_body=Sun)
    pos = inter.sv_planets(planet, date, time)

    print(f'ECI frame: position = {pos[0]} km, velocity = {pos[1]} km/s')




