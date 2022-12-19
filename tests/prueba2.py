import numpy as np
from src.Orbit import Orbit

if __name__ == '__main__':
    # Initial conditions
    r = np.array([2855, -35639, -24309])# km
    v = np.array([-1.656, -4.455, -3.2803])# km/s
    mu = 398600.4418 # km^3/s^2
    t0 = 0
    t = 3600

    # Orbit
    orbit = Orbit(r, v, mu, t0, t)

    # Print results
    print('Semi-major axis: {} km'.format(orbit.a))
    print('Eccentricity: {}'.format(orbit.e))
    print('Inclination: {} deg'.format(np.rad2deg(orbit.i)))
    print('Ascending node: {} deg'.format(np.rad2deg(orbit.Omega)))
    print('Argument of periapsis: {} deg'.format(np.rad2deg(orbit.omega)))
    print('True anomaly: {} deg'.format(np.rad2deg(orbit.nu)))
