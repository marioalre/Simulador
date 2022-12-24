import numpy as np

# Import classes
from src.rv2orb import Rv2orb
from src.orb2rv import Orb2rv
from src.CelestialBodies import CelestialBodies
from src.Orbit import Orbit

if __name__ == '__main__':
    # Initial conditions
    r = np.array([-4039.9, 4814.56, 3628.62])# km
    v = np.array([-10.386, -4.77192, 1.74388])# km/s
    t0 = 0
    t = 3600

    cuerpo = CelestialBodies()
    cuerpo.earth()
    # Orbit

    orbit = Rv2orb(r, v, cuerpo, t0, t)

    # Print results
    print('Semi-major axis: {} km'.format(orbit.a))
    print('Eccentricity: {}'.format(np.linalg.norm(orbit.e)))
    print('Inclination: {} deg'.format(np.rad2deg(orbit.i)))
    print('Ascending node: {} deg'.format(np.rad2deg(orbit.Omega)))
    print('Argument of periapsis: {} deg'.format(np.rad2deg(orbit.omega)))
    print('True anomaly: {} deg'.format(np.rad2deg(orbit.nu)))


    print('--------------------------------------------------')
    h = 80000 # km
    e = 1.4 # unitless
    i = 30 * np.pi / 180 # rad
    Omega = 40 * np.pi / 180 # rad
    omega = 60 * np.pi / 180 # rad
    nu = 30 * np.pi / 180 # rad

    orbit = Orb2rv(h, e, Omega, omega, i, nu, cuerpo)

    print('Position vector: {} km'.format(orbit.r))
    print('Velocity vector: {} km/s'.format(orbit.v))



