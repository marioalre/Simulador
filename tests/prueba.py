from src.Orbit import Orbit
from src.CelestialBodies import CelestialBodies
import numpy as np
from src.rv2orb import Rv2orb

Tierra = CelestialBodies()
Tierra.earth()

# Orbit
Orbita = Orbit(Tierra)

r0 = np.array([ 1131.340, -2282.343, 6672.423])# km
v0 = np.array([ -5.64305, 4.30333,  2.42879])# km/s
dt = 40*60 # s
# Print results
r, v = Orbita.r0v02rv(r0 = r0, v0 = v0, dt = dt)

print('Position vector: {} km'.format(r))
print('Velocity vector: {} km/s'.format(v))

# Propagate the orbit

tf = 3600*24*360 # s
dt = 120 # s

r, v = Orbita.propagate(r0 = r0, v0 = v0, tf = tf, dt = dt)


t0 = 0
t = 3600

# Orbit

orbit = Rv2orb(r[-1, :], v[-1, :], Tierra, t0, t)

# Print results
print('Semi-major axis: {} km'.format(orbit.a))
print('Eccentricity: {}'.format(np.linalg.norm(orbit.e)))
print('Inclination: {} deg'.format(np.rad2deg(orbit.i)))
print('Ascending node: {} deg'.format(np.rad2deg(orbit.Omega)))
print('Argument of periapsis: {} deg'.format(np.rad2deg(orbit.omega)))
print('True anomaly: {} deg'.format(np.rad2deg(orbit.nu)))


