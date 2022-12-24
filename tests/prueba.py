from src.Orbit import Orbit
from src.CelestialBodies import CelestialBodies
import numpy as np

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


