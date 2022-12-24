from src.Orbit import Orbit
from src.CelestialBodies import CelestialBodies
import numpy as np

Tierra = CelestialBodies()
Tierra.earth()

# Orbit
Orbita = Orbit(Tierra)

r0 = np.array([ 7000, -12124, 0])# km
v0 = np.array([2.6679, 4.6210, 0])# km/s
dt = 3600 # s
# Print results
r, v = Orbita.r0v02rv(r0 = r0, v0 = v0, dt = dt)

print('Position vector: {} km'.format(r))
print('Velocity vector: {} km/s'.format(v))


