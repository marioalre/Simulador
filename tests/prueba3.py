from src.CelestialBodies import CelestialBodies
from src.Lambert import Lambert
from src.Orbit import Orbit
from matplotlib import pyplot as plt
import numpy as np

# Initial position vector
r1 = np.array([5000, 10000, 2100])
r2 = np.array([-14600, 2500, 7000])

earth = CelestialBodies()
earth.earth()

lambert = Lambert(r1, r2, earth)
v1 = lambert.minimum_energy()

print('Initial velocity vector: {} km/s'.format(v1))

v11 , v22 = lambert.universal(3600, 'pro')
print('Initial velocity vector: {} km/s'.format(v11))
print('Final velocity vector: {} km/s'.format(v22))

# Propagador

orbit = Orbit(earth)

# Propagate the orbit
r, v, ax = orbit.propagate(r0 = r1, v0 = v1, tf = 3600*1.9, dt = 1)

# Añadir los vectores posicion r1 y r2
# Create a figure and an axes.
fig, ax = plt.subplots()
ax = fig.add_subplot(111, projection='3d')

# Plot the point in 3d.
ax.plot(r[:, 0], r[:, 1], r[:, 2], label = 'Orbit')

ax.plot([0], [0], [0], 'o', label='Earth', color='green', markersize=earth.radius/1000)

ax.plot([r1[0]], [r1[1]], [r1[2]], 'o', label='r1', color='red', markersize=5)

ax.plot([r2[0]], [r2[1]], [r2[2]], 'o', label='r2', color='blue', markersize=5)

ax.legend()

plt.show()


