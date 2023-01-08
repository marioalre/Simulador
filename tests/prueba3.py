from src.CelestialBodies import CelestialBodies
from src.Lambert import Lambert
from src.Orbit import Orbit
from src.orb2rv import Orb2rv
from matplotlib import pyplot as plt
import numpy as np

earth = CelestialBodies()
earth.earth()

# Molnya parameters

orb = Orb2rv(a = 26554, e = 0.72, Omega=0, i = 63.4, omega = -90, nu = 0, body=earth)
r0 = orb.position_eci()
v0 = orb.velocity_eci()

# Geo
orb = Orb2rv(a = 42165, e = 0, Omega=0, i = 0, omega = 0, nu = 0, body=earth)
r01 = orb.position_eci()
v01 = orb.velocity_eci()


lambert = Lambert(r0, r01, earth)
v1 = lambert.minimum_energy()

print('Initial velocity vector: {} km/s'.format(v1))

v11 , v22 = lambert.universal(3600, 'pro')
print('Initial velocity vector: {} km/s'.format(v11))
print('Final velocity vector: {} km/s'.format(v22))

# Propagador

orbit = Orbit(earth)

# Propagate the orbit
r, v = orbit.propagate(r0 = r0, v0 = v1, tf = 3600*4.6, dt = 1)


# Propagate the orbit

tf = 3600*24 # s
dt = 360 # s

r1, v1 = orbit.propagate(r0 = r0, v0 = v0, tf = tf, dt = dt)
r2, v2 = orbit.propagate(r0 = r01, v0 = v01, tf = tf, dt = dt)


# Añadir los vectores posicion r1 y r2
# Create a figure and an axes.
plt.style.use('classic')
fig = plt.figure()

ax = fig.add_subplot(111, projection='3d')

# Plot the point in 3d.
ax.plot(r[:, 0]/1000, r[:, 1]/1000, r[:, 2]/1000, label = 'Tansferrencia', color='orange', linewidth=1)

ax.plot([0], [0], [0], 'o', label='Tierra', color='green', markersize=earth.radius/1000)

ax.plot([r1[0, 0]/1000], [r1[0,1]/1000], [r1[0,2]/1000], 'o', color='red', markersize=5)
ax.plot([r2[0,0]/1000], [r2[0,1]/1000], [r2[0,2]/1000], 'o', color='blue', markersize=5)
ax.plot(r1[:, 0]/1000, r1[:, 1]/1000, r1[:, 2]/1000, label='Molniya', color='red', linewidth=1)
ax.plot(r2[:, 0]/1000, r2[:, 1]/1000, r2[:, 2]/1000, label='Geoestacionaria', color='blue', linewidth=1)

# Add legend outside the plot
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

ax.set_xlabel('x10e3 X [km]')
ax.set_ylabel('x10e3 Y [km]')
ax.set_zlabel('x10e3 Z [km]')
ax.set_title('Problema de  Lambert')

plt.show()


