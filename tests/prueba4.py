from src.Orbit import Orbit
from src.CelestialBodies import CelestialBodies
import numpy as np
import matplotlib.pyplot as plt
from src.rv2orb import Rv2orb
from src.orb2rv import Orb2rv
import src.Coordinates as co

Tierra = CelestialBodies()
Tierra.earth()

# Orbit
Orbita = Orbit(Tierra)

# Molnya parameters
# Specific energy with a and e

orb = Orb2rv(a = 26554, e = 0.72, Omega=0, i = 63.4, omega = -90, nu = 0, body=Tierra)
r0 = orb.position_eci()
v0 = orb.velocity_eci()

# Geo
orb = Orb2rv(a = 42165, e = 0, Omega=0, i = 0, omega = 0, nu = 0, body=Tierra)
r01 = orb.position_eci()
v01 = orb.velocity_eci()

# GPS
orb = Orb2rv(e = 0.00658, Omega=139.6, i = 54.9531, omega = 230.4953, M = 128.7809, n=0.000145852075 ,body=Tierra)
r02 = orb.position_eci()
v02 = orb.velocity_eci()



# Propagate the orbit

tf = 3600*24 # s
dt = 60*4 # s

r, v = Orbita.propagate(r0 = r0, v0 = v0, tf = tf, dt = dt)
r1, v1 = Orbita.propagate(r0 = r01, v0 = v01, tf = tf, dt = dt*2)
r2, v2 = Orbita.propagate(r0 = r02, v0 = v02, tf = tf, dt = dt)

##################################################


# plt.rcParams['grid.color'] = "black"

fig = plt.figure()

ax = fig.add_subplot(111, projection='3d')
# ax.w_xaxis.set_pane_color((0.0, 0.0, 0.0, 0.0))
# ax.w_yaxis.set_pane_color((0.0, 0.0, 0.0, 0.0))
# ax.w_zaxis.set_pane_color((0.0, 0.0, 0.0, 0.0))

ax.plot(r[:, 0]/1000, r[:, 1]/1000, r[:, 2]/1000, label='Molniya', color='red', linewidth=1)
ax.plot(r1[:, 0]/1000, r1[:, 1]/1000, r1[:, 2]/1000, label='Geoestacionaria', color='blue', linewidth=1)
ax.plot(r2[:, 0]/1000, r2[:, 1]/1000, r2[:, 2]/1000, label='GPS', color='orange', linewidth=1)
ax.plot([0], [0], [0], 'o', label='Earth', color='green', markersize=Tierra.radius/1000)

ax.set_xlabel('x10e3 X [km]')
ax.set_ylabel('x10e3 Y [km]')
ax.set_zlabel('x10e3 Z [km]')
ax.set_title('Órbitas Keplerianas')

ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

plt.show()

##################################################

r, v = co.cart2efix(r, v, tf, dt, 12720)
r1, v1 = co.cart2efix(r1, v1, tf, dt*2, 12720)
r2, v2 = co.cart2efix(r2, v2, tf, dt, 12720)

##################################################
lat, long = co.ecef2latlong(r)
lat1, long1 = co.ecef2latlong(r1)
lat2, long2 = co.ecef2latlong(r2)

##################################################
m = co.plot_ground_track(lat, long)
m = co.plot_ground_track(lat1, long1, m)
m = co.plot_ground_track(lat2, long2, m)

# Guarda el mapa en un archivo HTML
m.save('results/satellite_track.html')

# plt.rcParams['grid.color'] = "black"

fig1 = plt.figure()

ax1 = fig1.add_subplot(111, projection='3d')
# ax.w_xaxis.set_pane_color((0.0, 0.0, 0.0, 0.0))
# ax.w_yaxis.set_pane_color((0.0, 0.0, 0.0, 0.0))
# ax.w_zaxis.set_pane_color((0.0, 0.0, 0.0, 0.0))

ax1.plot(r[:, 0]/1000, r[:, 1]/1000, r[:, 2]/1000, label='Molniya', color='red', linewidth=1)
ax1.plot([r1[-1, 0]/1000], [r1[-1, 1]/1000], [r1[-1, 2]/1000], 'o' ,label='Geoestacionaria', color='blue', markersize=5)
ax1.plot(r2[:, 0]/1000, r2[:, 1]/1000, r2[:, 2]/1000, label='GPS', color='orange', linewidth=1)
ax1.plot([0], [0], [0], 'o', label='Earth', color='green', markersize=Tierra.radius/1000)


ax1.set_xlabel('x10e3 X [km]')
ax1.set_ylabel('x10e3 Y [km]')
ax1.set_zlabel('x10e3 Z [km]')
ax1.set_title('Órbitas Keplerianas')

ax1.legend(bbox_to_anchor=(1, 1), loc='upper left', borderaxespad=0.)

plt.show()
##################################################

t0 = 0
t = 3600

# Orbit

orbit = Rv2orb(r[-1, :], v[-1, :], Tierra, t0- t)

# Print results
print('Semi-major axis: {} km'.format(orbit.a))
print('Eccentricity: {}'.format(np.linalg.norm(orbit.e)))
print('Inclination: {} deg'.format(np.rad2deg(orbit.i)))
print('Ascending node: {} deg'.format(np.rad2deg(orbit.Omega)))
print('Argument of periapsis: {} deg'.format(np.rad2deg(orbit.omega)))
print('True anomaly: {} deg'.format(np.rad2deg(orbit.nu)))