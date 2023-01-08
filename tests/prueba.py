from src.Orbit import Orbit
from src.CelestialBodies import CelestialBodies
import numpy as np
import matplotlib.pyplot as plt
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

tf = 3600*24 # s
dt = 360 # s

r, v = Orbita.propagate(r0 = r0, v0 = v0, tf = tf, dt = dt)

##################################################
plt.style.use('dark_background')
# plt.rcParams['grid.color'] = "black"

fig = plt.figure()

ax = fig.add_subplot(111, projection='3d')
# ax.w_xaxis.set_pane_color((0.0, 0.0, 0.0, 0.8))
# ax.w_yaxis.set_pane_color((0.0, 0.0, 0.0, 0.6))
# ax.w_zaxis.set_pane_color((0.0, 0.0, 0.0, 1.0))

ax.plot(r[:, 0], r[:, 1], r[:, 2], label='orbit')

ax.plot([0], [0], [0], 'o', label='Earth', color='green', markersize=Tierra.radius/1000)


ax.set_xlabel('X [km]')

ax.set_ylabel('Y [km]')

ax.set_zlabel('Z [km]')

ax.legend()

plt.show()

##################################################33
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


'''
from src.utilities import Utilities

print('--------------------------------------------------')

util = Utilities()

times = util.findTOF(r0 ,r[-1, :], orbit.a)

print('Time of flight: {} s'.format(times))
'''
from src.utilities import Utilities

print('--------------------------------------------------')

util = Utilities()
r0 = [3419.85564, 6019.82602, 2784.60022]
t = 0
r1 = [2935.91195, 6326.18324 , 2660.59584]
t1 = 1*60 + 16.48
r2 = [ 2434.95202, 6597.38674, 2521.52311]
t2 = 2*60 + 33.04


v2, irr = util.Gibbs(r0, r1, r2)

v21 = util.HERRICK_GIBBS(r0, r1, r2, t0, t1, t2)

print('Velocity vector: {} km/s'.format(v2))
print('Velocity vector: {} km/s'.format(v21))

print('--------------------------------------------------')
# Tiempo de vuelo entre dos puntos

r0 = np.array([-2574.9533, 4267.0671, 4431.5026]) # km
r1 = np.array([2700.6738, -4303.5378, -4358.2499]) # km
p = 6681.571

times = util.findTOF(r0, r1, p)
print('Time of flight: {} s'.format(times))
