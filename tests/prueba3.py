import numpy as np

# gravitational constant (m^3/s^2)
G = 6.67430e-11

# masses of the bodies (kg)
masses = np.array([5.972e24, 7.34767309e22])

# positions of the bodies (m)
positions = np.array([[0.0, 0.0, 0.0], [384400000.0, 0.0, 0.0]])

# time step (s)
dt = 60.0

def propagate_orbit(r, v):
  # compute acceleration due to gravity
  a = np.zeros_like(r)
  for i in range(len(masses)):
    ri = positions[i] - r
    a += -G * masses[i] / np.linalg.norm(ri)**3 * ri
  
  # update velocity and position
  v += a * dt
  r += v * dt
  
  return r, v

# initial position (m)
r0 = np.array([-24364566.0, -47457512.0, 2101572.0])

# initial velocity (m/s)
v0 = np.array([5996.0, -3490.0, 7420.0])

# propagate orbit for 10 minutes
for i in range(10):
  r0, v0 = propagate_orbit(r0, v0)
  
print(r0)
print(v0)
