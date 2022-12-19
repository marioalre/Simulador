from astropy import units as u
import numpy as np

from poliastro.bodies import Earth, Mars, Sun
from poliastro.twobody import Orbit

r = [2855, -35639, -24309] << u.km
v = [-1.656, -4.455, -3.2803] << u.km / u.s

orb = Orbit.from_vectors(Earth, r, v)


print(orb.a)
print(orb.ecc)
print(orb.inc * 180 / np.pi)
print(orb.raan * 180 / np.pi)
print(orb.argp  * 180 / np.pi)
print(orb.nu * 180 / np.pi)

