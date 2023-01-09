from src.utilities import *
from src.CelestialBodies import CelestialBodies

earth = CelestialBodies()
earth.earth()

util = Utilities(earth)

lst = util.localSideralTime(2004, 3, 3, 4.5, 139.78333333333)

print(lst)

r, v = util.rv_from_observation(2551, 0, 90, 0.1130, 30, 0.05651, 300, 60, 0)

print(r)
print(v)