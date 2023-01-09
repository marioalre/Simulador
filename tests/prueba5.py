from src.utilities import *
from src.CelestialBodies import CelestialBodies

earth = CelestialBodies()
earth.earth()

util = Utilities(earth)

lst = util.localSideralTime(2004, 3, 3, 4.5, 139.78333333333)

print(lst)