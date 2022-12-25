from src.geopot import Geopot
from src.CelestialBodies import CelestialBodies

if __name__ == '__main__':
    Tierra = CelestialBodies()
    Tierra.earth()

    Potencial = Geopot(Tierra, 7000000, 1, 1)

    Potencial.gravitational_potential()