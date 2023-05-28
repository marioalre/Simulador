from src.Coordinates import ecef2latlongh
from src.geomag import Geomag
from src.geopot import Geopot
from src.kepler_propagator import KeplerPropagator

from src.CelestialBodies import CelestialBody

def dataFromTrajectory(tle_file, body, tf, dt):
    '''Magnetic and geopotential data from a trajectory given orbital parameters
    Parameters
    ----------
    orb_param : array_like
        Orbital parameters
    Returns
    -------

    '''

    propagator = KeplerPropagator(body = body)
    
    r0, v0 = propagator.RVFromTLE(tle_file)

    # Propagate orbit
    r, v = propagator.propagate(r0, v0, tf, dt, 0, 0)

    # ECI, necesito ECEF

    # Magnetic field
    mag = Geomag()
    data = mag.get_txt()

    # Geopotential
    geopot = Geopot(body=body) 

    # Change to lat, long, h
    lat = []
    long = []
    h = []

    magVal = []
    geoVal = []

    for r_i in r:
        lati, longi, hi = ecef2latlongh(r_i)
        lat.append(lati)
        long.append(longi)
        h.append(hi)

        # Magnetic field
        val = mag.magnetic_field(10, 10, 7000, N=4)
        magVal.append(val)

        # Geopotential
        valgoe = geopot.gravitational_potential(lati, longi, hi)
        geoVal.append(valgoe)

    return r, magVal, geoVal

