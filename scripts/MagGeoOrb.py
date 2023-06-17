from src.Coordinates import ecef2latlongh
from src.geomag import Geomag
from src.geopot import Geopot
from src.kepler_propagator import KeplerPropagator
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

from src.CelestialBodies import CelestialBodies

def dataFromTrajectory(tle_file, tf, dt, savedata = False):
    '''Magnetic and geopotential data from a trajectory given orbital parameters
    Parameters
    ----------
    orb_param : array_like
        Orbital parameters
    Returns
    -------

    '''
    earth = CelestialBodies()
    earth.earth()
    propagator = KeplerPropagator(body = earth)
    
    r0, v0 = propagator.RVFromTLE(tle_file)

    # Propagate orbit
    r, v = propagator.propagate(r0, v0, tf, dt, 0, 0)

    # ECI, necesito ECEF

    # Magnetic field
    mag = Geomag()

    # Geopotential
    geopot = Geopot() 

    # Change to lat, long, h
    lat = []
    long = []
    h = []

    magVal = []
    geoValpot = []
    geoValGrav = []

    for r_i in r:
        longgi, longi, lati, hi = ecef2latlongh(r_i)
        lat.append(np.degrees(lati)) # degrees
        long.append(np.degrees(longi)) # degrees
        h.append(hi) # km

        rn = np.linalg.norm(r_i)
        # Magnetic field
        val = mag.magnetic_field(rn, lati, longi, N=13)
        val = mag.transformation2NED(val[1])
        magVal.append(val)

        # Geopotential
        valgeopot, valgeograv  = geopot.gravitational_potential(rn, lati, longi, 30)
        geoValpot.append(valgeopot)
        geoValGrav.append(valgeograv)

    if savedata:
        # Create a pandas dataframe to save the data

        df = pd.DataFrame({'lat': lat, 'long': long, 'h': h, 'Bx': [x[0] for x in magVal], 'By': [x[1] for x in magVal], 'Bz': [x[2] for x in magVal], 'U': geoValpot, 'g1':[x[0] for x in geoValGrav], 'g2':[x[1] for x in geoValGrav], 'g3':[x[2] for x in geoValGrav]})
        df.to_csv('results/MagGeo.csv')


    return r, magVal, geoValpot, geoValGrav


def dataFromRV(R0, V0, tf, dt, savedata = False, printdata = True):
    '''Magnetic and geopotential data from a trajectory given orbital parameters
    Parameters
    ----------
    R0 : array_like
        Initial position vector in km
    V0 : array_like
        Initial velocity vector in km/s
    tf : float
        Final time in seconds
    dt : float
        Time step in seconds
    savedata : bool
        Save data in a csv file
    printdata : bool
        Print data
    Returns
    -------

    '''
    earth = CelestialBodies()
    earth.earth()
    propagator = KeplerPropagator(body = earth)

    # Propagate orbit
    r, v = propagator.propagate(R0, V0, tf, dt, 0, 0)

    # ECI, necesito ECEF

    # Magnetic field
    mag = Geomag()

    # Geopotential
    geopot = Geopot() 

    # Change to lat, long, h
    lat = []
    long = []
    h = []

    magVal = []
    geoValpot = []
    geoValGrav = []

    for r_i in r:
        longgi, longi, lati, hi = ecef2latlongh(r_i)
        lat.append(np.degrees(lati)) # degrees
        long.append(np.degrees(longi)) # degrees
        h.append(hi) # km

        rn = np.linalg.norm(r_i)
        # Magnetic field
        val = mag.magnetic_field(rn, lati, longi, N=13)
        val = mag.transformation2NED(val[1])
        magVal.append(val)

        # Geopotential
        valgeopot, valgeograv  = geopot.gravitational_potential(rn, lati, longi, 30)
        geoValpot.append(valgeopot)
        geoValGrav.append(valgeograv)

    if savedata:
        # Create a pandas dataframe to save the data

        df = pd.DataFrame({'lat': lat, 'long': long, 'h': h, 'Bx': [x[0] for x in magVal], 'By': [x[1] for x in magVal], 'Bz': [x[2] for x in magVal], 'U': geoValpot, 'g1':[x[0] for x in geoValGrav], 'g2':[x[1] for x in geoValGrav], 'g3':[x[2] for x in geoValGrav]})
        df.to_csv('results/MagGeo.csv')

    if printdata:
        print('Magnetic field')

        mag.plotMagneticFieldGeneral(np.arange(0, tf+dt, dt, dtype=float), magVal)

    return r, magVal, geoValpot, geoValGrav


if __name__ == "__main__":

    # dataFromTrajectory('data/paz.tle', 3600, 60, savedata = True)

    # R y V de una orbita circular de altura = 500km y un inclinacion de 60 grados

    dataFromRV([6878.137, 0, 0], [0,  7.61133671, 0.13912824], 3600*10, 60, savedata = True, printdata = True)



