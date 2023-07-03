from src.Coordinates import *
import numpy as np
import pandas as pd

# leer datos de results/propagate.csv
data = pd.read_csv('results/propagate.csv')
# print(data)

# Extraer las tres primeras columnas
r = np.array([data['x'], data['y'], data['z']]).T
v = np.array([data['vx'], data['vy'], data['vz']]).T
a = np.zeros_like(r)
# Clase para convertir coordenadas ECI a ECEF

ttt = 0.0426236319
jdut1 = 2453101.827406783
lod = 0.001556300
xp = -6.820455828585174e-07
yp = 1.615927632369383e-06
eqeterms = 2
ddpsi = -2.530485008551223e-7
ddeps = -1.878653014299452e-8

eci2ecef = ECI2ECEF(ttt, jdut1, lod, xp, yp, eqeterms, ddpsi, ddeps)
latt = []
lont = []
altt = []


# r, v = cart2efix(r, v, 3600, 180, 12720)

for i in range(len(r)):
    # Convertir a ECEF
    r[i], v[i], a[i] = eci2ecef.eci2ecef(r[i], v[i], a[i])


    # Cambio de coordenadas a latitud, longitud y altura
    longd, longc, lat, alt = ecef2latlongh(r[i])
    
    latt.append(lat)
    lont.append(longc)
    altt.append(alt)

# save to csv
data = {'x': r[:, 0], 'y': r[:, 1], 'z': r[:, 2], 'vx': v[:, 0], 'vy': v[:, 1], 'vz': v[:, 2]}
df = pd.DataFrame(data)
df.to_csv('results/propagate_ecef.csv')

data = {'lat': np.degrees(latt), 'long': np.degrees(lont), 'alt': altt}
df = pd.DataFrame(data)
df.to_csv('results/propagate_lla.csv')

