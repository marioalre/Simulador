from src.utilities import Utilities
from src.CelestialBodies import CelestialBodies
import numpy as np
import sys
import pandas as pd

util = Utilities()

def opciones1d():
    '''Elegir entre Gibbs o Henry Gibbs'''

    print("Seleccione una opcion:")
    print("\n")
    print("1. Metodo de Gibbs")
    print("2. Metodo de Henrry-Gibbs")

    while True:
        try:
            opcion = input("Opcion: ")

            if opcion == '1':
                print("Metodo de Gibbs")
                # Preguntar los vectores
                r1, r2, r3 = vectoress()
                # Calcular velocidad
                v2 = util.Gibbs(r1, r2, r3)
                print("Velocidad: ", v2)
                break
            elif opcion == '2':
                print("Metodo de Henry-Gibbs")
                # Preguntar los vectores
                r1, r2, r3 = vectoress()
                t1, t2, t3 = tiempos()
                # Calcular velocidad
                v2 = util.HERRICK_GIBBS(r1, r2, r3, t1, t2, t3)
                print("Velocidad: ", v2)
                break
            else:
                print("Opcion no valida")
                print("Seleccione una opcion:")
                print("\n")
                print("1. Metodo de Gibbs")
                print("2. Metodo de Henry-Gibbs")
        except ValueError:
            print("Introduzca un valor válido")

def opciones2d():
    '''This function calculates the position and velocity vectors from the observation of a satellite.
        Parameters
        ----------
        rho : float
            Range in km
        drho : float
            Range rate in km/s
        A : float
            Azimuth in degrees
        dA : float
            Azimuth rate in degrees/s
        a : float
            Elevation in degrees
        da : float
            Elevation rate in degrees/s
        theta : float
            Local sideral time in degrees
        phi : float
            Latitude in degrees
        H : float
            Height in km
        Returns
        ----------
        R : numpy array
            Position vector in km
        V : numpy array
            Velocity vector in km/s
        '''
    
    print("Introductir los datos de la observacion")

    while True:
        try:
            rho = float(input("Rango [km]: "))
            drho = float(input("Rango rate [km/s]: "))
            A = float(input("Azimuth [deg]: "))
            dA = float(input("Azimuth rate [deg/s]: "))
            a = float(input("Elevation [deg]: "))
            da = float(input("Elevation rate [deg/s]: "))
            theta = float(input("Local sideral time [deg]: "))
            phi = float(input("Latitude [deg]: "))
            H = float(input("Height [km]: "))

            # Calcular posicion y velocidad
            R, V = util.rv_from_observation(rho, drho, A, dA, a, da, theta, phi, H)
            print("Posicion: ", R)
            print("Velocidad: ", V)

            data = {'R': R, 'V': V}
            df = pd.DataFrame(data, index=['x', 'y', 'z'])
            df.to_csv('results/posvel_from_ell_rates.csv')
            print("Datos guardados en results/posvel_from_ell_rates.csv")
            break
        except ValueError:
            print("Introduzca un valor válido")

def opciones3d():
    print("Gauss' method with iterative improvement")

    while True:
        try:
            # Preguntar los vectores
            print("Introduzca los vectores posición")
            r1, r2, r3 = vectoress()
            print("Introduzca los vectores desde el punto del observador")
            rho1, rho2, rho3 = vectoress()
            print("Introduzca los tiempos")
            t1, t2, t3 = tiempos()
            # Calcular velocidad
            util.Gauss_POD(r1, r2, r3, rho1, rho2, rho3, t1, t2, t3)
            r2, v2 = util.iterative()
            print("Posicion: ", r2, "km")
            print("Velocidad: ", v2, "km/s")

            data = {'R': r2, 'V': v2}
            df = pd.DataFrame(data, index=['x', 'y', 'z'])
            df.to_csv('results/posvel_from_gauss.csv')
            print("Datos guardados en results/posvel_from_gauss.csv")
            break
        except ValueError:
            print("Introduzca un valor válido")


def vectoress():
    '''Preguntar por 3 vectores posición'''

    r1 = []
    r2 = []
    r3 = []

    while True:
        try:
            print("Introduzca los vectores posición")
            print("Vector 1:")
            x1 = float(input("x: "))
            y1 = float(input("y: "))
            z1 = float(input("z: "))
            r1 = np.array([x1, y1, z1])
            print("Vector 2:")
            x2 = float(input("x: "))
            y2 = float(input("y: "))
            z2 = float(input("z: "))
            r2 = np.array([x2, y2, z2])
            print("Vector 3:")
            x3 = float(input("x: "))
            y3 = float(input("y: "))
            z3 = float(input("z: "))
            r3 = np.array([x3, y3, z3])

            return r1, r2, r3
        except ValueError:
            print("Introduzca un valor válido")

def tiempos():
    '''preguntar por 3 tiempos'''
    print("Introduzca los tiempos")
    t1 = []
    t2 = []
    t3 = []

    while True:
        try:
            print("Introduzca los tiempos")
            print("Tiempo 1:")
            t1 = float(input("t: "))
            print("Tiempo 2:")
            t2 = float(input("t: "))
            print("Tiempo 3:")
            t3 = float(input("t: "))

            return t1, t2, t3
        except ValueError:
            print("Introduzca un valor válido")