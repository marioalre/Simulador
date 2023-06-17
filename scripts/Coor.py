import numpy as np
import sys
import pandas as pd
from src.Coordinates import *
from src.CelestialBodies import CelestialBodies
from src.rv2orb import Rv2orb
from src.orb2rv import Orb2rv

def opciones1c(type):
    '''Eci to ecef


        This function initializes the class with the input parameters.
    Inputs:
        ttt: Julian centuries of terrestrial time
        jdut1: Julian date of UT1
        lod: Length of day (sec)
        xp: Polar motion coefficient (arcsec)
        yp: Polar motion coefficient (arcsec)
        eqeterms: Boolean to select if the extra terms for nutation are used (0 or 2)
        ddpsi: Correction to delta psi (arcsec)
        ddeps: Correction to delta eps (arcsec)
    '''
    # print del comentario de arriba

    print('ECI to ECEF')
    print('This function initializes the class with the input parameters.')
    print('Inputs:')
    print('    ttt: Julian centuries of terrestrial time')
    print('    jdut1: Julian date of UT1')
    print('    lod: Length of day (sec)')
    print('    xp: Polar motion coefficient (arcsec)')
    print('    yp: Polar motion coefficient (arcsec)')
    print('    eqeterms: Boolean to select if the extra terms for nutation are used (0 or 2)')
    print('    ddpsi: Correction to delta psi (arcsec)')
    print('    ddeps: Correction to delta eps (arcsec)')

    # preguntar por los parametros de entrada o si se quiere usar los valores por defecto
    print('Do you want to use the default values?')
    print('    1. Yes')
    print('    2. No')

    while True:
        try:
            opcion = input('Option: ')
            if opcion == '1':
                print('Using default values')
                ttt = 0.0426236319
                jdut1 = 2453101.827406783
                lod = 0.001556300
                xp = -6.820455828585174e-07
                yp = 1.615927632369383e-06
                eqeterms = 2
                ddpsi = -2.530485008551223e-7
                ddeps = -1.878653014299452e-8
                break
            elif opcion == '2':
                print('Ask for the input parameters')
                ttt = float(input('Julian centuries of terrestrial time (centuries) = '))
                jdut1 = float(input('Julian date of UT1 (days from 4713 BC) = '))
                lod = float(input('Length of day (sec) = '))
                xp = float(input('x Polar motion coefficient (arcsec) = '))
                yp = float(input('y Polar motion coefficient (arcsec) = '))
                eqeterms = int(input('Boolean to select if the extra terms for nutation are used (0 or 2) = '))
                ddpsi = float(input('Correction to delta psi (arcsec) = '))
                ddeps = float(input('Correction to delta eps (arcsec) = '))
                break
            else:
                print('Option not valid')
                print('Do you want to use the default values?')
                print('    1. Yes')
                print('    2. No')
        except ValueError:
            print('Introduce a valid value')

    print('\n')

    # Initialize the class
    coord = ECI2ECEF(ttt, jdut1, lod, xp, yp, eqeterms, ddpsi, ddeps)

    # Calculate the transformation matrix
    print('\n')
    if type == '1':
        print('Ask for the vector position and velocity ECIs')
        r, v, a = rva()
        r, v, a = coord.eci2ecef(r, v, a)
        text = 'ECEF'
    elif type == '2':
        print('Ask for the vector position and velocity ECEF')
        r, v, a = rva()
        r, v, a = coord.ecef2eci(r, v, a)
        text = 'ECI'

    print('\n')
    print('Los vextores posición, velocidad y aceleración ' + text + ' son:')
    print('    r = ' + str(r) + ' km')
    print('    v = ' + str(v) + ' km/s')
    print('    a = ' + str(a) + ' km/s2')

    # Save the results in a pandas dataframe
    data = {'r': r, 'v': v, 'a': a}
    df = pd.DataFrame(data=data)
    df.to_csv('results/rva' + text +'.csv', index=False)

def opciones2c():
    '''Rv to coe'''
    print('RV to COE')

    # Ask for the vector position and velocity
    r, v, a = rva()

    cuerpo = elegir_cuerpo_central()
    # Initialize the class
    coord = Rv2orb(r, v, cuerpo)

    print('Resultados')
    print('    a = ' + str(coord.a) + ' km')
    print('    e = ' + str(coord.e))
    print('    i = ' + str(coord.i) + ' deg')
    print('    RAAN = ' + str(coord.Omega) + ' deg')
    print('    omega = ' + str(coord.omega) + ' deg')
    print('    nu = ' + str(coord.nu) + ' deg')
    print('    n = ' + str(coord.n) + ' rad/s')
    print('    T = ' + str(coord.T) + ' s')

    # Save the results in a pandas dataframe
    data = {'a': coord.a, 'e': coord.e, 'i': coord.i, 'RAAN': coord.raan, 'omega': coord.omega, 'nu': coord.nu, 'n': coord.n, 'T': coord.T}
    df = pd.DataFrame(data=data)
    df.to_csv('results/orbital_elements.csv', index=False)
    print('Results saved in results/orbital_elements.csv')

def opciones3c():
    '''Coe to rv'''
    print('COE to RV')

    # Ask for the orbital elements
    a = float(input('Semi-major axis (km) = '))
    e = float(input('Eccentricity = '))
    i = float(input('Inclination (deg) = '))
    Omega = float(input('RAAN (deg) = '))
    omega = float(input('Argument of perigee (deg) = '))
    nu = float(input('True anomaly (deg) = '))
    cuerpo = elegir_cuerpo_central()

    # Initialize the class
    coord = Orb2rv(a=a, e=e, i=i, Omega=Omega, omega=omega, nu=nu, body=cuerpo)

    # Calculate the position and velocity vectors
    r = coord.position_eci()
    v = coord.velocity_eci()

    print('Resultados')
    print('    r = ' + str(r) + ' km')
    print('    v = ' + str(v) + ' km/s')

    # Save the results in a pandas dataframe
    data = {'r': r, 'v': v}
    df = pd.DataFrame(data=data)
    df.to_csv('results/rv.csv', index=False)
    print('Results saved in results/rv.csv')

def opciones4c():
    '''ecef to lat, long and h'''

    print('ECEF to lat, long and h')

    # Ask for the vector position and velocity
    r, v, a = rva()

    cuerpo = elegir_cuerpo_central()
    # Initialize the class
    long_gd, long_gc, lambda_, hell = ecef2latlongh(r)

    print('Resultados')
    print('    Longitud geodésica = ' + str(long_gd) + ' deg')
    print('    Longitud geocéntrica = ' + str(long_gc) + ' deg')
    print('    Latitud = ' + str(lambda_) + ' deg')
    print('    Altura = ' + str(hell) + ' km')

    # Save the results in a pandas dataframe
    data = {'Longitud geodésica': long_gd, 'Longitud geocéntrica': long_gc, 'Latitud': lambda_, 'Altura': hell}
    df = pd.DataFrame(data=data)
    df.to_csv('results/latlongh.csv', index=False)
    print('Results saved in results/latlongh.csv')


def rva():
    # Ask for the vector position and velocity and acceleration
    while True:
        try:
            print('Ask for the vector position and velocity')
            rx = float(input('Position in x (km) = '))
            ry = float(input('Position in y (km) = '))
            rz = float(input('Position in z (km) = '))
            vx = float(input('Velocity in x (km/s) = '))
            vy = float(input('Velocity in y (km/s) = '))
            vz = float(input('Velocity in z (km/s) = '))
            ax = float(input('Acceleration in x (km/s2) = '))
            ay = float(input('Acceleration in y (km/s2) = '))
            az = float(input('Acceleration in z (km/s2) = '))

            r = np.array([rx, ry, rz])
            v = np.array([vx, vy, vz])
            a = np.array([ax, ay, az])

            return r, v ,a
        except ValueError:
            print('Introduce a un valor válido')


def elegir_cuerpo_central():
    '''Elegir el cuerpo central'''
    print("Seleccione una opcion:")
    print("\n")
    print("1. Tierra")
    print("2. Luna")
    print("3. Sol")
    print("4. Jupiter")
    print("5. Saturno")
    print("6. Marte")
    print("7. Venus")
    print("8. Mercurio")
    print("9. Pluton")
    print("10. Neptuno")
    print("11. Urano")
    print("12. Salir")

    cuerpo = CelestialBodies()
    while True:
        try:
            opcion = input("Opcion: ")
            if opcion == '1':
                print("Tierra")
                cuerpo.earth()
                return cuerpo
            elif opcion == '2':
                print("Luna")
                cuerpo.moon()
                return cuerpo
            elif opcion == '3':
                print("Sol")
                cuerpo.sun()
                return cuerpo
            elif opcion == '4':
                print("Jupiter")
                cuerpo.jupiter()
                return cuerpo
            elif opcion == '5':
                print("Saturno")
                cuerpo.saturn()
                return cuerpo
            elif opcion == '6':
                print("Marte")
                cuerpo.mars()
                return cuerpo
            elif opcion == '7':
                print("Venus")
                cuerpo.venus()
                return cuerpo
            elif opcion == '8':
                print("Mercurio")
                cuerpo.mercury()
                return cuerpo
            elif opcion == '9':
                print("Pluton")
                cuerpo.pluto()
                return cuerpo
            elif opcion == '10':
                print("Neptuno")
                cuerpo.neptune()
                return cuerpo
            elif opcion == '11':
                print("Urano")
                cuerpo.uranus()
                return cuerpo
            elif opcion == '12':
                print("Saliendo del programa")
                sys.exit()
            else:
                print("Opcion no valida")
                print("Seleccione una opcion:")
                print("\n")
                print("1. Tierra")
                print("2. Luna")
                print("3. Sol")
                print("4. Jupiter")
                print("5. Saturno")
                print("6. Marte")
                print("7. Venus")
                print("8. Mercurio")
                print("9. Pluton")
                print("10. Neptuno")
                print("11. Urano")
        except ValueError:
            print("Introduzca un valor válido")