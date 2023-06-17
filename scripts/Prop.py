import numpy as np
import sys
from src.kepler_propagator import KeplerPropagator
from src.orb2rv import Orb2rv
from src.Propagator import Propagator
from src.CelestialBodies import CelestialBodies


def opciones1b():
    '''Propagar un dt'''
    R, V, dt = valor_inicial()

    cuerpo = elegir_cuerpo_central()
    kepler = KeplerPropagator(cuerpo)
    # Propagar
    # Opciones de propagacion, kepler, y Pkepler, encke
    print("Seleccione una opcion:")
    print("\n")
    print("1. Metodo de Kepler")
    print("2. Metodo de Pkepler")
   #  print("3. Metodo de Encke")

    while True:
        try:
            opcion = input("Opcion: ")

            if opcion == '1':
                print("Metodo de Kepler")
                r, v = kepler.kepler(R, V, dt)
                break
            elif opcion == '2':
                print("Metodo de Pkepler")
                r, v = kepler.Pkepler(R, V, dt, savedata=True)
                break
            else:
                print("Opcion no valida")
                print("Seleccione una opcion:")
                print("\n")
                print("1. Metodo de Kepler")
                print("2. Metodo de Pkepler")
                print("3. Metodo de Encke")
        except ValueError:
            print("Introduzca un valor válido")
            '''
            elif opcion == '3':
                print("Metodo de Encke")
                cuerpo = elegir_cuerpo_central()
                prop_encke = Propagator(cuerpo)
                # Preguntar rangos de tiempo
                trange = []
                r, v = prop_encke.Encke(R, V, dt)
                break
            '''

def opciones2b():
    '''Propagar a lo largo del tiempo'''

    R0, V0, dt, cuerpo = valor_inicial()
    kepler = KeplerPropagator(cuerpo)
    # Propagar
    # Opciones de propagacion, kepler, y Pkepler, encke
    print("Seleccione una opcion:")
    print("\n")
    print("1. Metodo de Kepler")
    print("2. Metodo de Pkepler") 

    while True:
        try:
            opcion = input("Opcion: ")

            if opcion == '1':
                print("Metodo de Kepler")
                prop_kepler = Propagator(cuerpo)
                # Preguntar rangos de tiempo
                tr = trange() 
                kepler.propagate(R0, V0, tr[1], dt, type='kepler')
                break
            elif opcion == '2':
                print("Metodo de Pkepler")
                # Preguntar rangos de tiempo
                tr = trange() 
                kepler.propagate(R0, V0, tr[1], dt, type='Pkepler')
                break

            else:
                print("Opcion no valida")
                print("Seleccione una opcion:")
                print("\n")
                print("1. Metodo de Kepler")
                print("2. Metodo de Pkepler")

        except ValueError:
            print("Introduzca un valor válido")    

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

def trange():
    '''array con tinicial y tfinal'''
    print("introductir el rango de tiempo")
    t0 = float(input("t0: "))
    tf = float(input("tf: "))
    trange = [t0, tf]
    return trange

def rv():
    # Preguntar al usuario vector posicion y velocidad
    r = np.array([float(input("Posicion en x: ")), float(input("Posicion en y: ")), float(input("Posicion en z: "))])
    v = np.array([float(input("Velocidad en x: ")), float(input("Velocidad en y: ")), float(input("Velocidad en z: "))])
    return r, v

def valor_inicial():
    '''Orbita en un dt despues de un instante de tiempo'''

    cuerpo = elegir_cuerpo_central()
    kepler = KeplerPropagator(cuerpo)

    print("Seleccione una opcion:")
    print("\n")
    print("1. Introducir vector posicion y velocidad")
    print("2. Introducir elementos orbitales")
    print("3. Introducir un TLE")
    print("c. Exit. Salir")
    
    while True:
        try:
            opcion = input("Opcion: ")

            dt = float(input("dt: "))

            if opcion == '1':
                print("Introducir vector posicion y velocidad")
                r, v = rv()

                return r, v, dt, cuerpo
            
            elif opcion == '2':

                print("Introducir elementos orbitales")
                a = float(input("Semieje mayor: "))
                e = float(input("Excentricidad: "))
                i = float(input("Inclinacion: "))
                raan = float(input("RAAN: "))
                aop = float(input("Argumento del perigeo: "))
                nu = float(input("Anomalia verdadera: "))

                rrvv = Orb2rv(a=a, e=e, i=i, Omega=raan, omega=aop, nu=nu, body=cuerpo)
                r = rrvv.position_eci()
                v = rrvv.velocity_eci()
                return r, v, dt, cuerpo
            
            elif opcion == '3':
                print("Introducir un TLE")
                # preguntar al usuario el path del TLE
                path = input("Path del TLE: ")
                r, v = kepler.RVFromTLE(path)
                return r, v, dt, cuerpo
            
            elif opcion == 'c':
                print("Saliendo del programa")
                sys.exit()
            else:
                print("Opcion no valida")
                print("Seleccione una opcion:")
                print("\n")
                print("1. Introducir vector posicion y velocidad")
                print("2. Introducir elementos orbitales")
                print("3. Introducir un TLE")

        except ValueError:
            print("Introduzca un valor válido")
