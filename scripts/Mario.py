# Pantalla de bienvenida del simulador con un menu para las diferentes clases que hay en src
# Autor: Mario 

import os
import sys
from scripts.Mag import opciones1, opciones2, opciones3, opciones4
from scripts.Geo import opciones1a, opciones2a
from scripts.Prop import opciones1b, opciones2b
from scripts.Coor import opciones1c, opciones2c, opciones3c, opciones4c
from scripts.DPO import opciones1d, opciones2d, opciones3d


# pantalla de bienvenida del simulador con prints
def bienvenida():
    '''Pantalla de bienvenida del simulador 
    El trabajo tiene por objetivo el desarrollo de un software
    que sea capaz de resolver entre otras cosas, los principales
    parámetros orbitales, valores del campo magnético, estrellas
    visibles, ventanas de tiempo, etc. para su aplicación a una
    determinada misión espacial, así como la reproducción de las 
    condiciones del entorno espacial en la Tierra. 
    Para su elaboración se utilizará Python'''
    # Clear the terminal
    os.system('cls' if os.name == 'nt' else 'clear')

    print("Bienvenido al simulador de mecánica espacial")
    print("\n")
    print("El trabajo tiene por objetivo el desarrollo de un software")
    print("que sea capaz de resolver entre otras cosas, los principales")
    print("parámetros orbitales, valores del campo magnético, estrellas")
    print("visibles, ventanas de tiempo, etc. para su aplicación a una")
    print("determinada misión espacial, así como la reproducción de las ")
    print("condiciones del entorno espacial en la Tierra. ")
    print("Para su elaboración se utilizará Python")
    print('')


def menu():
    '''Menu principal del simulador'''
    print("\n")
    print("Seleccione una opcion:")
    print("\n")
    print("1. Campo magnetico")
    print("2. Campo goeopotencial")
    print("3. Propagadores orbitas")
    print("4. Cambio de coordenadas")
    print("5. Determinación preliminar de la órbita")
    print("6. Otras utilidades")
    print("7. Salir")

    while True:
        opcion = input("Opcion: ")
        if opcion == '1':
            print("Campo magnetico")
            menu_campo_magnetico()
            break
        elif opcion == '2':
            print("Campo goeopotencial")
            menu_geopot()
            break
        elif opcion == '3':
            print("Propagadores orbitas")
            meno_propagadores()
            break
        elif opcion == '4':
            print("Cambio de coordenadas")
            menu_cambio_coordenadas()
            break
        elif opcion == '5':
            print("Determinación preliminar de la órbita")
            menu_determinacion_preliminar()
            break
        elif opcion == '6':
            print("Otras utilidades")
            print("En proceso de desarrollo")
            # menu_otras_utilidades()
            break
        elif opcion == '7':
            print("Saliendo del programa")
            sys.exit()
        else:
            print("Opcion no valida")
            print("Seleccione una opcion:")
            print("\n")
            print("1. Campo magnetico")
            print("2. Campo goeopotencial")
            print("3. Propagadores orbitas")
            print("4. Cambio de coordenadas")
            print("5. Determinación preliminar de la órbita")
            print("6. Otras utilidades")

    return opcion

def menu_determinacion_preliminar():
    '''Menu de la determinacion preliminar de la orbita'''
    print("\n")
    print("Bienvenido al simulador de la determinacion preliminar de la orbita")
    print("\n")
    print("Se puede calcular la orbita de un satelite a partir de un conjunto de observaciones")
    print("Se utiliza el modelo SGP4")
    print("\n")
    print("Seleccione una opcion:")
    print("\n")
    print("1. 3 observaciones del vector posicion") 
    print("2. observaciones de latitud, longitud, altura y variacion de latitud, longitud y altura")
    print("3. observaciones desde la Tierra, desde el centro de la Tierra y tiempo")
    print("4. observaciones de latitud, longitud, altura")
    print("5 Salir")

    while True:
        opcion = input("Opcion: ")
        if opcion == '1':
            print("3 observaciones del vector posicion")
            print("En proceso de desarrollo")
            opciones1d()
            break
        elif opcion == '2':
            print("observaciones de latitud, longitud, altura y variacion de latitud, longitud y altura")
            print("En proceso de desarrollo")
            opciones2d()
            break
        elif opcion == '3':
            print("observaciones desde la Tierra, desde el centro de la Tierra y tiempo")
            print("En proceso de desarrollo")
            opciones3d()
            break
        elif opcion == '4':
            print("observaciones de latitud, longitud, altura")
            print("En proceso de desarrollo")
            # opciones1c('4')
            break
        elif opcion == '5':
            print("Saliendo del programa")
            sys.exit()
        else:
            print("Opcion no valida")
            print("Seleccione una opcion:")
            print("\n")
            print("1. 3 observaciones del vector posicion") 
            print("2. observaciones de latitud, longitud, altura y variacion de latitud, longitud y altura")
            print("3. observaciones desde la Tierra, desde el centro de la Tierra y tiempo")
            print("4. observaciones de latitud, longitud, altura")


def menu_cambio_coordenadas():
    '''Menu del cambio de coordenadas'''
    print("\n")
    print("Bienvenido al simulador del cambio de coordenadas")
    print("\n")
    print("Se puede calcular el cambio de coordenadas de un sistema de referencia a otro")
    print("Se utiliza el modelo IAU-76/FK5")
    print("\n")
    print("Seleccione una opcion:")
    print("\n")
    print("1. Cambio de coordenadas ECI a ECEF")
    print("2. Cambio de coordenadas ECEF a ECI")
    print("3. RV a COE")
    print("4. COE a RV")
    print("5. ECEF a latitud, longitud y altura")

    while True:
        opcion = input("Opcion: ")
        if opcion == '1':
            print("Cambio de coordenadas ECI a ECEF")
            print("En proceso de desarrollo")
            opciones1c('1')
            break
        elif opcion == '2':
            print("Cambio de coordenadas ECEF a ECI")
            print("En proceso de desarrollo")
            opciones1c('2')
            break
        elif opcion == '3':
            print("RV a COE")
            print("En proceso de desarrollo")
            opciones2c()
            break
        elif opcion == '4':
            print("COE a RV")
            print("En proceso de desarrollo")
            opciones3c()
            break
        elif opcion == '5':
            print("ECEF a latitud, longitud y altura")
            print("En proceso de desarrollo")
            opciones4c()
            break
        elif opcion == 'c':
            print("Saliendo del programa")
            sys.exit()
        else:
            print("Opcion no valida")
            print("Seleccione una opcion:")
            print("\n")
            print("1. Cambio de coordenadas ECI a ECEF")
            print("2. Cambio de coordenadas ECEF a ECI")
            print("3. RV a COE")
            print("4. COE a RV")
            print("5. ECEF a latitud, longitud y altura")


def meno_propagadores():
    '''Menu de los propagadores de orbitas'''
    print("\n")
    print("Bienvenido al simulador de propagadores de orbitas")
    print("\n")
    print("Se puede calcular la orbita de un satelite en un instante de tiempo o a lo largo del tiempo")
    print("Se utiliza el modelo SGP4")
    print("\n")
    print("Seleccione una opcion:")
    print("\n")
    print("1. Orbita en un instante de tiempo")
    print("2. Orbita a lo largo del tiempo")
    print("3. Ploteo de los resultados")
    print("c. Exit. Salir")

    while True:
        opcion = input("Opcion: ")
        if opcion == '1':
            print("Orbita en un dt despues de un instante de tiempo")
            opciones1b()
            break
        elif opcion == '2':
            print("Orbita a lo largo del tiempo")
            opciones2b()
            break
        elif opcion == '3':
            print("Ploteo de los resultados")
            print("En proceso de desarrollo")
            # opciones3b()
            break
        elif opcion == 'c':
            print("Saliendo del programa")
            sys.exit()
        else:
            print("Opcion no valida")
            print("Seleccione una opcion:")
            print("\n")
            print("1. Orbita en un instante de tiempo")
            print("2. Orbita a lo largo del tiempo")

def menu_geopot():
    '''Menu del campo geopotencial'''
    print("\n")
    print("Bienvenido al simulador del campo geopotencial")
    print("\n")
    print("Se puede calcular el campo geopotencial en un punto, en un conjunto de puntos o a lo largo del tiempo en un punto")
    print("Se utiliza el modelo EGM96")
    print("\n")
    print("Seleccione una opcion:")
    print("\n")
    print("1. Campo geopotencial en un punto")
    print("2. Campo geopotencial en un conjunto de puntos")
    print("3. Campo geopotencial a lo largo del tiempo en un punto")
    print("4. Ploteo de los resultados")
    print("c. Exit. Salir")

    while True:
        opcion = input("Opcion: ")
        if opcion == '1':
            print("Campo geopotencial en un punto")
            opciones1a()
            break
        elif opcion == '2':
            print("Campo geopotencial en un conjunto de puntos")
            opciones2a()
            break
        elif opcion == '3':
            print("En proceso de desarrollo")
            # opciones3a()
            break
        elif opcion == 'c':
            print("Saliendo del programa")
            sys.exit()
        else:
            print("Opcion no valida")
            print("Seleccione una opcion:")
            print("\n")
            print("1. Campo geopotencial en un punto")
            print("2. Campo geopotencial en un conjunto de puntos")
            print("3. Campo geopotencial a lo largo del tiempo en un punto")

def menu_campo_magnetico():
    '''Menu del campo magnetico'''
    print("\n")
    print("Bienvenido al simulador del campo magnetico")
    print("\n")
    print("Se puede calcular el campo magnetico en un punto, en un conjunto de puntos o a lo largo del tiempo en un punto")
    print("Se utiliza el modelo IGRF-13")
    print("\n")
    print("Seleccione una opcion:")
    print("\n")
    print("1. Campo magnetico en un punto")
    print("2. Campo magnetico en un conjunto de puntos")
    print("3. Campo magnético a lo largo del tiempo en un punto")
    print("4. Ploteo de los resultados")
    print("c. Exit. Salir")

    while True:
        opcion = input("Opcion: ")
        if opcion == '1':
            print("Campo magnetico en un punto")
            opciones1()
            break
        elif opcion == '2':
            print("Campo magnetico en un conjunto de puntos")
            opciones2()
            break
        elif opcion == '3':
            print("Campo magnético a lo largo del tiempo en un punto")
            opciones3()
            break
        elif opcion == '4':
            print("Ploteo de los resultados")
            opciones4()
            break
        elif opcion == 'c':
            print("Saliendo del programa")
            sys.exit()
        else:
            print("Opcion no valida")
            print("Seleccione una opcion:")
            print("\n")
            print("1. Campo magnetico en un punto")
            print("2. Campo magnetico en un conjunto de puntos")
            print("3. Campo magnético a lo largo del tiempo en un punto")

if __name__ =='__main__':
    bienvenida()
    opcion = menu()

