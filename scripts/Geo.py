from src.geopot import Geopot
import sys
import numpy as np

geopot = Geopot()

def opciones1a():
    '''Calcular el campo geopotencial en un punto'''
    while True:
        try:
            # Preguntar al usuario el punto
            lat, long, h = latlongh()

            # preguntar el orden, este tiene que ser de 1 a 360
            n = int(input("Orden del campo geopotencial: "))
            if n < 1 or n > 360:
                print("El orden debe ser de 1 a 360")
            else:
                # Calcular campo geopotencial
                geopot.gravitational_potential(h, lat, long, n)
                break

        except ValueError:
            print("Introduzca un valor válido")

def opciones2a():
    '''Calcular el campo geopotencial en un conjunto de puntos'''

    print("Seleccione una opcion:")
    print("\n")
    # Crear array de puntos
    print("1. Crear array de puntos")
    # Leer array de puntos
    print("2. Leer array de puntos")
    # Calcular campo geopotencial en los puntos
    print("3. Introductir puntos manualmente")

    while True:
        try:
            opcion = input("Opcion: ")
            if opcion == '1':
                print("Crear array de puntos")
                lat, long, h = opcion211a()

            elif opcion == '2':
                print("Leer array de puntos")
                lat, long, h = opcion212a()

            elif opcion == '3':
                print("Introductir puntos manualmente")
                lat = []
                long = []
                h = []
                while True:
                    try:
                        n = int(input("Numero de puntos: "))
                        for i in range(n):
                            lati, longi, hi = latlongh()
                            lat.append(lati)
                            long.append(longi)
                            h.append(hi)

                    except ValueError:
                        print("Introduzca un valor válido")

            else:
                print("Opcion no valida")
                print("Seleccione una opcion:")
                print("1. Crear array de puntos")
                print("2. Leer array de puntos")
                print("3. Introductir puntos manualmente")

            # preguntar el orden, este tiene que ser de 1 a 360
            n = int(input("Orden del campo geopotencial: "))
            if n < 1 or n > 360:
                print("El orden debe ser de 1 a 360")
            else:
                # Calcular campo geopotencial
                geopot.arraylatlong(lat, long, h, n, savedata=True)
                break
        
        except ValueError:
            print("Introduzca un valor válido")

def opciones3a():
    '''Plotear el campo geopotencial en un conjunto de puntos'''

    print("Obtencion de una representacion grafica del campo geopotencial")
    print("Seleccione una opcion:")
    print("\n")
    print("1. grafica 2D")
    print('2. grafica 3D')

    while True:
        try:
            opcion = input("Opcion: ")
            if opcion == '1':
                print("grafica 2D") 
                type_ = '2D' 
                break
            elif opcion == '2':
                print("grafica 3D")
                type_ = '3D'
                break
            else:
                print("Opcion no valida")
                print("Seleccione una opcion:")
                print("1. grafica 2D")
                print('2. grafica 3D')
        except ValueError:
            print("Introduzca un valor válido")

    # Preguntar el orden
    while True:
        try:
            n = int(input("Orden del campo geopotencial: "))
            if n > 1 or n < 360:
                break
            else:
                print("El orden debe ser de 1 a 360")
        except ValueError:
            print("Introduzca un valor válido")



    print("Elija la resolución del mallado")
    print("1. Alta")
    print("2. Media")
    print("3. Baja")
    while True:
        try:
            res = input("Opcion: ")
            if res == '1':
                res = 20
                break
            elif res == '2':
                res = 40
                break
            elif res == '3':
                res = 80
                break
            else:
                print("Opcion no valida")
                print("Seleccione una opcion:")
                print("1. Alta")
                print("2. Media")
                print("3. Baja")
        except ValueError:
            print("Introduzca un valor válido")

    # Preguntar el potencial o gravedad
    print("Seleccione el campo a representar")
    print("1. Potencial")
    print("2. Gravedad")
    while True:
        try:
            campo = input("Opcion: ")
            if campo == '1':
                campo = 'potential'
                break
            elif campo == '2':
                campo = 'gravity'
                break
            else:
                print("Opcion no valida")
                print("Seleccione una opcion:")
                print("1. Potencial")
                print("2. Gravedad")
        except ValueError:
            print("Introduzca un valor válido")

    d1, d2, d3, = geopot.calculate(res, n, campo)
    data = [d1, d2, d3]

    geopot.plot_potential(data, dtype=type_)


def latlongh():
    '''Función que pide al usuario en grados la latitud, longitud y altura y año'''

    while True:
        try:
            lat = float(input("Latitud: [grados]"))
            long = float(input("Longitud: [grados]"))
            h = float(input("Altura (desde el centro de la Tierra): [km]"))
            break
        except ValueError:
            print("Introduzca un valor válido")
    return lat, long, h


def timerange():
    # Preguntar al usuario el tiempo
    while True:
        try:
            tmin = float(input("Tiempo mínimo: [año] "))
            tmax = float(input("Tiempo máximo: [año] "))
            n = int(input("Número de puntos: "))
            if tmin > tmax:
                print("El tiempo mínimo debe ser menor que el tiempo máximo")
            else:
                break
        except ValueError:
            print("Introduzca un valor válido")
    return tmin, tmax, n

def opcion212a():
    '''Esta función lee un array de puntos'''
    print("Leer array de puntos")
    print("Los datos deben estar en la carpeta data")
    print("Los datos deben estar en formato .txt")
    print("Los datos deben estar en el siguiente orden: lat.txt, long.txt, h.txt")
    print("Los datos deben en la carpeta data")

    # Leer array de puntos
    lat = np.loadtxt('data/lat.txt')
    long = np.loadtxt('data/long.txt')
    h = np.loadtxt('data/h.txt')

    return lat, long, h

def opcion211a():
    '''Esta función crea un array de puntos'''

    # preguntar al usuario el numero de puntos
    while True:
        try:
            n = int(input("Numero de puntos: "))
            break
        except ValueError:
            print("Introduzca un valor válido")
    # Preguntar los limites de latitud, longitud y altura
    while True:
        try:
            latmin = float(input("Latitud mínima: [grados] "))
            latmax = float(input("Latitud máxima: [grados] "))
            longmin = float(input("Longitud mínima: [grados] "))
            longmax = float(input("Longitud máxima: [grados] "))
            hmin = float(input("Altura mínima: [km] "))
            hmax = float(input("Altura máxima: [km] "))
            break
        except ValueError:
            print("Introduzca un valor válido")

    # Crear array de puntos
    lat = np.linspace(latmin, latmax, n)
    long = np.linspace(longmin, longmax, n)
    h = np.linspace(hmin, hmax, n)

    return lat, long, h