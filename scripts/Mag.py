from src.geomag import Geomag
import sys
import numpy as np

geomag = Geomag()

def opciones1():
    '''Diferentes modelos de campo magnetico,
    dipolo, dipolo centrado, cuadripolo, octupolo, modelo completo'''

    print("Seleccione una opcion:")
    print("1. Dipolo")
    print("2. Dipolo centrado")
    print("3. Cuadripolo")
    print("4. Octupolo")
    print("5. Modelo completo")
    print("6. Salir")

    while True:
        try:
            opcion = input("Opcion: ")

            lat, long, h = latlongh()
            year = yearf()

            if opcion == '1':
                print("Dipolo")
                # Se necesita el radio de la tierra, la latitud y longitud
                val = geomag.dipole(h, lat, long, year)[1]
                break
            elif opcion == '2':
                print("Dipolo centrado")
                val = geomag.centered_dipole(h, lat, long, year)[1]
                break
            elif opcion == '3':
                print("Cuadripolo")
                val = geomag.quadrupole(h, lat, long, year)[1]
                break
            elif opcion == '4':
                print("Octupolo")
                val = geomag.octupole(h, lat, long, year)[1]
                break
            elif opcion == '5':
                print("Modelo completo")
                val = geomag.magnetic_field(h, lat, long, year, N=13)[1]
                break
            elif opcion == '6':
                print("Saliendo del programa")
                sys.exit()
            else:
                print("Opcion no valida")
                print("Seleccione una opcion:")
                print("1. Dipolo")
                print("2. Dipolo centrado")
                print("3. Cuadripolo")
                print("4. Octupolo")
                print("5. Modelo completo")

            print("Guardando datos en results/magPoint.txt")
            geomag.printresults(val, savedata=True)
        except ValueError:
            print("Introduzca un valor válido")


def opciones2():
    '''Campo magnetico en un conjunto de puntos'''

    print("Seleccione una opcion:")
    # Crear array de puntos
    print("1. Crear array de puntos")
    # Leer array de puntos
    print("2. Leer array de puntos")
    # Calcular campo magnetico en los puntos
    print("3. Introductir puntos manualmente")

    while True:
        try:
            opcion = input("Opcion: ")
            if opcion == '1':
                print("Crear array de puntos")
                lat, long, h, year = opcion211()

            elif opcion == '2':
                print("Leer array de puntos")
                lat, long, h, year = opcion212()

            elif opcion == '3':
                print("Introductir puntos manualmente")
                lat = []
                long = []
                h = []
                while True:
                    try:
                        n = int(input("Numero de puntos: "))
                        for i in range(n):
                            lati, longi, hi, yeari = latlongh()
                            lat.append(lati)
                            long.append(longi)
                            h.append(hi)
                        year = float(input("Año: [decimal]"))

                    except ValueError:
                        print("Introduzca un valor válido")
            else:
                print("Opcion no valida")
                print("Seleccione una opcion:")
                print("1. Crear array de puntos")
                print("2. Leer array de puntos")
                print("3. Introductir puntos manualmente")
        except ValueError:
            print("Introduzca un valor válido")

        geomag.arrayofpoints(h, lat, long, year, N=13, savedata=True)[1]
        print("Datos guardados en results/geomag.csv")
        break

def opcion211():
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
            year = float(input("Año: [decimal]"))
            break
        except ValueError:
            print("Introduzca un valor válido")

    # Crear array de puntos
    lat = np.linspace(latmin, latmax, n)
    long = np.linspace(longmin, longmax, n)
    h = np.linspace(hmin, hmax, n)

    return lat, long, h, year

def opcion212():
    '''Esta función lee un array de puntos'''
    print("Leer array de puntos")
    print("Los datos deben estar en la carpeta data")
    print("Los datos deben estar en formato .txt")
    print("Los datos deben estar en el siguiente orden: lat, long, h, year")
    print("Los datos deben estar separados por espacios")
    print("Los datos deben estar en columnas")

    # Leer array de puntos
    lat = np.loadtxt('data/lat.txt')
    long = np.loadtxt('data/long.txt')
    h = np.loadtxt('data/h.txt')
    year = np.loadtxt('data/year.txt')

    return lat, long, h, year

def opciones3():
    '''Campo magnético a lo largo del tiempo en un punto'''
    while True:
        try:
            # Preguntar al usuario el punto
            lat, long, h = latlongh()

            tmin, tmax, n = timerange()
            
            # Crear array de tiempos
            t = np.linspace(tmin, tmax, n)

            # Calcular campo magnético
            geomag.arrayDatesAtLocation([h, lat, long] , t, N=13, savedata=True)
            break
        except ValueError:
            print("Introduzca un valor válido")

def opciones4():

    while True:
        try:
            lat, long, h, = latlongh()

            tmin, tmax, n = timerange()
            # Crear array de tiempos
            
            geomag.plotMagneticField([h, lat, long], tmin, N=13, modelos=[1, 1, 1, 1,1], timeRange= [tmin, tmax])
            break
        except ValueError:
            print("Introduzca un valor válido")
# Función que pide al usuario en grados la latitud, longitud y altura

def opciones5():
    '''Cambio de coordenadas'''

    print("Seleccione una opcion:")
    print("1. Geodésicas a geocentricas")
    print("2. Geocentricas a geodésicas")
    print("3. XYZ a DHIF")
    print("4. Coordenadas XYZ a inerciales")
    print("5. Coordenadas inerciales a orbitales")
    print("6. Coordenadas orbitales a ejes cuerpo")
    print("c para salir")

    while True:
        try:
            opcion = input("Opcion: ")

            if opcion == '1':
                print("Geodésicas a geocentricas")
                try:
                    h = float(input("Altura geodética: [km]"))
                    gdcolat = float(input("Colatitud geodésica: [grados]"))
                except ValueError:
                    print("Introduzca un valor válido")
                geomag.gg_to_geo(h, gdcolat)
                break
            elif opcion == '2':
                print("Geocentricas a geodésicas")
                try:
                    h = float(input("Altura geocéntrica: [km]"))
                    gdcolat = float(input("Colatitud geocéntrica: [grados]"))
                except ValueError:
                    print("Introduzca un valor válido")
                geomag.geo_to_gg(h, gdcolat)
                break
            elif opcion == '3':
                print("XYZ a DHIF")
                try:
                    x = float(input("Bx: [nT]"))
                    y = float(input("By: [nT]"))
                    z = float(input("Bz: [nT]"))
                except ValueError:
                    print("Introduzca un valor válido")
                geomag.xyz2dhif(x, y, z)
                break
            elif opcion == '4':
                print("Coordenadas XYZ a inerciales")
                try:
                    Btheta = float(input("Btheta : [nT]"))
                    Bphi = float(input("Bphi: [nT]"))
                    Br = float(input("Br: [nT]"))
                    lat = float(input("Latitud: [grados]"))
                    long = float(input("Longitud: [grados]"))
                    thetag = float(input("Tiempo celestial en Greenxich: [grados]"))
                except ValueError:
                    print("Introduzca un valor válido")
                geomag.transformation2inertial([Br, Btheta, Bphi],lat, long, thetag)
                break
            elif opcion == '5':
                print("Coordenadas inerciales a orbitales")
                try:
                    Btheta = float(input("Btheta : [nT]"))
                    Bphi = float(input("Bphi: [nT]"))
                    Br = float(input("Br: [nT]"))
                    lat = float(input("Latitud: [grados]"))
                    long = float(input("Longitud: [grados]"))
                    thetag = float(input("Tiempo celestial en Greenxich: [grados]"))
                except ValueError:
                    print("Introduzca un valor válido")
                geomag.transformation2orb([Br, Btheta, Bphi], lat, long, thetag)
                break
            elif opcion == '6':
                print("Coordenadas orbitales a ejes cuerpo")
                try:
                    Bxorb = float(input("Bx orb: [nT]: "))
                    Byorb = float(input("By orb: [nT]: "))
                    Bzorb = float(input("Bz orb: [nT]: "))
                    roll = float(input("Roll: [grados]: "))
                    pitch = float(input("Pitch: [grados]: "))
                    yaw = float(input("Yaw: [grados]: "))
                except ValueError:
                    print("Introduzca un valor válido")
                geomag.orb2body([Bxorb, Byorb, Bzorb],[roll, pitch, yaw])
                geomag.linearOrb2body([Bxorb, Byorb, Bzorb],[roll, pitch, yaw])
            elif opcion == 'c':
                print("Saliendo del programa")
                sys.exit()
            else:
                print("Opcion no valida")
                print("Seleccione una opcion:")
                print("1. Geodésicas a geocentricas")
                print("2. Geocentricas a geodésicas")
                print("3. XYZ a DHIF")
                print("4. Coordenadas XYZ a inerciales")
                print("5. Coordenadas inerciales a orbitales")
                print("6. Coordenadas orbitales a ejes cuerpo")
                print("c para salir")
        except ValueError:
            print("Introduzca un valor válido")
            print("Opcion no valida")
            print("Seleccione una opcion:")
            print("1. Geodésicas a geocentricas")
            print("2. Geocentricas a geodésicas")
            print("3. XYZ a DHIF")
            print("4. Coordenadas XYZ a inerciales")
            print("5. Coordenadas inerciales a orbitales")
            print("6. Coordenadas orbitales a ejes cuerpo")
            print("c para salir")



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

def yearf():
    while True:
        try:
            year = float(input("Año: [decimal]"))
            break
        except ValueError:
            print("Introduzca un valor válido")

    return year

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