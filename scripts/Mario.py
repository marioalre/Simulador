# Pantalla de bienvenida del simulador con un menu para las diferentes clases que hay en src
# Autor: Mario 

import os
import sys
from scripts.Mag1pto import opciones1, opciones2, opciones3, opciones4
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
    print("Bienvenido al simulador de redes")
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
    print("Bienvenido al simulador de redes")
    print("Seleccione una opcion:")
    print("1. Campo magnetico")
    print("2. Campo goeopotencial")
    print("3. Propagadores orbitas")
    print("4. Cambio de coordenadas")
    print("5. Determinación preliminar de la órbita")
    print("6. Salir")

    while True:
        opcion = input("Opcion: ")
        if opcion == '1':
            print("Campo magnetico")
            menu_campo_magnetico()
            break
        elif opcion == '2':
            print("Campo goeopotencial")
            break
        elif opcion == '3':
            print("Propagadores orbitas")
            break
        elif opcion == '4':
            print("Cambio de coordenadas")
            break
        elif opcion == '5':
            print("Determinación preliminar de la órbita")
            break
        elif opcion == '6':
            print("Saliendo del programa")
            sys.exit()
        else:
            print("Opcion no valida")
            print("Seleccione una opcion:")
            print("1. Campo magnetico")
            print("2. Campo goeopotencial")
            print("3. Propagadores orbitas")
            print("4. Cambio de coordenadas")
            print("5. Determinación preliminar de la órbita")

    return opcion

def menu_campo_magnetico():
    '''Menu del campo magnetico'''

    print("Bienvenido al simulador del campo magnetico")
    print("Se puede calcular el campo magnetico en un punto, en un conjunto de puntos o a lo largo del tiempo en un punto")
    print("Se utiliza el modelo IGRF-13")

    print("Seleccione una opcion:")

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
            print("1. Campo magnetico en un punto")
            print("2. Campo magnetico en un conjunto de puntos")
            print("3. Campo magnético a lo largo del tiempo en un punto")




if __name__ =='__main__':
    bienvenida()
    opcion = menu()

