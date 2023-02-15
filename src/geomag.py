# This script allows the user to calculate the magnetic field at a given location
# on the Earth's surface. The user can specify the location using either geographic
# coordinates (latitude, longitude, and altitude) or geocentric coordinates (X, Y, Z).
# The user can also specify the date and time for which to calculate the magnetic
# field components. 

import numpy as np
import urllib.request
import requests

class Geomat:
    '''This class contains the functions necessary to calculate the magnetic field'''

    def __init__(self):
        '''This function initializes the class'''
        # We define the constants
        self.a = 6378.137
  
# Esta función extrae el archivo txt de la url indicada
# La fuente es: https://www.ngdc.noaa.gov/IAGA/vmod/coeffs/igrf13coeffs.txt


    def extraer_txt(self):
        '''This function extracts the txt file from the indicated url'''
        
        url = 'https://www.ngdc.noaa.gov/IAGA/vmod/coeffs/igrf13coeffs.txt'
        # Hacemos una petición GET a la url
        response = requests.get(url)
        # Comprobamos que la respuesta sea exitosa
        if response.status_code == 200:
            # Obtenemos el contenido de la respuesta como texto
            txt = response.text
            print(txt)

            # We use the split method to separate the string by the character '/'
            parts = url.split('/')
            # We get the last part of the resulting list
            last_part = parts[-1]
            # We print the last part
            print(last_part)

            filename = 'data/' + last_part
            # Abrimos el archivo en modo escritura
            file = open(filename, "w")
            # Escribimos el texto en el archivo
            file.write(txt)
            # Cerramos el archivo
            file.close()
            # Imprimimos un mensaje de confirmación
            print(f"El archivo {filename} se ha guardado correctamente.")

            # Devolvemos el archivo txt
            return txt
        else:
            # Si la respuesta no es exitosa, devolvemos None
            return None

    def leer_txt(self):
        '''This function reads the txt file'''
        # Abrimos el archivo txt
        with open('igrf13coeffs.txt', 'r') as f:
            # Leemos el archivo
            data = f.read()
            # Devolvemos el archivo
            return data

if __name__ == '__main__':
    geomag = Geomat()
    data = geomag.extraer_txt()

