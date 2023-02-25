# This script allows the user to calculate the magnetic field at a given location
# on the Earth's surface. The user can specify the location using either geographic
# coordinates (latitude, longitude, and altitude) or geocentric coordinates (X, Y, Z).
# The user can also specify the date and time for which to calculate the magnetic
# field components. 

import numpy as np
import requests
from bs4 import BeautifulSoup
import re
from pathlib import Path

class Geomat:
    '''This class contains the functions necessary to calculate the magnetic field'''

    def __init__(self):
        '''This function initializes the class'''
        # We define the constants
        self.a = 6378.137
        self.url = self.get_all_txt_from_url()
  
# Esta función extrae el archivo txt de la url indicada
# La fuente es: https://www.ngdc.noaa.gov/IAGA/vmod/coeffs/igrf13coeffs.txt


    def extraer_txt(self):
        '''This function extracts the txt file from the indicated url'''

        if self.url is None:
            'https://www.ngdc.noaa.gov/IAGA/vmod/coeffs/igrf13coeffs.txt'

        # Hacemos una petición GET a la url
        response = requests.get(self.url)
        # Comprobamos que la respuesta sea exitosa
        if response.status_code == 200:
            # Obtenemos el contenido de la respuesta como texto
            txt = response.text

            # We use the split method to separate the string by the character '/'
            parts = self.url.split('/')
            # We get the last part of the resulting list
            last_part = parts[-1]
            # We print the last part
            print(last_part)

            filename = 'data/' + last_part

            fileObj = Path(filename)
            
            if fileObj.is_file():
                print('El archivo ya existe')
                # We read the txt file
                txt = self.leer_txt(filename)
            else:
                print('El archivo no existe')
                # Save the txt file
                self.write_txt(filename, txt)

            # Return the txt file
            return txt
        else:
            # If the response is not successful, we return None
            return None

    def write_txt(self, filename, data):
        
        # We open the file in write mode
        file = open(filename, "w")
        # We write the data to the file
        file.write(data)
        # We close the file
        file.close()
        # We print a message to the user
        print(f"El archivo {filename} se ha guardado correctamente.")

    def leer_txt(self, filemane):
        '''This function reads the txt file'''
        # We open the file in read mode
        with open(filemane, 'r') as f:
            # We read the file
            data = f.read()
            # We return the data
            return data

    def get_all_txt_from_url(self, url="https://www.ncei.noaa.gov/products/international-geomagnetic-reference-field"):

        # Get the content of the web page
        response = requests.get(url)
        soup = BeautifulSoup(response.text, "html.parser")

        # Find all the links to .txt files
        txt_files = soup.find_all("a", href=re.compile("\.txt$"))

        # Find the substring in the content of the .txt files
        substring = "igrf"

        urls = []

        for txt_file in txt_files:
            txt_url = txt_file["href"] # Get the url of the .txt file
            urls.append(txt_url) # Add the url to the list

        print(urls) # Print the list of urls

        num_max = 0
        for cadena in urls:
            idx = cadena.index('igrf') + 4
            idx_f = idx
  
            for ii in cadena[idx:]:
                idx_f += 1
                if ii.isdigit() == False:
                    idx_f -= 1
                    break

            valor = int(cadena[idx:idx_f])
            if valor > num_max:
                num_max = valor

        print('igrf' + str(num_max))

        # Find the num_max between the urls

        for url_txt in urls:
            if 'igrf' + str(num_max) in url_txt:
                return url_txt
            else:
                return None
            

if __name__ == '__main__':
    geomag = Geomat()
    data = geomag.extraer_txt()
    print(data)
    # geomag.get_all_txt_from_url()

