# This script allows the user to calculate the magnetic field at a given location
# on the Earth's surface. The user can specify the location using either geographic
# coordinates (latitude, longitude, and altitude) or geocentric coordinates (X, Y, Z).
# The user can also specify the date and time for which to calculate the magnetic
# field components. 

import numpy as np
import pandas as pd
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


    def get_txt(self):
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
                data = self.read_txt(filename)
            else:
                print('El archivo no existe')
                # Save the txt file
                self.write_txt(filename, txt)
                # We read the txt file
                data = self.read_txt(filename)

            self.data = data

            # Return the txt file
            return data
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

    def read_txt(self, filemane):
        '''This function reads the txt file'''
        # We open the file in read mode
        
        data = pd.read_csv(filemane, sep='\s+', skiprows=3, header=0, engine='python')

        # We print a message to the user
        print(f"El archivo {filemane} se ha leído correctamente.")
        print(data.head())

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

    def get_gh_data(self, n, m, year, coeff):
        '''This function returns the g and h coefficients for a given n, m, and year
        Parameters
        ----------
        n : int
            The order of the spherical harmonic function
        m : int
            The degree of the spherical harmonic function
        year : int
            The year for which to calculate the coefficients
        coeff : str
            The type of coefficient to return. Can be 'b' (both), 'g' or 'h'
        Returns
        -------
        gh : dict
            The value of the g or h coefficient
        '''

        # We filter the data by g and h coefficients

        index = []
        c = 0
        for ii in range(len(self.data)):
            if self.data['n'][ii] == n and self.data['m'][ii] == m:
                if coeff == self.data['g/h'][ii]:
                    index.append(ii)
                    break
                elif coeff == 'b':
                    index.append(ii)
                    c += 1
            if c == 2:
                break
        
        gh = {} # Dictionary to store the coefficients

        # check if the year is between the range of years in the data
        if float(self.data.columns[3]) <= float(year) <= float(self.data.columns[-2]):
            # Interpolate year data to get the coefficients
            for jj in index:
                if float(year) % 5 == 0:
                    gh[self.data['g/h'][jj]] = self.data[str(float(year))][jj]
                else:
                    year = float(year)
                    year1 = year - year % 5
                    year2 = year1 + 5
                    gh1 = self.data[str(year1)][jj]
                    gh2 = self.data[str(year2)][jj]
                    gh[self.data['g/h'][jj]] = (gh1 + (gh2 - gh1) * (year - year1) / 5)
        elif float(self.data.columns[-2]) < float(year) <= float(self.data.columns[-2]) + 5.0:
            # Extrapolate year data to get the coefficients
            for jj in index:
                year = float(year)
                year1 = self.data.columns[-2] # Last year in the data
                value_year1 = self.data[year1][jj] # Value of the coefficient for the last year in the data
                # The rate per year is the value of the last column
                rate = self.data.iloc[jj, -1] 

                gh[self.data['g/h'][jj]] = value_year1 + rate * (year - float(year1))

        return gh


    def Snm(self, n, m):
        '''This function calculates the Snm function
        Parameters
        ----------
        n : int
            The order of the spherical harmonic function
        m : int
            The degree of the spherical harmonic function
        Returns
        -------
        Snm : float
            The value of the Snm function
        '''
        
        if m==0 and n==0:
            Snm = 1

        pass

    def legendre_poly(nmax, theta):
        """
        Returns associated Legendre polynomials `P(n,m)` (Schmidt quasi-normalized)
        and the derivative :math:`dP(n,m)/d\\theta` evaluated at :math:`\\theta`.

        Parameters
        ----------
        nmax : int, positive
            Maximum degree of the spherical expansion.
        theta : ndarray, shape (...)
            Colatitude in degrees :math:`[0^\\circ, 180^\\circ]`
            of arbitrary shape.

        Returns
        -------
        Pnm : ndarray, shape (n, m, ...)
            Evaluated values and derivatives, grid shape is appended as trailing
            dimensions. `P(n,m)` := ``Pnm[n, m, ...]`` and `dP(n,m)` :=
            ``Pnm[m, n+1, ...]``

        """

        costh = np.cos(np.radians(theta))
        sinth = np.sqrt(1-costh**2)

        Pnm = np.zeros((nmax+1, nmax+2) + costh.shape)
        Pnm[0, 0] = 1  # is copied into trailing dimenions
        Pnm[1, 1] = sinth  # write theta into trailing dimenions via broadcasting

        rootn = np.sqrt(np.arange(2 * nmax**2 + 1))

        # Recursion relations after Langel "The Main Field" (1987),
        # eq. (27) and Table 2 (p. 256)
        for m in range(nmax):
            Pnm_tmp = rootn[m+m+1] * Pnm[m, m]
            Pnm[m+1, m] = costh * Pnm_tmp

            if m > 0:
                Pnm[m+1, m+1] = sinth*Pnm_tmp / rootn[m+m+2]

            for n in np.arange(m+2, nmax+1):
                d = n * n - m * m
                e = n + n - 1
                Pnm[n, m] = ((e * costh * Pnm[n-1, m] - rootn[d-e] * Pnm[n-2, m])
                            / rootn[d])

        # dP(n,m) = Pnm(m,n+1) is the derivative of P(n,m) vrt. theta
        Pnm[0, 2] = -Pnm[1, 1]
        Pnm[1, 2] = Pnm[1, 0]
        for n in range(2, nmax+1):
            Pnm[0, n+1] = -np.sqrt((n*n + n) / 2) * Pnm[n, 1]
            Pnm[1, n+1] = ((np.sqrt(2 * (n*n + n)) * Pnm[n, 0]
                        - np.sqrt((n*n + n - 2)) * Pnm[n, 2]) / 2)

            for m in np.arange(2, n):
                Pnm[m, n+1] = (0.5*(np.sqrt((n + m) * (n - m + 1)) * Pnm[n, m-1]
                            - np.sqrt((n + m + 1) * (n - m)) * Pnm[n, m+1]))

            Pnm[n, n+1] = np.sqrt(2 * n) * Pnm[n, n-1] / 2

        return Pnm

    def dipole(self, phi, theta, r):
        '''This function calculates the magnetic field due to a dipole
        Parameters
        ----------
        phi : float
            The co-elevation of the point where to calculate the magnetic field
        theta : float
            The  East longitude of the point where to calculate the magnetic field
        r : float
            The distance from the center of the Earth to the point where to calculate the magnetic field
        Returns
        -------
        B : float
            The magnetic field at the given location
        '''
        
        pass

            
if __name__ == '__main__':
    geomag = Geomat()
    data = geomag.get_txt()
    print(data)

    print(geomag.get_gh_data(10, 4, '2022.5', 'b'))
    # geomag.get_all_txt_from_url()


