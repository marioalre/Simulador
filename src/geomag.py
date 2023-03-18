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
        self.Re = 6378.137
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


    def Snm(self, nm_max=13):
        '''This function calculates the Snm function
        Parameters
        ----------
        nm_max : int
            The maximum value of n and m
            n is the order of the spherical harmonic function
            m is the degree of the spherical harmonic function
        Returns
        -------
        Snm : float
            The value of the Snm function
        '''
        
        Snm = np.zeros((nm_max + 1, nm_max + 1))

        for m in range(nm_max + 1): # m is the degree of the spherical harmonic function
            for n in range(nm_max + 1): # n is the order of the spherical harmonic function
                if m==0 and n==0:
                    Snm[n, m] = 1
                elif m==0 and n>=1:
                    Snm[n, m] = np.sqrt((2*n - 1)/(n)) * Snm[n-1, m]
                elif m>=1 and n>=m:
                    factor = np.sqrt(2)
                    if m==1:
                        Snm[n, m] = factor * np.sqrt((n - m + 1)/(n + m)) * Snm[n, m-1]
                    else:
                        Snm[n, m] = np.sqrt((n - m + 1)/(n + m)) * Snm[n, m-1] 
                
        return Snm
    
    def gauss_norm_ass_leg_poly(self, nm_max, theta):
        '''This function calculates the associated Legendre polynomials with the normalization factor
        Gauss's normalization factor is used. 
        Recurrence relation for the associated Legendre polynomials is used.
        Parameters
        ----------
        nm_max : int
            The maximum value of n and m
            n is the order of the spherical harmonic function
            m is the degree of the spherical harmonic function
        theta : float
            The latitude in radians
        Returns
        -------
        Pnm : array 
            The value of the associated Legendre polynomials with the normalization factor
        dPnm : array
            The derivative of the associated Legendre polynomials with the normalization factor
        '''

        Pnm = np.zeros((nm_max + 1, nm_max + 1))
        dPnm = np.zeros((nm_max + 1, nm_max + 1))

        for m in range(nm_max + 1): # m is the degree of the spherical harmonic function
            for n in range(nm_max + 1): # n is the order of the spherical harmonic function
                if m==0 and n==0:
                    Pnm[n, m] = 1
                    dPnm[n, m] = 0
                elif m==n:
                    Pnm[n, m] = np.sin(theta) * Pnm[n-1, m-1]
                    dPnm[n, m] = np.cos(theta) * Pnm[n-1, m-1] + np.sin(theta) * dPnm[n-1, m-1]
                elif m>=0 and n>=1 and n>m:
                    if n==1:
                        Knm = 0   
                    elif n > 1:
                        Knm = ((n - 1)^2 - m^2)/((2*n-1)*(2*n - 3))
                    
                    Pnm[n, m] = np.cos(theta) * Pnm[n-1, m] - Knm * Pnm[n-2, m] #Posible error cuando n=m=1

                    dPnm[n, m] = -np.sin(theta) * Pnm[n-1, m] + np.cos(theta) * dPnm[n-1, m] - Knm * dPnm[n-2, m]   

        return Pnm, dPnm

    def get_gh_norm(self, n, m, year, coeff):
        '''This function returns the normalization factor for a given n and m
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

        Snm = self.Snm()

        data = self.get_gh_data(n, m, year, coeff)

        if coeff == 'b':
            data['g'] = data['g'] * Snm[n, m]
            data['h'] = data['h'] * Snm[n, m]
        else:
            data[coeff] = data[coeff] * Snm[n, m]

        return data

    def dipole(self, phi, theta, r, year=2020):
        '''This function calculates the magnetic field due to a dipole
        Parameters
        ----------
        phi : float
            The co-elevation of the point where to calculate the magnetic field (in degrees)
        theta : float
            The  East longitude of the point where to calculate the magnetic field (in degrees)
        r : float
            The distance from the center of the Earth to the point where to calculate the magnetic field (in km)
        Returns
        -------
        B : float
            The magnetic field at the given location
        '''
        gh = self.get_gh_norm(n = 1, m = 1, year = year, coeff = 'b')
        h11 = gh['h']
        g11 = gh['g']
        g10 = self.get_gh_norm(n = 1, m = 0, year = year, coeff = 'g')['g']

        # Convert to radians
        phi = np.radians(phi)
        theta = np.radians(theta)

        # Verify the limits of the input
        if r < 0:
            raise ValueError('The distance must be positive')
        if r < self.Re:
            raise ValueError('The distance must be higher than the radius of the Earth')
        if phi <= -np.pi/2 or phi >= np.pi/2:
            raise ValueError('The co-elevation must be between -90 and 90')
        if theta <= -np.pi or theta >= np.pi:
            raise ValueError('The longitude must be between -180 and 180')
        
        # Calculate the magnetic field

        Br = 2 * (r/self.Re)**3 * (g10 * np.cos(theta) + (g11 * np.cos(phi) + h11 * np.sin(phi)) * np.sin(theta))
        Btheta = (r/self.Re)**3 * (g10 * np.sin(theta) - (g11 * np.cos(phi) + h11 * np.sin(phi)) * np.cos(theta))
        Bphi = (r/self.Re)**3 * (-h11 * np.cos(phi) + g11 * np.sin(phi))

        Bvector = np.array([Br, Btheta, Bphi])
        Bmodule = np.linalg.norm(Bvector)
        
        return Bmodule, Bvector
    
    def centered_dipole(self, phi, theta, r, year=2020):
        '''This function calculates the magnetic field due to a dipole
        Parameters
        ----------
        phi : float
            The co-elevation of the point where to calculate the magnetic field (in degrees)
        theta : float
            The  East longitude of the point where to calculate the magnetic field (in degrees)
        r : float
            The distance from the center of the Earth to the point where to calculate the magnetic field (in km)
        Returns
        -------
        B : float
            The magnetic field at the given location
        '''

        g10 = self.get_gh_norm(n = 1, m = 0, year = year, coeff = 'g')['g']

        # Convert to radians
        phi = np.radians(phi)
        theta = np.radians(theta)

        # Verify the limits of the input
        if r < 0:
            raise ValueError('The distance must be positive')
        if r < self.Re:
            raise ValueError('The distance must be higher than the radius of the Earth')
        if phi <= -np.pi/2 or phi >= np.pi/2:
            raise ValueError('The co-elevation must be between -90 and 90')
        if theta <= -np.pi or theta >= np.pi:
            raise ValueError('The longitude must be between -180 and 180')
        
        # Calculate the magnetic field

        Br = 2 * (r/self.Re)**3 * g10 * np.cos(theta)
        Btheta = (r/self.Re)**3 * g10 * np.sin(theta)
        Bphi = 0

        Bvector = np.array([Br, Btheta, Bphi])
        Bmodule = np.linalg.norm(Bvector)
        
        return Bmodule, Bvector
    
    def quadrupole(self, phi, theta, r, year=2020):
        '''This function calculates the magnetic field due to a dipole
        Parameters
        ----------
        phi : float
            The co-elevation of the point where to calculate the magnetic field (in degrees)
        theta : float
            The  East longitude of the point where to calculate the magnetic field (in degrees)
        r : float
            The distance from the center of the Earth to the point where to calculate the magnetic field (in km)
        Returns
        -------
        B : float
            The magnetic field at the given location
        '''
        gh = self.get_gh_norm(n = 2, m = 2, year = year, coeff = 'b')
        h22 = gh['h']
        g22 = gh['g']
        gh = self.get_gh_norm(n = 2, m = 1, year = year, coeff = 'b')
        h21 = gh['h']
        g21 = gh['g']
        g20 = self.get_gh_norm(n = 2, m = 0, year = year, coeff = 'g')['g']


        # Convert to radians
        phi = np.radians(phi)
        theta = np.radians(theta)

        # Verify the limits of the input
        if r < 0:
            raise ValueError('The distance must be positive')
        if r < self.Re:
            raise ValueError('The distance must be higher than the radius of the Earth')
        if phi <= -np.pi/2 or phi >= np.pi/2:
            raise ValueError('The co-elevation must be between -90 and 90')
        if theta <= -np.pi or theta >= np.pi:
            raise ValueError('The longitude must be between -180 and 180')
        
        # Calculate the magnetic field

        Bm , B = self.dipole(phi, theta, r, year)

        Br_dip = B[0]
        Btheta_dip = B[1]
        Bphi_dip = B[2]

        Br = Br_dip + 3 * (self.Re/r)**4 * (0.5*g20*(np.cos(2*theta) + 1/3) +0.5*(g21 * np.cos(phi) + h21 * np.sin(phi)) * np.sin(2*theta) + \
            0.5*(g22 * np.cos(2*phi) + h22 * np.sin(2*phi)) * (1- np.cos(2*theta)))
        
        Btheta = Btheta_dip + (self.Re/r)**4 * (g20*np.sin(2*theta) + (g21*np.cos(phi) + h21*np.sin(phi))*np.cos(2*theta) - \
            (g22*np.cos(2*phi) + h22*np.sin(2*phi))*np.sin(2*theta))
        
        Bphi = Bphi_dip + (self.Re/r)**4 * ((g21*np.sin(phi) - h21*np.cos(phi)) * np.cos(theta) - \
            2 * (g22*np.sin(2*phi) - h22*np.cos(2*phi)) * np.sin(theta))

        Bvector = np.array([Br, Btheta, Bphi])
        Bmodule = np.linalg.norm(Bvector)

        return Bmodule, Bvector
            
if __name__ == '__main__':
    geomag = Geomat()
    data = geomag.get_txt()
    print(data)
    
    print(geomag.Snm())
    print(geomag.get_gh_data(1, 1, '2015.0', 'b'))
    print(geomag.dipole(0, 0, 7000))
    print(geomag.centered_dipole(0, 0, 7000))
    print(geomag.quadrupole(0, 0, 7000))

    # geomag.gauss_norm_ass_leg_poly(2, np.pi/4)


