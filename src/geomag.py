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
import matplotlib.pyplot as plt
import os
import pandas as pd

class Geomag:
    '''This class contains the functions necessary to calculate the magnetic field'''

    def __init__(self):
        '''This function initializes the class'''
        # We define the constants
        self.Re = 6371.2
        self.url = self.get_all_txt_from_url()
        self.data = self.get_txt()
  
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
        '''This function reads the txt file
        VERIFIED
        '''
        # We open the file in read mode
        
        data = pd.read_csv(filemane, sep='\s+', skiprows=3, header=0, engine='python')

        # We print a message to the user
        print(f"El archivo {filemane} se ha leído correctamente.")
        print(data.head())

        return data

    def get_all_txt_from_url(self, url="https://www.ncei.noaa.gov/products/international-geomagnetic-reference-field"):
        '''
        This function returns the url of the .txt file with the most recent data
        VERIFIED
        '''

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
        VERIFIED
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
        
        self.max_mn = self.data['n'].max()

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
        VERIFIED
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
                    Snm[n, m] = (2*n - 1)/(n) * Snm[n-1, m]
                elif m>=1 and n>=m:
                    factor = np.sqrt(2)
                    if m==1:
                        Snm[n, m] = factor * np.sqrt((n - m + 1)/(n + m)) * Snm[n, m-1]
                    else:
                        Snm[n, m] = np.sqrt((n - m + 1)/(n + m)) * Snm[n, m-1] 
                
        return Snm
    
    def legendre_poly(self, nmax, theta):
        """
        Returns associated Legendre polynomials `P(n,m)` (Schmidt quasi-normalized)
        and the derivative :math:`dP(n,m)/d\\theta` evaluated at :math:`\\theta`.
        Parameters
        ----------
        nmax : int, positive
            Maximum degree of the spherical expansion.
        theta : ndarray, shape (...)
            Colatitude in radians :math:`[0^\\circ, 180^\\circ]`
            of arbitrary shape.
        Returns
        -------
        Pnm : ndarray, shape (n, m, ...)
            Evaluated values and derivatives, grid shape is appended as trailing
            dimensions. `P(n,m)` := ``Pnm[n, m, ...]`` and `dP(n,m)` :=
            ``Pnm[m, n+1, ...]``
        References
        ----------
        Based on Equations 26-29 and Table 2 in:
        Langel, R. A., "Geomagnetism - The main field", Academic Press, 1987,
        chapter 4
        """

        costh = np.cos(theta)
        sinth = np.sqrt(1-costh**2)

        Pnm = np.zeros((nmax+1, nmax+2) + costh.shape)
        Pnm[0, 0] = 1.  # is copied into trailing dimensions
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
    
    def gauss_norm_ass_leg_poly(self, nm_max, theta):
        '''This function calculates the associated Legendre polynomials with the normalization factor
        Gauss's normalization factor is used. 
        Recurrence relation for the associated Legendre polynomials is used.
        VERIFIED
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

        Snm = self.Snm(nm_max)

        Pnm[0, 0] = 1
        
        for n in range(1, nm_max +1): # m is the degree of the spherical harmonic function
            for m in range(0, min([n + 1, nm_max + 1])): # n is the order of the spherical harmonic function
                # do the legendre polynomials and derivatives
                if m==n:
                    Pnm[n, m] = np.sin(theta) * Pnm[n-1, m-1]
                    dPnm[n, m] = np.cos(theta) * Pnm[n-1, m-1] + np.sin(theta) * dPnm[n-1, m-1]
                else:
                    if n==1:
                        Knm = 0
                        Pnm[n, m] = np.cos(theta) * Pnm[n-1, m]
                        dPnm[n, m] = -np.sin(theta) * Pnm[n-1, m] + np.cos(theta) * dPnm[n-1, m]
                    elif n > 1:
                        Knm = ((n - 1)**2 - m**2)/((2*n-1)*(2*n - 3))                   
                        Pnm[n, m] = np.cos(theta) * Pnm[n-1, m] - Knm * Pnm[n-2, m] #Posible error cuando n=m=1
                        dPnm[n, m] = -np.sin(theta) * Pnm[n-1, m] + np.cos(theta) * dPnm[n-1, m] - Knm * dPnm[n-2, m]   
        '''     
        for n in range(1, nm_max +1): # m is the degree of the spherical harmonic function
            for m in range(0, min([n + 1, nm_max + 1])): # n is the order of the spherical harmonic function
                Pnm[n, m] = Pnm[n, m] * Snm[n, m]
                dPnm[n, m] = dPnm[n, m] * Snm[n, m]
        ''' 
        return Pnm, dPnm
    

    def get_gh_norm(self, n, m, year, coeff):
        '''This function returns the normalization factor for a given n and m
        VERIFIED
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
            if m != 0:
                data['h'] = data['h'] * Snm[n, m]
        else:
            data[coeff] = data[coeff] * Snm[n, m]

        return data

    def dipole(self, r, theta, phi, year=2020):
        '''This function calculates the magnetic field due to a dipole
        VERIFIED
        Parameters
        ----------
        phi : float
            Longitude measured in degrees positive east from Greenwich (in degrees)
        theta : float
            The latitude measured in degrees positive from equator (in degrees)
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
        theta = np.radians(90 - theta)

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

        Br = 2 * (self.Re/r)**3 * (g10 * np.cos(theta) + (g11 * np.cos(phi) + h11 * np.sin(phi)) * np.sin(theta))
        Btheta = (self.Re/r)**3 * (g10 * np.sin(theta) - (g11 * np.cos(phi) + h11 * np.sin(phi)) * np.cos(theta))
        Bphi = (self.Re/r)**3 * (-h11 * np.cos(phi) + g11 * np.sin(phi))

        Bvector = np.array([Br, Btheta, Bphi])
        Bmodule = np.linalg.norm(Bvector)
        
        return Bmodule, Bvector
    
    def centered_dipole(self, r, theta, phi, year=2020):
        '''This function calculates the magnetic field due to a dipole
        Parameters
        ----------
        phi : float
            Longitude measured in degrees positive east from Greenwich (in degrees)
        theta : float
            The latitude measured in degrees positive from equator (in degrees)
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
        theta = np.radians(90 -theta)

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

        Br = 2 * (self.Re/r)**3 * g10 * np.cos(theta)
        Btheta = (self.Re/r)**3 * g10 * np.sin(theta)
        Bphi = 0

        Bvector = np.array([Br, Btheta, Bphi])
        Bmodule = np.linalg.norm(Bvector)
        
        return Bmodule, Bvector
    
    def quadrupole(self, r, theta, phi, year=2020):
        '''This function calculates the magnetic field quadrupole
        Parameters
        ----------
        phi : float
            Longitude measured in degrees positive east from Greenwich (in degrees)
        theta : float
            The latitude measured in degrees positive from equator (in degrees)
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

        
        # Calculate the magnetic field

        Bm , B = self.dipole( r, theta, phi, year)

        # Convert to radians
        phi = np.radians(phi)
        theta = np.radians(90 - theta)

        # Verify the limits of the input
        if r < 0:
            raise ValueError('The distance must be positive')
        if r < self.Re:
            raise ValueError('The distance must be higher than the radius of the Earth')
        if phi <= -np.pi/2 or phi >= np.pi/2:
            raise ValueError('The co-elevation must be between -90 and 90')
        if theta <= -np.pi or theta >= np.pi:
            raise ValueError('The longitude must be between -180 and 180')

        Br_dip = B[0]
        Btheta_dip = B[1]
        Bphi_dip = B[2]

        Br = Br_dip + 3 * (self.Re/r)**4 * (0.5*g20*(np.cos(2*theta) + 1/3) +0.5*(g21 * np.cos(phi) + h21 * np.sin(phi)) * np.sin(2*theta) + \
            0.5*(g22 * np.cos(2*phi) + h22 * np.sin(2*phi)) * (1- np.cos(2*theta)))
        
        Btheta = Btheta_dip + (self.Re/r)**4 * (g20*np.sin(2*theta) - (g21*np.cos(phi) + h21*np.sin(phi))*np.cos(2*theta) - \
            (g22*np.cos(2*phi) + h22*np.sin(2*phi))*np.sin(2*theta))
        
        Bphi = Bphi_dip + (self.Re/r)**4 * ((g21*np.sin(phi) - h21*np.cos(phi)) * np.cos(theta) + \
            2 * (g22*np.sin(2*phi) - h22*np.cos(2*phi)) * np.sin(theta))

        Bvector = np.array([Br, Btheta, Bphi])
        Bmodule = np.linalg.norm(Bvector)

        return Bmodule, Bvector
    
    def octupole(self, r, theta, phi, year=2020):
        '''This function calculates the magnetic field octupole
        Parameters
        ----------
        phi : float
            Longitude measured in degrees positive east from Greenwich (in degrees)
        theta : float
            The latitude measured in degrees positive from equator (in degrees)
        r : float
            The distance from the center of the Earth to the point where to calculate the magnetic field (in km)
        Returns
        -------
        B : float
            The magnetic field at the given location
        '''
        gh = self.get_gh_norm(n = 3, m = 3, year = year, coeff = 'b')
        h33 = gh['h']
        g33 = gh['g']
        gh = self.get_gh_norm(n = 3, m = 2, year = year, coeff = 'b')
        h32 = gh['h']
        g32 = gh['g']
        gh = self.get_gh_norm(n = 3, m = 1, year = year, coeff = 'b')
        h31 = gh['h']
        g31 = gh['g']
        g30 = self.get_gh_norm(n = 3, m = 0, year = year, coeff = 'g')['g']
        
        # Calculate the magnetic field

        Bm , B = self.quadrupole(r, theta, phi, year)

        Br_quad = B[0]
        Btheta_quad = B[1]
        Bphi_quad = B[2]

        # Convert to radians
        phi = np.radians(phi)
        theta = np.radians(90-theta)

        # Verify the limits of the input
        if r < 0:
            raise ValueError('The distance must be positive')
        if r < self.Re:
            raise ValueError('The distance must be higher than the radius of the Earth')
        if phi <= -np.pi/2 or phi >= np.pi/2:
            raise ValueError('The co-elevation must be between -90 and 90')
        if theta <= -np.pi or theta >= np.pi:
            raise ValueError('The longitude must be between -180 and 180')

        Br = Br_quad + 4 * (self.Re/r)**5 * (0.5*g30*(np.cos(2*theta) - 1/5) +0.5*(g31 * np.cos(phi) + \
            h31 * np.sin(phi)) * np.sin(theta)*(np.cos(2*theta) - 3/5)+ 0.5*(g32 * np.cos(2*phi) + h32 * np.sin(2*phi)) * \
            np.cos(theta)*(1-np.cos(2*theta)) + 0.5*(g33 * np.cos(3*phi) + h33 * np.sin(3*phi)) * np.sin(theta)*(1-np.cos(2*theta)))
        
        Btheta = Btheta_quad - 3*(self.Re/r)**5 * (0.5*g30*np.sin(theta)*(np.cos(2*theta)-3/5) + 0.5*(g31*np.cos(phi) + \
            h31*np.sin(phi))*np.cos(theta)*(np.cos(2*theta)-2/5) + 0.5*(g32*np.cos(2*phi) + h32*np.sin(2*phi))*np.sin(theta)*\
            (np.cos(2*theta) - 1/3) + 0.5*(g33*np.cos(3*phi) + h33*np.sin(3*phi))*np.cos(theta)*(1-np.cos(2*theta)))
        
        Bphi = Bphi_quad + (self.Re/r)**5 * (0.5*(g31*np.sin(phi) - h31*np.cos(phi))*(np.cos(2*theta)+3/5) + \
            (g32*np.sin(2*phi) - h32*np.cos(2*phi))*np.sin(2*theta) + 3/2 * (g33*np.sin(3*phi) + \
            h33*np.cos(3*phi))*(1-np.cos(2*theta)))
        
        Bvector = np.array([Br, Btheta, Bphi])
        Bmodule = np.linalg.norm(Bvector)

        return Bmodule, Bvector
    
    def magnetic_field(self, r, theta, phi, year=2021, N=3):

        '''This function calculates the magnetic field at a given location
        Parameters
        ----------
        phi : float
            Longitude measured in degrees positive east from Greenwich (in degrees)
        theta : float
            The latitude measured in degrees positive from equator (in degrees)
        r : float
            The distance from the center of the Earth to the point where to calculate the magnetic field (in km)
        N : int
            The order of the magnetic field to calculate
        Returns
        -------
        B : float
            The magnetic field at the given location
        '''
        # Checks to see if located at either pole to avoid singularities
        theta = 90 - theta

        if -0.00000001 < theta < 0.00000001:
            theta = 0.00000001
        elif 179.99999999 < theta < 180.00000001:
            theta = 179.99999999

       
        # Convert to radians
        phi = np.radians(phi)
        theta = np.radians(theta)

        # Verify the limits of the input
        if r < 0:
            raise ValueError('The distance must be positive')
        if r < self.Re:
            raise ValueError('The distance must be higher than the radius of the Earth')
        if phi < -np.pi or phi > np.pi:
            raise ValueError('The co-elevation must be between -180 and 180')
        if theta < 0 or theta > np.pi:
            raise ValueError('The longitude must be between -90 and 90')

        Snm = self.Snm(nm_max=N)
        Pnm, dPnm = self.gauss_norm_ass_leg_poly(nm_max=N, theta=theta)
        # Pnm = self.legendre_poly(N, theta)


        if not 0<N<=13:
            raise ValueError('The order must be between 1 and 13')
         
        Br = 0
        Btheta = 0
        Bphi = 0

        for m in range(N+1):
            for n in range(1, N+1):
                if m <= n:
                    gh = self.get_gh_norm(n, m, year, 'b')
                    g = gh['g']

                    if m != 0:
                        h = gh['h']                        
                    else:
                        h = 0

                    Br += (self.Re / r)**(n+2)* (n + 1) * Pnm[n, m] * (g * np.cos(m*phi) + h*np.sin(m*phi))
                    Btheta += -(self.Re/r)**(n+2) * dPnm[n, m] * (g * np.cos(m*phi) + h * np.sin(m*phi))

                    # Pnm[m, n+1] dPnm[n, m]

                    if theta == 0 or theta == np.pi:
                        Bphi = np.inf
                    else:
                        Bphi -= (self.Re/r)**(n+2) * (-m*g*np.sin(m*phi) + h*m*np.cos(m*phi)) * Pnm[n, m] / np.sin(theta)

                    # print(n, m, Br, Btheta, Bphi)

        Bvector = np.array([Br, Btheta, Bphi])
        Bmodule = np.linalg.norm(Bvector)

        return Bmodule, Bvector

    def gg_to_geo(self, h, gdcolat):
        """
        Compute geocentric colatitude and radius from geodetic colatitude and
        height.

        Alken, P., Thebault, E., Beggan, C. et al, (2020) International Geomagnetic Reference Field: the thirteenth generation, Earth Planets and Space, vol. XX, XX-XX, doi:10.XXX

        Parameters
        ----------
        h : ndarray, shape (...)
            Altitude in kilometers.
        gdcolat : ndarray, shape (...)
            Geodetic colatitude

        Returns
        -------
        radius : ndarray, shape (...)
            Geocentric radius in kilometers.
        theta : ndarray, shape (...)
            Geocentric colatitude in degrees.
        
        sd : ndarray shape (...) 
            rotate B_X to gd_lat 
        cd :  ndarray shape (...) 
            rotate B_Z to gd_lat 

        References
        ----------
        Equations (51)-(53) from "The main field" (chapter 4) by Langel, R. A. in:
        "Geomagnetism", Volume 1, Jacobs, J. A., Academic Press, 1987.
        
        Malin, S.R.C. and Barraclough, D.R., 1981. An algorithm for synthesizing 
        the geomagnetic field. Computers & Geosciences, 7(4), pp.401-405.

        """
        # Use WGS-84 ellipsoid parameters

        eqrad = 6378.137 # equatorial radius
        flat  = 1/298.257223563
        plrad = eqrad*(1-flat) # polar radius
        ctgd  = np.cos(np.deg2rad(gdcolat))
        stgd  = np.sin(np.deg2rad(gdcolat))
        a2    = eqrad*eqrad
        a4    = a2*a2
        b2    = plrad*plrad
        b4    = b2*b2
        c2    = ctgd*ctgd
        s2    = 1-c2
        rho   = np.sqrt(a2*s2 + b2*c2)
        
        rad   = np.sqrt(h*(h+2*rho) + (a4*s2+b4*c2)/rho**2)

        cd    = (h+rho)/rad
        sd    = (a2-b2)*ctgd*stgd/(rho*rad)
        
        cthc  = ctgd*cd - stgd*sd           # Also: sthc = stgd*cd + ctgd*sd
        thc   = np.rad2deg(np.arccos(cthc)) # arccos returns values in [0, pi]

        # Show the results
        print('Geodetic colatitude: %8.3f' % gdcolat)
        print('Geocentric colatitude: %8.3f' % thc)
        print('Geocentric radius: %8.3f' % rad)
        print('Geodetic height: %8.3f' % h)
        
        return rad, thc, sd, cd


    def geo_to_gg(self, radius, theta):
        """
        Compute geodetic colatitude and vertical height above the ellipsoid from
        geocentric radius and colatitude.

        Alken, P., Thebault, E., Beggan, C. et al, (2020) International Geomagnetic Reference Field: the thirteenth generation, Earth Planets and Space, vol. XX, XX-XX, doi:10.XXX

        Parameters
        ----------
        radius : ndarray, shape (...)
            Geocentric radius in kilometers.
        theta : ndarray, shape (...)
            Geocentric colatitude in degrees.

        Returns
        -------
        height : ndarray, shape (...)
            Altitude in kilometers.
        beta : ndarray, shape (...)
            Geodetic colatitude

        Notes
        -----
        Round-off errors might lead to a failure of the algorithm especially but
        not exclusively for points close to the geographic poles. Corresponding
        geodetic coordinates are returned as NaN.

        References
        ----------
        Function uses Heikkinen's algorithm taken from:

        Zhu, J., "Conversion of Earth-centered Earth-fixed coordinates to geodetic
        coordinates", IEEE Transactions on Aerospace and Electronic Systems}, 1994,
        vol. 30, num. 3, pp. 957-961

        """
        
        # Use WGS-84 ellipsoid parameters
        a =  6378.137  # equatorial radius
        b =  6356.752  # polar radius
        
        a2 = a**2
        b2 = b**2

        e2 = (a2 - b2) / a2  # squared eccentricity
        e4 = e2*e2
        ep2 = (a2 - b2) / b2  # squared primed eccentricity

        r = radius * np.sin(np.radians(theta))
        z = radius * np.cos(np.radians(theta))

        r2 = r**2
        z2 = z**2

        F = 54*b2*z2

        G = r2 + (1. - e2)*z2 - e2*(a2 - b2)

        c = e4*F*r2 / G**3

        s = (1. + c + np.sqrt(c**2 + 2*c))**(1./3)

        P = F / (3*(s + 1./s + 1.)**2 * G**2)

        Q = np.sqrt(1. + 2*e4*P)

        r0 = -P*e2*r / (1. + Q) + np.sqrt(
            0.5*a2*(1. + 1./Q) - P*(1. - e2)*z2 / (Q*(1. + Q)) - 0.5*P*r2)

        U = np.sqrt((r - e2*r0)**2 + z2)

        V = np.sqrt((r - e2*r0)**2 + (1. - e2)*z2)

        z0 = b2*z/(a*V)

        height = U*(1. - b2 / (a*V))

        beta = 90. - np.degrees(np.arctan2(z + ep2*z0, r))

        # Show the results
        print('Geocentric colatitude: %8.3f' % theta)
        print('Geodetic colatitude: %8.3f' % beta)
        print('Geocentric radius: %8.3f' % radius)
        print('Geodetic height: %8.3f' % height)

        return height, beta
    
    def transformation2NED(self, Bvector, epsilon=0):
        '''Transformation of the earth's magnetic field from local tangent plane coordinate to NED
        Parameters
        ----------
        Bvector : array
            The magnetic field vector r, theta, phi
        epsilon : float
            epsilon = lambdaa - delta # Elevation (geodetic latitude)
                delta = 90 - theta # Declination
                lambdaa = 90 - phi # Co-elevation
                (degrees)
        Returns
        -------
        B : float
            The magnetic field at the given location in NED
        '''

        if epsilon > 0.2:
            print('Warning: The elevation is higher than 0.2 degrees')

        epsilon = np.radians(epsilon)

        Br = Bvector[0]
        Btheta = Bvector[1]
        Bphi = Bvector[2]

        BN = -Btheta * np.cos(epsilon) - Br * np.sin(epsilon)
        BE = Bphi
        BD = Btheta * np.sin(epsilon) - Br * np.cos(epsilon)

        # geomagnetic field components in NED coordinates
        print('Coordenadas en el sistema de referencia local (NED)')
        print('BN [nT] =', BN)
        print('BE [nT] =', BE)
        print('BD [nT] =', BD)

        return np.array([BN, BE, BD])
    
    def transformation2inertial(self, Bvector, lat, long, thetag = 0):
        '''Transformation of the earth's magnetic field from local tangent plane coordinate to inertial coordinates
        Geocentric inertial components used in satellite work
        Based on Spacecraft Attitude Determination and Control by James R. Wertz
        Parameters
        ----------
        Bvector : array
            The magnetic field vector r, theta, phi
        lat : float
            The latitude (declination) of the point where to calculate the magnetic field (in degrees)
        long : float
            The longitude from Grenwiwh of the point where to calculate the magnetic field (in degrees)
        thetag : float
            The right ascension of the Greenwich meridian or the sidereal time at 
            Greenwich (in degrees)
        Returns
        -------
        B : float
            The magnetic field at the given location in inertial coordinates
        '''
        
        # thetag indicate declination and celestial time in Greenwich

        delta = lat # declination
        alpha = long + thetag # right ascension

        delta = np.radians(delta)
        alpha = np.radians(alpha)


        Br = Bvector[0]
        Btheta = Bvector[1]
        Bphi = Bvector[2]

        BxI = (Br * np.cos(delta) + Btheta * np.sin(delta)) * np.cos(alpha) - Bphi * np.sin(alpha)
        ByI = (Br * np.cos(delta) + Btheta * np.sin(delta)) * np.sin(alpha) + Bphi * np.cos(alpha)
        BzI = Br * np.sin(delta) - Btheta * np.cos(delta)

        # geomagnetic field components in inertial coordinates
        print('Coordenadas en el sistema de referencia inercial (Inertial)')
        print('BxI [nT] =', BxI)
        print('ByI [nT] =', ByI)
        print('BzI [nT] =', BzI)

        return np.array([BxI, ByI, BzI])
    
    def transformation2orb(self, Bvector, orbparam, delta, long, thetag = 0):
        '''Transformation of the earth's magnetic field from inertial coordinates to orbital coordinates
        Ecuaciones comprobadas
        Parameters
        ----------
        Bvector : array
            The magnetic field vector r, theta, phi
        orbparam : array
            The orbital parameters
            - true anomaly (in degrees)
            - Right ascension of ascending node (in degrees)
            - inclination (in degrees)
        delta : float
            The latitude (declination) of the point where to calculate the magnetic field (in degrees)
        long : float
            The longitude from Grenwiwh of the point where to calculate the magnetic field (in degrees)
        thetag : float
            The declination of the Greenwich meridian (in degrees)
        Returns
        -------
        B : float
            The magnetic field at the given location in orbital coordinates
        '''

        alpha = long + thetag 

        # to radians
        delta = np.radians(delta)
        alpha = np.radians(alpha) 

        ta = np.radians(orbparam[0])
        RAAN = np.radians(orbparam[1])
        inclination = np.radians(orbparam[2])

        Br = Bvector[0]
        Btheta = Bvector[1]
        Bphi = Bvector[2]

        # Transformation to orbital frame
        BxO = (-np.sin(ta)*np.cos(RAAN) - np.cos(RAAN)*np.sin(RAAN)*np.cos(inclination))* \
            (Br*np.cos(delta)*np.cos(alpha) + Btheta*np.sin(delta)*np.cos(alpha) - Bphi*np.sin(alpha)) + \
            (-np.sin(ta)*np.sin(RAAN) + np.cos(ta)*np.cos(RAAN)*np.cos(inclination))* \
            (Br*np.cos(delta)*np.sin(alpha) + Btheta*np.sin(delta)*np.sin(alpha) + Bphi*np.cos(alpha)) + \
            (np.cos(ta)*np.sin(inclination))*(Br*np.sin(delta) - Btheta*np.cos(delta))
        
        ByO = -np.sin(inclination)*np.sin(RAAN)*(Br*np.cos(delta)*np.cos(alpha) + Btheta*np.sin(delta)*np.cos(alpha) - Bphi*np.sin(alpha)) + \
                np.sin(inclination)*np.cos(RAAN)*(Br*np.cos(delta)*np.sin(alpha) + Btheta*np.sin(delta)*np.sin(alpha) + Bphi*np.cos(alpha)) - \
                np.cos(inclination)*(Br*np.sin(delta) - Btheta*np.cos(delta))
        
        BzO = (-np.cos(ta)*np.cos(RAAN) + np.sin(ta)*np.sin(RAAN)*np.cos(inclination))* \
            (Br*np.cos(delta)*np.cos(alpha) + Btheta*np.sin(delta)*np.cos(alpha) - Bphi*np.sin(alpha)) + \
            (-np.cos(ta)*np.sin(RAAN) - np.sin(ta)*np.cos(RAAN)*np.cos(inclination))* \
            (Br*np.cos(delta)*np.sin(alpha) + Btheta*np.sin(delta)*np.sin(alpha) + Bphi*np.cos(alpha)) - \
            (np.sin(ta)*np.sin(inclination))*(Br*np.sin(delta) - Btheta*np.cos(delta))

        print('Coordenadas en el sistema de referencia de la órbita (Orbital)') 
        print('BxO [nT] =', BxO)
        print('ByO [nT] =', ByO)
        print('BzO [nT] =', BzO)
        
        return np.array([BxO, ByO, BzO])
    
    def orb2body(self, BvectorOrb, attiparam):
        '''Transformation of the earth’s magnetic field from inertial coordinates to body frame
        Parameters
        ----------
        BvectorOrb : array
            The magnetic field vector in orbital coordinates
        attiparam : array
            The attitude parameters
            - roll (in degrees)
            - pitch (in degrees)
            - yaw (in degrees)
        Returns
        -------
        B : float
            The magnetic field at the given location in body frame
        '''
        # to radians
        roll = np.radians(attiparam[0])
        pitch = np.radians(attiparam[1])
        yaw = np.radians(attiparam[2])

        BxO = BvectorOrb[0]
        ByO = BvectorOrb[1]
        BzO = BvectorOrb[2]
        
        # transforming the magnetic field to the satellite body frame

        BxB = BxO*np.cos(yaw)*np.cos(pitch) + ByO*(np.cos(yaw)*np.sin(pitch)*np.sin(roll) - np.sin(yaw)*np.cos(roll)) + \
            BzO*(np.cos(yaw)*np.sin(pitch)*np.cos(roll) + np.sin(yaw)*np.sin(roll))
        ByB = BxO*np.sin(yaw)*np.cos(pitch) + ByO*(np.sin(yaw)*np.sin(pitch)*np.sin(roll) + np.cos(yaw)*np.cos(roll)) + \
            BzO*(np.sin(yaw)*np.sin(pitch)*np.cos(roll) - np.cos(yaw)*np.sin(roll))
        BzB = -BxO*np.sin(pitch) + ByO*np.cos(pitch)*np.sin(roll) + BzO*np.cos(pitch)*np.cos(roll)

        print('Coordenadas en el sistema de referencia del satélite (Cuerpo)')
        print('BxORB [nT] =', BxB)
        print('ByORB [nT] =', ByB)
        print('BzORB [nT] =', BzB)

        return np.array([BxB, ByB, BzB])
    
    def linearOrb2body(self, BvectorOrb, attiparam):
        '''Transformation of the earth’s magnetic field from inertial coordinates to body frame
        Parameters
        ----------
        Bvector : array

        '''
        # to radians
        roll = np.radians(attiparam[0])
        pitch = np.radians(attiparam[1])
        yaw = np.radians(attiparam[2])

        BxO = BvectorOrb[0]
        ByO = BvectorOrb[1]
        BzO = BvectorOrb[2]
        
        # transforming the magnetic field to the satellite body frame

        BxB = BxO - yaw*ByO - pitch*BzO
        ByB = yaw*BxO + ByO - roll*BzO
        BzB = -pitch*BxO + roll*ByO + BzO

        print('Coordenadas en el sistema de referencia del satélite (Cuerpo)')
        print('Versión lineariada')
        print('BxORB [nT] =', BxB)
        print('ByORB [nT] =', ByB)
        print('BzORB [nT] =', BzB)

        return np.array([BxB, ByB, BzB])
    
    def xyz2dhif(self, x, y, z):
        '''Calculate D, H, I and F from (X, Y, Z)
        
        Based on code from D. Kerridge, 2019
        
        Parameters
        ---------------
        X: north component (nT) : float
        Y: east component (nT) : float
        Z: vertical component (nT) : float
        
        Returns
        ------
        A tuple: (D, H, I, F) : float
        D: declination (degrees) : float
        H: horizontal intensity (nT) : float
        I: inclination (degrees) : float
        F: total intensity (nT) : float
        '''
        hsq = x*x + y*y
        hoz  = np.sqrt(hsq)
        eff = np.sqrt(hsq + z*z)
        dec = np.arctan2(y,x)
        inc = np.arctan2(z,hoz)

        print('D = ', np.degrees(dec))
        print('H = ', hoz)
        print('I = ', np.degrees(inc))
        print('F = ', eff)
        
        return np.degrees(dec), hoz, np.degrees(inc), eff
    
    def arrayofpoints(self, r, theta, phi, year=2020, N=3, savedata=False):
        '''This function calculates the magnetic field at multiple locations
        Parameters
        ----------
        phi : array
            The longitude measured in degrees positive east from Greenwich (in degrees)
        theta : array
            The latitude measured in degrees positive from equator (in degrees)
        r : array
            The distance from the center of the Earth to the point where to calculate the magnetic field (in km)
        N : int
            The order of the magnetic field to calculate
        savedata : bool
            Save the data to a file
        Returns
        -------
        B : float
            The magnetic field at the given location
        '''

        # Checks the length of the arrays
        if len(phi) != len(theta) or len(phi) != len(r) or len(theta) != len(r):
            raise ValueError('The arrays must have the same length')
        
        # Verify the limits of the input
        for r_i in r:
            if r_i < 0:
                raise ValueError('The distance must be positive')
            if r_i < self.Re:
                raise ValueError('The distance must be higher than the radius of the Earth')
            
        for phi_i in phi:
            if phi_i < -180 or phi_i > 180:
                raise ValueError('The longitude must be between -180 and 180')
        
        for theta_i in theta:
            if theta_i < -90 or theta_i > 90: # latitude
                raise ValueError('The latitude must be between -90 and 90')
        
        B_array = np.zeros((len(phi), 9))
        for i in range(len(phi)):
            B = self.magnetic_field(r[i], theta[i], phi[i], year=year, N=N)[1]

            xyz = self.transformation2NED(B)

            # xyz = self.transformation2orb(xyz, [0, 0, 60], phi[i], theta[i], thetag=0)

            B_array[i, 0] = round(r[i], 2) # distance
            B_array[i, 1] = round(theta[i], 2) # latitude
            B_array[i, 2] = round(phi[i], 2) # longitude
            B_array[i, 3] = round(B[0], 2) # Br
            B_array[i, 4] = round(B[1], 2) # Btheta
            B_array[i, 5] = round(B[2], 2) # Bphi
            B_array[i, 6] = round(xyz[0], 2) # BN
            B_array[i, 7] = round(xyz[1], 2) # BE
            B_array[i, 8] = round(xyz[2], 2) # BD


        if savedata:
            path = os.getcwd() + '/results/B_array_latlon.csv'
            data = pd.DataFrame(B_array, columns=['r (km)', 'theta (deg)', 'phi (deg)', 'Br (nT)', 'Btheta (nT)', 'Bphi (nT)', 'BN (nT)', 'BE (nT)', 'BD (nT)'])
            data.to_csv(path)

            print('Datos guardados en ' + path)

        return B_array
    
    def arrayDatesAtLocation(self, val, years, N=3, savedata=False):
        '''This function calculates the magnetic field at one location for multiple dates
        Parameters
        ----------
        val : array
            The location where to calculate the magnetic field
            [r, theta, phi]
        years : array
            The year for which to calculate the coefficients
        N : int
            The order of the magnetic field to calculate
        savedata : bool
            Save the data to a file
        Returns
        -------
        B : float
            The magnetic field at the given location
        '''

        # Checks the bounds of the input
        if len(val) != 3:
            raise ValueError('The input must be an array of length 3')

        # Verify the limits of the input
        if val[0] < 0:
            raise ValueError('The distance must be positive')
        if val[0] < self.Re:
            raise ValueError('The distance must be higher than the radius of the Earth')
        if val[2] < -180 or val[2] > 180:
            raise ValueError('The longitude must be between -180 and 180')
        if val[1] < -90 or val[1] > 90:
            raise ValueError('The latitude must be between -90 and 90')
        
        # create a dataframes to store the results
        B_array = np.zeros((len(years), 11))

        for i, year in enumerate(years):
            B = self.magnetic_field(val[0], val[1], val[2], year=year, N=N)[1]

            xyz = self.transformation2NED(B)

            dhif = self.xyz2dhif(xyz[0], xyz[1], xyz[2])

            xyz = self.transformation2orb(xyz, [0, 0, 60], val[2], val[1], thetag=0)

            B_array[i, 0] = year

            B_array[i, 1] = round(B[0], 2)
            B_array[i, 2] = round(B[1], 2)
            B_array[i, 3] = round(B[2], 2)
            B_array[i, 4] = round(xyz[0], 2)
            B_array[i, 5] = round(xyz[1], 2)
            B_array[i, 6] = round(xyz[2], 2)
            B_array[i, 7] = round(dhif[0], 2)
            B_array[i, 8] = round(dhif[1], 2)
            B_array[i, 9] = round(dhif[2], 2)
            B_array[i, 10] = round(dhif[3], 2)


        if savedata:
            path = os.getcwd() + '/results/B_array_dates.csv'
            data = pd.DataFrame(B_array, columns=['year', 'Br (nT)', 'Btheta (nT)', 'Bphi (nT)', 'BN (nT)', 'BE (nT)', 'BD (nT)', 'D (deg)', 'H (nT)', 'I (deg)', 'F (nT)'])
            data.to_csv(path)

        return B_array

    
    def printresults(self, rtp, savedata=False):
        '''This function prints the magnetic field at a given location
        Parameters
        ----------
        val : array
            The location where to calculate the magnetic field
            [theta, phi, r]
        year : float or int
            The year for which to calculate the coefficients
        N : int
            The order of the magnetic field to calculate
        savedata : bool
            Save the data to a file      
        Returns
        -------
        Geomagnetic field values at: 10.0� / 10.0�, at altitude 7000.0, for 2020.0
        Declination (D):-0.63�
        Inclination (I):-1.06�
        Horizontal intensity (H): 24435.2nT
        Total intensity (F)     : 24439.4nT
        North component (X)     : 24433.7nT
        East component (Y)      :  -268.7nT
        Vertical component (Z)  :  -454.0nT
        '''

        xyz = self.transformation2NED(rtp)

        dhif = self.xyz2dhif(xyz[0], xyz[1], xyz[2])

        print('Geomagnetic field values at: ')
        print('Declination (D): ' + str(dhif[0]) + ' degrees')
        print('Inclination (I): ' + str(dhif[2]) + ' degrees')
        print('Horizontal intensity (H): ' + str(dhif[1]) + ' nT')
        print('Total intensity (F): ' + str(dhif[3]) + ' nT')
        print('North component (X): ' + str(xyz[0]) + ' nT')
        print('East component (Y): ' + str(xyz[1]) + ' nT')
        print('Vertical component (Z): ' + str(xyz[2]) + ' nT')

        if savedata:
            path = os.getcwd() + '/results/magPoint.csv'

            # Write to file
            with open(path, 'w') as f:
                f.write('Geomagnetic field values \n')
                f.write('Declination (D): ' + str(dhif[0]) + ' degrees\n')
                f.write('Inclination (I): ' + str(dhif[2]) + ' degrees\n')
                f.write('Horizontal intensity (H): ' + str(dhif[1]) + ' nT\n')
                f.write('Total intensity (F): ' + str(dhif[3]) + ' nT\n')
                f.write('North component (X): ' + str(xyz[0]) + ' nT\n')
                f.write('East component (Y): ' + str(xyz[1]) + ' nT\n')
                f.write('Vertical component (Z): ' + str(xyz[2]) + ' nT\n')

                print('Results saved to /results/magPoint.csv')
    
    def plotMagneticField(self, val, year=2020, N=13, modelos=[1, 0, 0, 0, 1], timeRange=[1900, 2025], absolute_value=False):
        '''This function plots the magnetic field at a given location
        Parameters
        ----------
        val : array
            The location where to calculate the magnetic field
            [r, theta phi]
        year : float or int
            The year for which to calculate the coefficients
        N : int
            The order of the magnetic field to calculate
        modelos : array
            The models to plot 
            If 1 the model will be plotted, if 0 the model will not be plotted
            The array position corresponds to the following models
            [dup,cdip, quad, octo, general]
            [0/1, 0/1,  0/1,  0/1, 0/1]
        timeRange : array
            The time range to plot the magnetic field
            [start, end]
            Limits: 1900 - 2025
        absolute_value : bool
            If True the absolute value of the magnetic field will be plotted
        Returns
        -------
        Br : array
            The radial component of the magnetic field
        Btheta : array
            The theta component of the magnetic field
        Bphi : array
            The phi component of the magnetic field

        Examples
        --------
        geomag.plotMagneticField(val=[7000, 10, 10], 
                            year=2020, 
                            N=13, 
                            modelos= [1, 0, 0, 0, 0], 
                            timeRange=[2020, 2025],
                            absolute_value=True)
        '''

        r = val[0]
        theta = val[1]
        phi = val[2]

        start = timeRange[0]
        end = timeRange[1]

        ntotal = 1000
        time = np.linspace(start, end, ntotal)

        Br = np.zeros([len(modelos), len(time)])
        Btheta = np.zeros([len(modelos), len(time)])
        Bphi = np.zeros([len(modelos), len(time)])

        list_modelos = ['Dipolo', 'Dipolo centrado', 'Cuadrupolo', 'Octupolo', 'Modelo completo de Orden '+ str(N)]
        list_colors = ['b', 'r', 'g', 'y', 'k']
        legend = []

        if absolute_value:
            fig, axs = plt.subplots(1, 1, figsize=(10, 10))
            for idx, model in enumerate(modelos):
                if model==1:
                    for idxt, t in enumerate(time):
                        if idx == 0:
                            val = self.dipole(r, theta, phi, year=t)[0]
                        elif idx == 1:
                            val = self.centered_dipole(r, theta, phi, year=t)[0]
                        elif idx == 2:
                            val = self.quadrupole(r, theta, phi, year=t)[0]
                        elif idx == 3:
                            # val = geomag.octupole(phi, theta, r, year=t)[0]
                            val = self.magnetic_field(r, theta, phi, year=t, N=3)[0]
                        elif idx == 4:
                            val = self.magnetic_field(r, theta, phi, year=t, N=N)[0]

                        Br[idx, idxt] = val


                    axs.plot(time, Br[idx, :], list_colors[idx], label=list_modelos[idx])

                    legend.append(list_modelos[idx])

            axs.set_ylabel('B (nT)')
            axs.legend(legend, loc='upper right')
            axs.set_title('Campo magnético')
            axs.set_xlabel('Año')

            axs.grid( which='major', axis='both')

        else:
            fig, axs = plt.subplots(3, 1, figsize=(10, 10))
            for idx, model in enumerate(modelos):
                if model==1:
                    for idxt, t in enumerate(time):
                        if idx == 0:
                            val = self.dipole(r, theta, phi, year=t)[1]
                        elif idx == 1:
                            val = self.centered_dipole(r, theta, phi, year=t)[1]
                        elif idx == 2:
                            val = self.quadrupole(r, theta, phi, year=t)[1]
                        elif idx == 3:
                            # val = geomag.octupole(phi, theta, r, year=t)[1]
                            val = self.magnetic_field(r, theta, phi, year=t, N=3)[1]
                        elif idx == 4:
                            val = self.magnetic_field(r, theta, phi, year=t, N=N)[1]
                            # val = geomag.magnet(r, theta, phi, year=t, N=N)

                        Br[idx, idxt] = val[0]
                        Btheta[idx, idxt] = val[1]
                        Bphi[idx, idxt] = val[2]


                    axs[0].plot(time, Br[idx, :], list_colors[idx], label=list_modelos[idx])
                    axs[0].set_ylabel('B_r (nT)')
                    axs[1].plot(time, Btheta[idx, :], list_colors[idx], label=list_modelos[idx])
                    axs[1].set_ylabel('B_theta (nT)')
                    axs[2].plot(time, Bphi[idx, :], list_colors[idx], label=list_modelos[idx])
                    axs[2].set_ylabel('B_phi (nT)')
                    

                    legend.append(list_modelos[idx])

            axs[0].legend(legend, loc='upper right')
            axs[0].set_title('Campo magnético')
            axs[2].set_xlabel('Año')

            axs[0].grid( which='major', axis='both')
            axs[1].grid( which='major', axis='both')
            axs[2].grid( which='major', axis='both')

        plt.show()

        return Br, Btheta, Bphi
    
    def plotMagneticFieldGeneral(self, time, Bvector):
        '''This function plots the magnetic field at a given location
        Parameters
        ----------
        time : array
            The time for which to calculate the magnetic field
        Bvector : array
            The magnetic field vector r, theta, phi
        Returns
        -------
        None
        '''

        Bx = [x[0] for x in Bvector]
        By = [x[1] for x in Bvector]
        Bz = [x[2] for x in Bvector]

        fig, axs = plt.subplots(3, 1, figsize=(10, 10))

        axs[0].plot(time, Bx, 'b', label='Bx')
        axs[0].set_ylabel('B_x (nT)')
        axs[1].plot(time, By, 'r', label='By')
        axs[1].set_ylabel('B_y (nT)')
        axs[2].plot(time, Bz, 'g', label='Bz')
        axs[2].set_ylabel('B_z (nT)')

        axs[0].set_title('Campo magnético')
        axs[2].set_xlabel('Tiempo')

        axs[0].grid( which='major', axis='both')
        axs[1].grid( which='major', axis='both')
        axs[2].grid( which='major', axis='both')

        plt.show()

        return fig, axs
                                                 

if __name__ == '__main__':
    geomag = Geomag()

    print(geomag.magnetic_field(7000, 10, 10,year=2020, N=3))
    print('dipole')
    print(geomag.dipole(7000, 10, 10, year=2020))
    print(geomag.centered_dipole(7000, 10, 10, year=2020))
    print('quadrupole')
    print(geomag.quadrupole(7000, 10, 10, year=2020))
    print('octupole')
    print(geomag.octupole(7000, 10, 10, year=2020))

    print(geomag.transformation2NED(geomag.magnetic_field(7000, 10, 10, N=12, year=2020)[1]))
    print(geomag.transformation2orb(geomag.magnetic_field(7000, 10, 10, N=12, year=2020)[1], [0, 0, 60], 10, 10))

    # geomag.printresults(7000, 10, 10, year=2020, N=13)

    geomag.plotMagneticField(val=[7000, 10, 10], 
                    year=2020, 
                    N=13, 
                    modelos= [0, 0, 0, 0, 1], 
                    timeRange=[1900, 2025],
                    absolute_value=False)

    B = geomag.arrayDatesAtLocation([7000, 10, 10], years=[2020, 2021, 2022, 2023, 2024, 2025], N=13, savedata=True)
    plt.figure()
    plt.plot(B[:, 3], 'b', label='Br')
    plt.plot(B[:, 4], 'r', label='Btheta')
    plt.plot(B[:, 5], 'g', label='Bphi')
    plt.show()
    '''

    r = 7000 * np.ones(11)
    theta = np.arange(-90, 91, 36/2)
    phi = np.arange(-180, 181, 18*2)
    
    geomag.arrayofpoints(r, theta, phi, year=2020, N=13, savedata=True)
'''
    


