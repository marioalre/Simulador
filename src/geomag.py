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
            if m != 0:
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

        Br = 2 * (self.Re/r)**3 * (g10 * np.cos(theta) + (g11 * np.cos(phi) + h11 * np.sin(phi)) * np.sin(theta))
        Btheta = (self.Re/r)**3 * (g10 * np.sin(theta) - (g11 * np.cos(phi) + h11 * np.sin(phi)) * np.cos(theta))
        Bphi = (self.Re/r)**3 * (-h11 * np.cos(phi) + g11 * np.sin(phi))

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

        Br = 2 * (self.Re/r)**3 * g10 * np.cos(theta)
        Btheta = (self.Re/r)**3 * g10 * np.sin(theta)
        Bphi = 0

        Bvector = np.array([Br, Btheta, Bphi])
        Bmodule = np.linalg.norm(Bvector)
        
        return Bmodule, Bvector
    
    def quadrupole(self, phi, theta, r, year=2020):
        '''This function calculates the magnetic field quadrupole
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
    
    def octupole(self, phi, theta, r, year=2020):
        '''This function calculates the magnetic field octupole
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

        Bm , B = self.quadrupole(phi, theta, r, year)

        Br_quad = B[0]
        Btheta_quad = B[1]
        Bphi_quad = B[2]

        Br = Br_quad + 4 * (self.Re/r)**5 * (0.5*g30*(np.cos(2*theta) + 1/5) +0.5*(g31 * np.cos(phi) + \
            h31 * np.sin(phi)) * np.sin(theta)*(np.cos(2*theta) - 3/5)+ 0.5*(g32 * np.cos(2*phi) + h32 * np.sin(2*phi)) * \
            np.cos(theta)*(1-np.cos(2*theta)) + 0.5*(g33 * np.cos(3*phi) + h33 * np.sin(3*phi)) * np.sin(theta)*(1-np.cos(2*theta)))
        
        Btheta = Btheta_quad - 3*(self.Re/r)**5 * (0.5*g30*np.sin(theta)*(np.cos(2*theta)-3/5) + 0.5*(g31*np.cos(phi) + \
            h31*np.sin(phi))*np.cos(theta)*(np.cos(2*theta)-2/5) - 0.5*(g32*np.cos(2*phi) + h32*np.sin(2*phi))*np.sin(theta)*\
            (np.cos(2*theta) - 1/3) + 0.5*(g33*np.cos(3*phi) + h33*np.sin(3*phi))*np.cos(theta)*(1-np.cos(2*theta)))
        
        Bphi = Bphi_quad + (self.Re/r)**5 * (0.5*(g31*np.sin(phi) - h31*np.cos(phi))*(np.cos(2*theta)+3/5) + \
            (g32*np.sin(2*phi) - h32*np.cos(2*phi))*np.sin(2*theta) + 3/2 * (g33*np.sin(3*phi) + \
            h33*np.cos(3*phi))*(1-np.cos(2*theta)))
        
        Bvector = np.array([Br, Btheta, Bphi])
        Bmodule = np.linalg.norm(Bvector)

        return Bmodule, Bvector
    
    def magnetic_field(self, phi, theta, r, year=2020, N=3):

        '''This function calculates the magnetic field at a given location
        Parameters
        ----------
        phi : float
            The co-elevation of the point where to calculate the magnetic field (in degrees)
        theta : float
            The  East longitude of the point where to calculate the magnetic field (in degrees)
        r : float
            The distance from the center of the Earth to the point where to calculate the magnetic field (in km)
        N : int
            The order of the magnetic field to calculate
        Returns
        -------
        B : float
            The magnetic field at the given location
        '''
       
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


        Pnm, dPnm = self.gauss_norm_ass_leg_poly(nm_max=N, theta=theta)


        if not 0<N<self.max_mn:
            raise ValueError('The order must be between 1 and 13')
         

        Br = 0
        Btheta = 0
        Bphi = 0

        for n in range(1, N+1):
            for m in range(0, N+1):
                if n >= m:
                    gh = self.get_gh_norm(n, m, year, 'b')
                    g = gh['g']

                    if m != 0:
                        h = gh['h']                        
                    else:
                        h = 0

                    Br += (self.Re / r)**(n+2)* (n + 1) * Pnm[n, m] * (g * np.cos(m*phi) + h*np.sin(m*phi))
                    Btheta -= (self.Re/r)**(n+2) * dPnm[n, m] * (g * np.cos(m*phi) + h * np.sin(m*phi))

                    if theta == 0 or theta == np.pi:
                        Bphi = np.inf
                    else:
                        Bphi -= (self.Re/r)**(n+2) * (-m*g*np.sin(m*phi) + h*m*np.cos(m*phi)) * Pnm[n, m] / np.sin(theta)

        Bvector = np.array([Br, Btheta, Bphi])
        Bmodule = np.linalg.norm(Bvector)

        return Bmodule, Bvector
    
    def transformation2NED(self, r, theta, phi, lambdaa, Bvector):
        '''Transformation of the earth’s magnetic field from local tangent plane coordinate to NED
        Parameters
        ----------
        r : float
            The distance from the center of the Earth to the point where to calculate the magnetic field (in km)
        theta : float
            The  East longitude of the point where to calculate the magnetic field (in degrees)
        phi : float
            The co-elevation of the point where to calculate the magnetic field (in degrees)
        lambda : float
            Geodetic longitude of the point where to calculate the magnetic field (in degrees)
        Returns
        -------
        B : float
            The magnetic field at the given location in NED
        '''

        # Convert to radians


        delta = 90 - theta
        epsilon = lambdaa - delta

        delta = np.radians(delta)
        epsilon = np.radians(epsilon)

        Br = Bvector[0]
        Btheta = Bvector[1]
        Bphi = Bvector[2]

        BN = -Btheta * np.cos(delta) - Br * np.sin(delta)
        BE = Bphi
        BD = Btheta * np.sin(delta) - Br * np.cos(delta)

        return np.array([BN, BE, BD])
    
    def transformation2inertial(self, Bvector, r, theta, phi, thetag = 0):
        '''Transformation of the earth’s magnetic field from local tangent plane coordinate to inertial coordinates
        Parameters
        ----------
        r : float
            The distance from the center of the Earth to the point where to calculate the magnetic field (in km)
        theta : float
            The  East longitude of the point where to calculate the magnetic field (in degrees)
        phi : float
            The co-elevation of the point where to calculate the magnetic field (in degrees)
        thetag : float
            The declination of the Greenwich meridian (in degrees)
        Bvector : array
            The magnetic field vector r, phi, theta
        Returns
        -------
        B : float
            The magnetic field at the given location in inertial coordinates
        '''
        
        # thetag indicate declination and celestial time in Greenwich

        delta = 90 - phi
        alpha = theta + thetag

        delta = np.radians(delta)
        alpha = np.radians(alpha)


        Br = Bvector[0]
        Btheta = Bvector[1]
        Bphi = Bvector[2]

        BxI = (Br * np.cos(delta) + Btheta * np.sin(delta)) * np.cos(alpha) - Bphi * np.sin(alpha)
        ByI = (Br * np.cos(delta) + Btheta * np.sin(delta)) * np.sin(alpha) + Bphi * np.cos(alpha)
        BzI = Br * np.sin(delta) - Btheta * np.cos(delta)

        # geomagnetic field components in inertial coordinates

        return np.array([BxI, ByI, BzI])
    
    def transformation2orb(self, Bvector, orbparam, phi, theta, thetag = 0):
        '''Transformation of the earth’s magnetic field from inertial coordinates to body frame
        Parameters
        ----------
        Bvector : array
            The magnetic field vector r, phi, theta
        orbparam : array
            The orbital parameters
            - true anomaly (in degrees)
            - Right ascension of ascending node (in degrees)
            - inclination (in degrees)
        attiparam : array
            The attitude parameters
            - roll (in degrees)
            - pitch (in degrees)
            - yaw (in degrees)
        phi : float
            The co-elevation of the point where to calculate the magnetic field (in degrees)
        theta : float
            The  East longitude of the point where to calculate the magnetic field (in degrees)
        thetag : float
            The declination of the Greenwich meridian (in degrees)
        Returns
        -------
        B : float
            The magnetic field at the given location in body frame
        '''

        delta = 90 - phi
        alpha = theta + thetag

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
                np.sin(inclination)*np.cos(RAAN)*(Br*np.cos(delta)*np.sin(alpha) + Btheta*np.sin(delta)*np.sin(alpha) + Bphi*np.cos(alpha)) + \
                np.cos(inclination)*(Br*np.sin(delta) - Btheta*np.cos(delta))
        
        BzO = (-np.cos(ta)*np.cos(RAAN) - np.sin(RAAN)*np.cos(RAAN)*np.sin(inclination))* \
            (Br*np.cos(delta)*np.cos(alpha) + Btheta*np.sin(delta)*np.cos(alpha) - Bphi*np.sin(alpha)) + \
            (-np.cos(ta)*np.sin(RAAN) - np.sin(ta)*np.cos(RAAN)*np.cos(inclination))* \
            (Br*np.cos(delta)*np.sin(alpha) + Btheta*np.sin(delta)*np.sin(alpha) + Bphi*np.cos(alpha)) + \
            (np.sin(ta)*np.sin(inclination))*(Br*np.sin(delta) - Btheta*np.cos(delta))
        
        return np.array([BxO, ByO, BzO])
    
    def orb2body(self, BvectorOrb, attiparam):
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

        BxB = BxO*np.cos(yaw)*np.cos(pitch) + ByO*(np.cos(yaw)*np.sin(pitch)*np.sin(roll) - np.sin(yaw)*np.cos(roll)) + \
            BzO*(np.cos(yaw)*np.sin(pitch)*np.cos(roll) + np.sin(yaw)*np.sin(roll))
        ByB = BxO*np.sin(yaw)*np.cos(pitch) + ByO*(np.sin(yaw)*np.sin(pitch)*np.sin(roll) + np.cos(yaw)*np.cos(roll)) + \
            BzO*(np.sin(yaw)*np.sin(pitch)*np.cos(roll) - np.cos(yaw)*np.sin(roll))
        BzB = -BxO*np.sin(pitch) + ByO*np.cos(pitch)*np.sin(roll) + BzO*np.cos(pitch)*np.cos(roll)

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

        return np.array([BxB, ByB, BzB])
                                              

if __name__ == '__main__':
    geomag = Geomat()
    data = geomag.get_txt()
    print(data)
    
    print(geomag.Snm())
    print(geomag.get_gh_data(2, 0, '2000.0', 'b'))
    print(geomag.dipole(10, 10, 7000))
    print(geomag.centered_dipole(10, 10, 7000))
    print(geomag.quadrupole(10, 10, 7000))
    print(geomag.octupole(10, 10, 7000))


    print(geomag.magnetic_field(10, 10, 7000, N=1))
    # a, b = geomag.gauss_norm_ass_leg_poly(2, np.pi/4)




