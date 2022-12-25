import numpy as np
from scipy.special import lpmn, factorial

class Geopot():
    '''Class to calculate the potential arround the Earth'''

    def __init__(self, body, r, elevation, azimuth):
        self.mu = body.mu * 10**3  # Gravitational parameter
        self.a = 6378136.3         # Earth equatorial radius
        self.r = r                 # Radius
        self.elev = elevation      # Elevation angle
        self.azi = azimuth         # Azimuth angle
        
        self.grav = np.array([-self.mu / self.r**2, 0, 0])
        self.pot = -self.mu / self.r

    def gravitational_potential(self):
        '''Gravitational potential 4.42'''
        
        # Read the data
        data = np.loadtxt('data\egm96_to360.ascii.txt')
        n, m = np.shape(data)

        # Max column 1 and 2
        nn = np.max(data[:, 0])
        mm = np.max(data[:, 1])
        N_max = int(np.max([nn, mm]))

        # Legendre polynomials
        P, Pd = self.legendre(N_max, np.sin(self.elev))

        for i in range(N_max):

            n = int(data[i, 0])
            m = int(data[i, 1])
            C = data[i, 2]
            S = data[i, 3]

            Pn = P[m, n]
            Pnd = Pd[m, n]
            
            self.pot = -self.pot + self.mu /self.r * (self.a / self.r)**n * Pn * (C * np.cos(m * self.azi) + S * np.sin(m * self.azi))

            self.grav[0] = self.grav[0] - self.mu * (n + 1) / self.r**(n +2) * self.a**n * Pn * (C * np.cos(m * self.azi) + S * np.sin(m * self.azi))
            self.grav[1] = self.mu / self.r**(n+2) * self.a**n * Pnd * (C * np.cos(m * self.azi) + S * np.sin(m * self.azi)) *  np.cos(self.elev) # derivada
            self.grav[2] = self.mu / self.r**(n+2) / np.sin(self.elev) * self.a**n * Pn * m *(S * np.cos(m * self.azi) - C * np.sin(m * self.azi))
            
            print(self.pot)
            print(self.grav)

    def legendre(self, N_max, x):
        '''Legendre polynomials
        n: degree
        m: order
        x: argument
        Returns: 
        Legendre polynomial normalized evaluated at x
        '''
        P, Pd= lpmn(N_max, N_max, x)

        # Normalization

        for m in range(N_max):
            for n in range(N_max):
                if m==0:
                    P[m, n] = P[m, n] * np.sqrt(((2*n+1)*factorial(n-m))/factorial(n+m))
                else:
                    P[m, n] = P[m, n] * np.sqrt(((2*n+1)*factorial(n-m))/factorial(n+m)) * np.sqrt(2)
        
        return P, Pd




