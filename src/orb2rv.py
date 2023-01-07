# Class to transform orbital parameters to position and velocity vectors

import numpy as np
from src.Orbit import Orbit

class Orb2rv(Orbit):
    '''Orbital parameters to position and velocity vectors'''

    def __init__(self,a=None, h=None, e=None, Omega=None, omega=None, i=None, nu=None, M=None, p=None, n=None, body=None):
        '''Initialize the class
        Parameters
        ----------
        a : float
            Semi-major axis in km
        h : float
            Specific angular momentum in km^2/s
        e : float
            Eccentricity
        Omega : float
            Right ascension of the ascending node in degrees
        omega : float
            Argument of periapsis in degrees
        i : float
            Inclination in degrees
        nu : float
            True anomaly in degrees
        M : float
            Mean anomaly
        p : float
            Semi-latus rect
        n : float
            Mean motion in rad/s
        body : Body
            Body object
        '''

        # Verify that at least 6 parameters are defined
        ver =(a is None) + (h is None) + (e is None) + (Omega is None)
        ver += (omega is None) + (i is None) + (nu is None) + (M is None) + (p is None) + (n is None)
        if  ver < 5:
            Warning('At least 6 parameters must be defined')

        self.a = a
        self.h = h
        self.n = n
        self.e = e

        if a is not None:
            self.a = a
        elif self.h is not None:
            self.a = self.h**2 / body.mu
        elif self.n is not None:
            self.a = (body.mu / self.n**2)**(1/3)
        else:
            raise ValueError('Either a, h or e must be defined')

        self.p = p

        if h is not None:
            self.h = h
        elif self.p is not None:
            self.h = self.hfromp(self.p)
        elif (self.e is not None) and (self.a is not None):
            self.h = (self.a * body.mu * (1 - self.e**2))**0.5
        else:
            raise ValueError('Either h, p or a and e must be defined')
        
        self.Omega = Omega * np.pi / 180
        self.omega = omega * np.pi / 180
        self.i = i * np.pi / 180
        self.M = M

        if nu is not None:
            self.nu = nu * np.pi / 180
        elif self.M is not None:
            self.M = self.M * np.pi / 180
            self.nu = self.nufromM()
        else:
            raise ValueError('Either nu or M must be defined')

        self.mu = body.mu
        self.radius = body.radius

        self.r = self.position_eci()
        self.v = self.velocity_eci()


    def position(self):
        '''Position vector 4.37'''
        return (self.h**2 / self.mu) * (1 / (1 + self.e * np.cos(self.nu))) * np.array([np.cos(self.nu), np.sin(self.nu), 0])

    def velocity(self):
        '''Velocity vector 4.38'''
        return (self.mu/self.h) * np.array([-np.sin(self.nu), self.e + np.cos(self.nu), 0])

    def rotation_matrix(self):
        '''Rotation matrix 4.39 4.40 4.41 4.44'''
        R1 = np.array([[np.cos(self.Omega), np.sin(self.Omega), 0],
                       [-np.sin(self.Omega), np.cos(self.Omega), 0],
                       [0, 0, 1]])

        R2 = np.array([[1, 0, 0],
                       [0, np.cos(self.i), np.sin(self.i)],
                       [0, -np.sin(self.i), np.cos(self.i)]])

        R3 = np.array([[np.cos(self.omega), np.sin(self.omega), 0],
                       [-np.sin(self.omega), np.cos(self.omega), 0],
                       [0, 0, 1]])

        return np.transpose(R1) @ np.transpose(R2) @ np.transpose(R3)

    def position_eci(self):
        '''Position vector ECI 4.46'''
        return self.rotation_matrix() @ self.position()

    def velocity_eci(self):
        '''Velocity vector ECI 4.46'''
        return self.rotation_matrix() @ self.velocity()

    def hfromp(self):
        '''Get h from p'''
        return np.sqrt(self.mu * self.p)

    def nufromM(self):
        '''Get nu from M'''
        e = self.e
        M = self.M

        E = self.EfromM()
        return 2 * np.arctan(np.sqrt((1 + e) / (1 - e)) * np.tan(E / 2))

    def EfromM(self):
        '''Get E from M'''
        e = self.e
        M = self.M

        E = M + e * np.sin(M)
        E_old = E + 1
        while np.abs(E - E_old) > 1e-8:
            E_old = E
            E = E - (E - e * np.sin(E) - M) / (1 - e * np.cos(E))
        return E


    





