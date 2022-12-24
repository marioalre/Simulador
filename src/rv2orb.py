# Class to convert position and velocity vectors to orbital parameters

import numpy as np
from src.Orbit import Orbit

class Rv2orb(Orbit):
    '''Position and velocity vectors to orbital parameters'''
    def __init__(self, r, v, body, t0, t):
        '''Initialize the class
        Parameters
        ----------
        r : numpy array
            Position vector ECI
        v : numpy array
            Velocity vector ECI
        body : Body
            Body object
        t0 : float
            Initial time
        t : float
            Final time
        '''
        self.r = r   # position vector ECI 
        self.v = v   # velocity vector ECI
        self.mu = body.mu
        self.h = np.cross(self.r, self.v)
        self.dt = t - t0
        self.e = self.eccentricity()
        self.a = self.semi_major_axis()
        self.i = self.inclination()
        self.Omega = self.ascending_node()
        self.omega = self.argument_of_periapsis()
        self.nu = self.true_anomaly()
        self.n = self.mean_motion()
        self.T = self.period()
    
    def eccentricity(self):
        '''Eccentricity vector'''
        return (np.cross(self.v, self.h) / self.mu) - (self.r / np.linalg.norm(self.r))

    def vector_n(self):
        return np.cross(np.array([0, 0, 1]), self.h)

    def semi_major_axis(self):
        '''Semi-major axis'''
        return - self.mu / (2 * self.energy())

    def energy(self):
        '''Specific energy'''
        return np.linalg.norm(self.v)**2 / 2 - self.mu / np.linalg.norm(self.r)

    def inclination(self):
        '''Inclination in radians'''
        return np.arccos(np.dot(self.h / np.linalg.norm(self.h) , np.array([0, 0, 1])))

    def ascending_node(self):
        '''Ascending node in radians'''
        Omega = np.arccos(np.dot(np.array([1, 0, 0]), self.vector_n() / np.linalg.norm(self.vector_n())))
        
        if self.vector_n()[1] < 0:
            Omega = 2 * np.pi - Omega

        return Omega

    def argument_of_periapsis(self):
        '''Argument of periapsis in radians'''
        omega = np.arccos(np.dot(self.vector_n(),self.e) /
         np.linalg.norm(self.vector_n())/ np.linalg.norm(self.e))
        
        if self.e[2] < 0:
            omega = 2 * np.pi - omega

        return omega

    def true_anomaly(self):
        '''True anomaly in radians'''
        self.nu = np.arccos(np.dot(self.r, self.e) / np.linalg.norm(self.r) / np.linalg.norm(self.e))
        
        if np.dot(self.r, self.v) < 0:
            self.nu = 2 * np.pi - self.nu

        return self.nu

    def period(self):
        '''Period in seconds'''
        return 2 * np.pi / self.n

    def mean_motion(self):
        '''Mean motion in radians/second'''
        if self.a < 0:
            return - np.sqrt(-self.mu / self.a**3)
        return np.sqrt(self.mu / self.a**3)

    def mean_anomaly(self):
        '''Mean anomaly in radians'''

        e = np.linalg.norm(self.eccentricity())
        self.M = self.nu + e * np.sin(self.nu)
 
        if self.M <= 0:
            self.M += 2 * np.pi
        elif self.M >= 2 * np.pi:
            self.M -= 2 * np.pi

        return self.M


    def eccentric_anomaly(self):
        '''Eccentric anomaly in radians'''
        return np.arctan2(np.dot(self.r, self.v) / (self.e * np.sqrt(self.mu * self.a)), 
        (self.a - np.linalg.norm(self.r) / (self.a*self.e)))

    