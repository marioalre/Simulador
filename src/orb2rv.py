# Class to transform orbital parameters to position and velocity vectors

import numpy as np
from src.Orbit import Orbit

class Orb2rv(Orbit):
    '''Orbital parameters to position and velocity vectors'''

    def __init__(self, h, e, Omega, omega, i, nu, body):
        '''Initialize the class'''
        self.h = h
        self.e = e
        self.Omega = Omega
        self.omega = omega
        self.i = i
        self.nu = nu
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


    





