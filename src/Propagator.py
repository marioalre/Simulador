import numpy as np
from src.Orbit import Orbit 

class Propagator(Orbit):
    '''Class for perturbations propagation and computation''' 

    def __init__(self, body, R0, V0):
        self.mu = body.mu
        self.radius = body.radius
        self.R0 = np.array(R0)
        self.V0 = np.array(V0)
        
    def Encke(self, t0, tf, dt):
        '''Encke propagation method
        Parameters
        ----------
        t0 : float
            Initial time in seconds
        tf : float
            Final time in seconds
        dt : float
            Time step in seconds
        Returns
        -------
        r : ndarray
            Position vector in km
        v : ndarray
            Velocity vector in km/s
        '''
        # Initial conditions
        r0 = np.linalg.norm(self.R0)
        v0 = np.linalg.norm(self.V0)

        dr = np.array([0, 0, 0])
        epsilon = 0
        f = 0
        Rp = self.R0
        Vp = self.V0

        t = t0

        while t < tf:
            R_osc, V_osc = self.r0v02rv(Rp, Vp, dt)

            r_osc = np.linalg.norm(R_osc)
            v_osc = np.linalg.norm(V_osc)
            
            # Compute the perturbation
            epsilon = 1
            pass

    def a_third_body(self, body):
        '''Compute the perturbation due to a third body
        Parameters
        ----------
        body : Body
            Third body
        Returns
        -------
        dr : ndarray
            Perturbation in position vector
        dv : ndarray
            Perturbation in velocity vector
        '''
        pass



