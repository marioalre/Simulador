import numpy as np
from src.Orbit import Orbit 
import perturbations as pert

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
            R_osc, V_osc = self.r0v02rv(Rp, Vp, t-t0)

            r_osc = np.linalg.norm(R_osc)
            v_osc = np.linalg.norm(V_osc)
            
            dr = R_osc - Rp # Vector perturbation of position
            dv = V_osc - Vp # Vector perturbation of velocity

            # Compute the perturbation
            epsilon = R_osc * dr / r_osc**2

            f = 1/epsilon * (1 - 1/(1-2*epsilon)**(3/2))

            ad = self.acceleration_perturbation()

            d2dr = ad + self.mu / r_osc**3 * (f/epsilon*Rp - dr)

            if (np.linalg.norm(dr)/np.linalg.norm(Rp)) > 0.01:
                R_osc = Rp
                V_osc = Vp
            else:
                R_osc = Rp + dr
                V_osc = Vp + dv
                t += dt
            
    def acceleration_perturbation(self):
        pass

    def acc_third_body(self, body_3rd):
        '''Compute the acceleration due to a third body
        Parameters
        ----------
        body_3rd : CelestialBody
            Third body
        Returns
        ------- 
        acc : ndarray
            Acceleration due to third body'''

        pass

