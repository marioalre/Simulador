# Import modules

import numpy as np
from scipy.integrate import RK45
# Define class
class nbody:

    def __init__(self, bodies, dt):

        self.mu = np.zeros(len(bodies))
        self.radius = np.zeros(len(bodies))


        for ii, body in enumerate(bodies):
            self.radius[ii] = body.radius
            self.mu[ii] = body.mu
    

    def acceleration(self, t, f):
        '''Compute the acceleration of a n-body system
        Parameters
        ----------
        t : float
            Time in seconds
        f : ndarray
            State vector
        Returns
        -------
        a : ndarray
            Acceleration vector in km/s^2
        '''

        # If f has not a pair number of elements, raise an error
        if len(f) % 2 != 0:
            raise ValueError('f must have a pair number of elements')
        else:
            r = f[:len(f)//2, :]
            v = f[len(f)//2:, :]

        a = np.zeros((len(r), 3))

        for ii in range(len(r)-1):
            for jj in range(len(r)-1):
                if ii != jj and ii < jj:
                    r_ij = r[ii] - r[jj]
                    r_ij_norm = np.linalg.norm(r_ij)
                    a[ii] += -self.mu[jj] * r_ij / r_ij_norm**3

        return np.concatenate((v, a, [0 , 0, 0]), axis=0)



    def propagate(self, r0, v0, a):
        '''Propagate the n-body system
        Parameters
        ----------
        v : ndarray
            Velocity vector in km/s
            n x 3 array
        a : ndarray
            Acceleration vector in km/s^2
            n x 3 array
        Returns
        -------
        v : ndarray
            Velocity vector in km/s
            n x 3 array
        '''

        Rg = np.dot(self.mu, r0) / np.sum(self.mu)
        Vg = np.dot(self.mu, v0) / np.sum(self.mu)

        f0 = np.concatenate((r0, Rg, v0, Vg), axis=0)

        t0 = 0
        tf = 3600

        t , f = RK45(self.acceleration, t0, f0, tf)
        


# Path: scripts\Propagator.py
if __name__ == '__main__':
    
    pass