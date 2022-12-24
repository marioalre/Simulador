import numpy as np

class Orbit():
    def __init__(self, body):
        self.mu = body.mu
        self.radius = body.radius

    def kepler_elliptic(self, M, e):
        '''Kepler equation for elliptic orbits by Newton-Raphson method
        Parameters
        ----------
        M : float 
            Mean anomaly 
        e : float
            Eccentricity
        '''

        if -np.pi<M<0 or np.pi<M:
            E = M - e
        else:
            E = M + e

        while abs(E - e * np.sin(E) - M) > 1e-10:
            E = E - (E - e * np.sin(E) - M) / (1 - e * np.cos(E))

        return E

    def kepler_parabolic(self, dt, p):
        '''Kepler equation for parabolic orbits by Newton method

        B^3 /3 + B - np * dt = 0

        Parameters
        ----------
        dt : float
            Time since periapsis
        p : float

        e : float
            Eccentricity
        Returns
        -------
        B : float
            Parabolic anomaly
        '''

        n_p = 2* np.sqrt(p**3 / self.mu)

        # Initial guess
        B = 0

        while abs(B**3 / 3 + B - n_p * dt) > 1e-10:
            B = B - (B**3 / 3 + B - n_p * dt) / (B**2 + 1)
            print(B)

        # No converge, corregir
        return B

    def kepler_hyperbolic(self, M, e):
        '''Kepler equation for hyperbolic orbits by Newton-Raphson method

        Parameters
        ----------
        M : float
            Mean anomaly
        e : float
            Eccentricity
        Returns 
        -------
        F : float
            Hyperbolic anomaly
        '''

        if e < 1.6:
            if -np.pi<M<0 or np.pi<M:
                F = M - e
            else:
                F = M + e
        else:
            if e<3.6 and np.abs(M)>np.pi:
                F = M - np.sign(M) * e
            else:
                F = M / (e - 1)

        while abs(F - e * np.sinh(F) - M) > 1e-10:
            F = F - (F - e * np.sinh(F) - M) / (e * np.cosh(F) - 1)

        return F


        

    

    