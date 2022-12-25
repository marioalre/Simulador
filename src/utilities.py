import numpy as np

class Utilities:
    def __init__(self, body=None):
        if body is not None:
            self.mu = body.mu
        else:
            self.mu = 398600.4418 # km^3/s^2

    def findTOF(self, r0, r1, p):
        '''
        Find time of flight
        Parameters
        ----------
        r0 : numpy array
            Initial position vector
        r1 : numpy array
            Final position vector
        p : float
            Semi-latus
        Returns
        -------
        TOF : float
            Time of flight


        ATTENTION: This function is not working properly
        '''

        cosdv = np.dot(r0, r1) / (np.linalg.norm(r0) * np.linalg.norm(r1))

        k = np.linalg.norm(r0) * np.linalg.norm(r1) * (1 - cosdv)

        l = np.linalg.norm(r0) + np.linalg.norm(r1)

        m = np.linalg.norm(r0) * np.linalg.norm(r1) * (1 + cosdv)

        a = (m * k * p) / ((2*m-l**2)*p**2 + 2*k*l*p - k**2)

        f = 1- np.linalg.norm(r1) / p * (1 - cosdv)

        g = np.linalg.norm(r0) * np.linalg.norm(r1) * np.sin(np.arccos(cosdv)) / np.sqrt(p * self.mu)

        if a > 0:
            # Elliptic
            df = np.sqrt(self.mu / p) * np.tan(np.arccos(cosdv)/2) *((1 - cosdv) / p - 1 / np.linalg.norm(r0) - 1 / np.linalg.norm(r1)) 

            cosde = 1 - np.linalg.norm(r0) / a * (1 - f)

            sinde = - np.linalg.norm(r0) * np.linalg.norm(r1) * df / (np.sqrt(p * self.mu))

            TOF = g + np.sqrt(a**3 / self.mu) * (np.arccos(cosde) - sinde)

        elif a == np.inf:
            # Parabolic
            c = np.sqrt(np.linalg.norm(r0)**2 + np.linalg.norm(r1)**2 - 2 * np.linalg.norm(r0) * np.linalg.norm(r1) * cosdv)

            s = (np.linalg.norm(r0) + np.linalg.norm(r1) + c) / 2

            TOF = 2/3 * np.sqrt(s**3 / (2*self.mu)) * (1 - ((s - c) / s)**(3/2))

        else:
            # Hyperbolic
            cosdh = 1 + np.linalg.norm(r0) / a * (f - 1)

            TOF = g + np.sqrt(-a**3 / self.mu) * (np.sinh(np.arccosh(cosdh)) - np.arccosh(cosdh))

        return TOF


    def Gibbs(self, R0, R1, R2):
        '''
        Gibs method to find v2. the position vectors are given in the geometric frame
        Parameters
        ----------
        r0 : numpy array
            Initial position vector
        r1 : numpy array
            Intermediate position vector
        r2 : numpy array
            Final position vector
        Returns
        ----------
        v2 : numpy array
            Final velocity vector'''

        # Parameters
        irr = 0
        tol = 1e-10

        R0 = np.array(R0)
        R1 = np.array(R1)
        R2 = np.array(R2)
        
        # Normalize position vectors
        r0 = np.linalg.norm(R0)
        r1 = np.linalg.norm(R1)
        r2 = np.linalg.norm(R2)

        # Cross product
        r0xr1 = np.cross(R0, R1)
        r1xr2 = np.cross(R1, R2)
        r2xr0 = np.cross(R2, R0)

        # Coplanarity
        if np.abs(np.dot(R0, r1xr2) / r0 / np.linalg.norm(r1xr2)) < 1e-10:
            print('Coplanarity condition not satisfied')
            irr = 1

        N = r0 * r1xr2 + r1 * r2xr0 + r2 * r0xr1

        D = r0xr1 + r1xr2 + r2xr0

        S = R0 * (r1 - r2) + R1 * (r2 - r0) + R2 * (r0 - r1)

        return np.sqrt(self.mu / np.linalg.norm(N) / np.linalg.norm(D)) * (np.cross(D, R1) / r1 + S), irr




