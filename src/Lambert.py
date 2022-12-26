import numpy as np 
from src.Orbit import Orbit

class Lambert(Orbit):
    '''Lambert's problem solver class for two point transfers in space with a given body's gravitational parameter mu
    R0 and R1 are the initial and final position vectors'''
    def __init__(self, R0, R1, body):
        self.R0 = np.array(R0)
        self.R1 = np.array(R1)
        self.mu = body.mu

    def minimum_energy(self):
        '''Minimum energy transfer between two points in space'''

        r0 = np.linalg.norm(self.R0)
        r1 = np.linalg.norm(self.R1)

        cosdv = np.dot(self.R0, self.R1) / (r0 * r1)
        c = np.sqrt(r0**2 + r1**2 - 2 * r0 * r1 * cosdv)
        s = (r0 + r1 + c) / 2
        a_min = s / 2
        p_min = r0 * r1 / c * (1 - cosdv)
        e_min = np.sqrt(1 - 2*p_min / s)
        alpha_e = np.pi
        sinb_2 = np.sqrt((s-c)/s)

        t_min_a_min = np.sqrt(a_min**3 / self.mu) * (alpha_e - (2*np.arcsin(sinb_2) - sinb_2))

        t_min_abs = 1/3 * np.sqrt(2/self.mu) * (s**(3/2) - (s -c)**(3/2))

        v0  =  np.sqrt(self.mu * p_min) / (r0*r1*np.sin(np.arccos(cosdv))) * (self.R1 - (1 - r1/p_min * (1 - cosdv)) * self.R0)

        return v0

    def Gauss(self, dt, tm):
        pass

    def universal(self, dt, str):
        '''Universal variable solution for Lambert's problem
        Parameters
        ----------
        dt : float
            Time of flight in seconds
        str : string
            Transfer type: 'pro' or 'retro'
        Returns
        -------
        v0 : float
            Initial velocity vector in km/s
        v1 : float
            Final velocity vector in km/s
        '''
        r0 = np.linalg.norm(self.R0)
        r1 = np.linalg.norm(self.R1)

        cross01 = np.cross(self.R0, self.R1)

        theta = np.arccos(np.dot(self.R0, self.R1) / (r0 * r1))

        if str == 'pro':
            if cross01[2] < 0:
                theta = 2 * np.pi - theta
        elif str == 'retro':
            if cross01[2] >= 0:
                theta = 2 * np.pi - theta
        else:
            print('We will assume a prograde transfer')
            
        A = np.sin(theta) * np.sqrt(r0 * r1 / (1 - np.cos(theta)))

        # Strarting guess for z
        z = -100
        while self.f(z, A, dt) < 0:           
            z += 0.1

        # Newton-Raphson method
        tol = 1e-10
        nMax = 1000

        ratio = 1
        iter = 0

        while (np.abs(ratio) > tol) and (iter < nMax):
            iter += 1
            ratio = self.f(z, A, dt) / self.df(z, A)
            z -= ratio

        if iter == nMax:
            print('Maximum number of iterations reached')

        f = 1 - self.y(z, A) / r0
        g = A * np.sqrt(self.y(z, A) / self.mu)

        gdot = 1 - self.y(z, A) / r1

        V0 = (1/g) * (self.R1 - f * self.R0)
        V1 = (1/g) * (gdot * self.R1 - self.R0)

        return V0.real, V1.real

    def y(self, z, A):
        r1 = np.linalg.norm(self.R1)
        r0 = np.linalg.norm(self.R0)
        S = complex(self.stumpff_s(z))
        C = complex(self.stumpff_c(z))

        return r0 + r1 + A * (z * S - 1) / np.sqrt(C)
    
    def f(self, z, A, dt):
        S = complex(self.stumpff_s(z))
        C = complex(self.stumpff_c(z))
        y = complex(self.y(z, A))
        return (y / C)**(3/2) * S + A * np.sqrt(y) - np.sqrt(self.mu) * dt

    def df(self, z, A):

        S = complex(self.stumpff_s(z))
        C = complex(self.stumpff_c(z))
        y = complex(self.y(z, A))

        if z == 0:
            return np.sqrt(2)/40*y**(3/2) + 1/8 * A * (np.sqrt(y) + A *np.sqrt(1/(2*y)))
        else:
            return (y / C)**(3/2) * (1/(2*z) * (C - 3/2* S / C) + 3/4 *S**2 / C) + A / 8 * (3* S / C * np.sqrt(y) + A * np.sqrt(C/y))