import numpy as np 

class Lambert():
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

    def universal(self, dt, alpha):
        pass
