# pasar de v y r a los parametros orbitales
import numpy as np

class Orbit:
    def __init__(self, r, v, mu, t0, t):
        self.r = r   # position vector ECI 
        self.v = v   # velocity vector ECI
        self.mu = mu
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

    def semi_major_axis(self):
        '''Semi-major axis'''
        return 1 / (2/np.linalg.norm(self.r) - np.linalg.norm(self.v)**2 / self.mu)

    def inclination(self):
        '''Inclination'''
        return np.arccos(self.h[2] / np.linalg.norm(self.h))

    def ascending_node(self):
        '''Ascending node'''
        return np.arctan2(self.h[0], -self.h[1])

    def argument_of_periapsis(self):
        '''Argument of periapsis'''
        return np.arctan2(np.dot(np.cross([0, 0, 1], self.h), self.e), np.dot(self.e, self.h))

    def true_anomaly(self):
        '''True anomaly'''
        return np.arctan2(np.linalg.norm(np.cross(self.h, self.e)), np.dot(self.e, self.r))

    def mean_anomaly(self):
        '''Mean anomaly'''
        E = 2 * np.arctan(np.sqrt((1 - np.linalg.norm(self.e)) / (1 + np.linalg.norm(self.e))) * np.tan(self.nu / 2))
        return E - self.e * np.sin(E)

    def mean_motion(self): 
        '''Mean motion'''
        return np.sqrt(self.mu / self.a**3)

    def period(self):
        '''Period'''
        return 2 * np.pi / self.n

    def position(self):
        '''Position vector at time t'''
        M = self.mean_anomaly() + self.n * self.dt
        E = self.kepler(M)
        nu = 2 * np.arctan(np.sqrt((1 + np.linalg.norm(self.e)) / (1 - np.linalg.norm(self.e))) * np.tan(E / 2))
        r = self.a * (1 - self.e**2) / (1 + self.e * np.cos(nu))
        return r * np.array([np.cos(nu), np.sin(nu), 0])

    def kepler(self, M):
        '''Solve Kepler's equation iteratively'''
        E = M
        for i in range(100):
            E = M + self.e * np.sin(E)
        return E

    def velocity(self):
        '''Velocity vector at time t'''
        M = self.mean_anomaly() + self.n * self.dt
        E = self.kepler(M)
        nu = 2 * np.arctan(np.sqrt((1 + np.linalg.norm(self.e)) / (1 - np.linalg.norm(self.e))) * np.tan(E / 2))
        r = self.a * (1 - self.e**2) / (1 + self.e * np.cos(nu))
        v = np.sqrt(self.mu * self.a) / r * np.array([-np.sin(nu), np.sqrt(1 - self.e**2) * np.cos(nu), 0])
        return v

