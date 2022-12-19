# Esta clase se encarga de pasar de los vectores de posición y velocidad a los vectores de posición y velocidad
# en el sistema de referencia de la Tierra a los parametros de la orbita (a, e, i, O, w, M) para posteriormente
# pasarlos a los vectores de posición y velocidad en el sistema de referencia de la Tierra en un instante de 
# tiempo posterior.
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
        nu = np.arccos(np.dot(self.r, self.e) / np.linalg.norm(self.r) / np.linalg.norm(self.e))
        
        if np.dot(self.r, self.v) < 0:
            nu = 2 * np.pi - nu

        return nu

    def period(self):
        '''Period in seconds'''
        return 2 * np.pi / self.n
    def eccentric_anomaly(self):
        '''Eccentric anomaly in radians'''
        return 2 * np.arctan(np.sqrt((1 - self.e) / (1 + self.e)) * np.tan(self.nu / 2))

    def mean_anomaly(self):
        '''Mean anomaly in radians'''
        E = self.eccentric_anomaly()
        return E - self.e * np.sin(E)

    def mean_motion(self):
        '''Mean motion in radians per second'''
        return np.sqrt(self.mu / self.a**3)

    def orbital_elements(self):
        '''Orbital elements'''
        return self.a, self.e, self.i, self.Omega, self.omega, self.nu

    
    # Cambios de coordenadas

    def eci2pqw(self, t):
        '''Transforms ECI to PQW'''
        r_pqw = np.dot(self.rotation_matrix(self.Omega, 'z'), self.position())
        r_pqw = np.dot(self.rotation_matrix(self.i, 'x'), r_pqw)
        r_pqw = np.dot(self.rotation_matrix(self.omega, 'z'), r_pqw)

        v_pqw = np.dot(self.rotation_matrix(self.Omega, 'z'), self.velocity())
        v_pqw = np.dot(self.rotation_matrix(self.i, 'x'), v_pqw)
        v_pqw = np.dot(self.rotation_matrix(self.omega, 'z'), v_pqw)

        return r_pqw, v_pqw

    def pqw2eci(self, t):
        '''Transforms PQW to ECI'''
        r_eci = np.dot(self.rotation_matrix(-self.omega, 'z'), self.position())
        r_eci = np.dot(self.rotation_matrix(-self.i, 'x'), r_eci)
        r_eci = np.dot(self.rotation_matrix(-self.Omega, 'z'), r_eci)

        v_eci = np.dot(self.rotation_matrix(-self.omega, 'z'), self.velocity())
        v_eci = np.dot(self.rotation_matrix(-self.i, 'x'), v_eci)
        v_eci = np.dot(self.rotation_matrix(-self.Omega, 'z'), v_eci)

        return r_eci, v_eci

    def rotation_matrix(self, angle, axis):
        '''Rotation matrix'''
        if axis == 'x':
            return np.array([[1, 0, 0], [0, np.cos(angle), -np.sin(angle)], [0, np.sin(angle), np.cos(angle)]])
        elif axis == 'y':
            return np.array([[np.cos(angle), 0, np.sin(angle)], [0, 1, 0], [-np.sin(angle), 0, np.cos(angle)]])
        elif axis == 'z':
            return np.array([[np.cos(angle), -np.sin(angle), 0], [np.sin(angle), np.cos(angle), 0], [0, 0, 1]])

























