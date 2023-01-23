import numpy as np
from src.orb2rv import Orb2rv
from src.rv2orb import Rv2orb

class KeplerPropagator:

    def __init__(self, body):
        self.mu = body.mu
        self.radius = body.radius
        self.body = body

    def c2c3fromchi(self, chi):
        '''Compute c2 and c3 from chi
        A.1
        Parameters
        ----------
        chi : float
            Chi parameter
        Returns
        -------
        c2 : float
            c2 parameter
        c3 : float
            c3 parameter
        '''

        # comprobar si el cubo va dentro o fuera de la raiz

        if chi > 1e-6:
            c2 = (1 - np.cos(np.sqrt(chi))) / chi
            c3 = (np.sqrt(chi) - np.sin(np.sqrt(chi))) / np.sqrt(chi)**3

        else:
            if chi < -1e-6:
                c2 = (1 - np.cosh(np.sqrt(-chi))) / chi
                c3 = (np.sinh(np.sqrt(-chi)) - np.sqrt(-chi)) / np.sqrt(-chi)**3
            else:
                c2 = 1 / 2
                c3 = 1 / 6

        return c2, c3


    def KepEqtnE(self, M, e):
        '''Solve Kepler's equation for eccentric anomaly
        A.2
        Parameters
        ----------
        M : float
            Mean anomaly in radians
        e : float
            Eccentricity
        Returns
        -------
        E : float
            Eccentric anomaly in radians
        '''

        if -np.pi < M < 0 or np.pi < M:
            E = M - e
        else:
            E = M + e

        while abs(E - e * np.sin(E) - M) > 1e-10:
            E = E - (E - e * np.sin(E) - M) / (1 - e * np.cos(E))

        return E

    def KepEqtnP(self, dt, p):
        '''Solve Kepler's equation for parabolic anomaly
        A.3
        Parameters
        ----------
        dt : float
            Time since periapsis in seconds
        p : float

        Returns
        -------
        B : float
            Parabolic anomaly in radians
        '''

        n_p = 2 * np.sqrt(self.mu / p ** 3)

        # Solve Backer equation
        cot2s = 3 / 2 * dt * n_p

        s = 0.5 * np.arctan(1 / cot2s)

        w = np.arctan(np.cbrt(np.tan(s)))

        return 2 * 1 / np.tan(2 * w)

    def KepEqtnH(self, M, e):
        '''Solve Kepler's equation for hyperbolic anomaly
        A.4
        Parameters
        ----------
        M : float
            Mean anomaly in radians
        e : float
            Eccentricity
        Returns
        -------
        H : float
            Hyperbolic anomaly in radians
        '''

        if e < 1.6:
            if -np.pi < M < 0 or np.pi < M:
                H = M - e
            else:
                H = M + e
        else:
            if e < 3.6 and np.abs(M) > np.pi:
                H = M - np.sign(M) * e
            else:
                H = M / (e - 1)


        while abs(M + H - e * np.sinh(H)) > 1e-10:
            H = H + (M + H - e * np.sinh(H)) / (e * np.cosh(H) - 1)

        return H

    def nu2anomaly(self, e, nu):
        '''Compute true anomaly from eccentric anomaly
        A.5
        Parameters
        ----------
        e : float
            Eccentricity
        nu : float
            True anomaly in radians
        Returns
        ------- 
        E : float
            Eccentric anomaly in radians
        '''

        if e < 1:
            sinE = np.sqrt(1 - e ** 2) * np.sin(nu) / (1 + e * np.cos(nu))
            cosE = (e + np.cos(nu)) / (1 + e * np.cos(nu))

            E = np.arcsin(sinE)

        elif e == 1:
            E = np.tan(nu/2)

        else:
            sinhH = np.sqrt(e ** 2 - 1) * np.sin(nu) / (1 + e * np.cos(nu))
            coshH = (e + np.cos(nu)) / (1 + e * np.cos(nu))

            E = np.arcsinh(sinhH)

        return E

    def anomaly2nu(self, e, E, extra=None):
        '''Compute true anomaly from eccentric anomaly
        A.6
        Parameters
        ----------
        e : float
            Eccentricity
        E : float
            Eccentric anomaly in radians
        entra : list
            Extra parameters (For parabolic and hyperbolic orbits)
            - p
            - r
        Returns
        -------
        nu : float
            True anomaly in radians
        '''
        if extra is not None:
            p = extra[0]
            r = extra[1]

        if e < 1:
            sinnu = np.sqrt(1 - e ** 2) * np.sin(E) / (1 - e * np.cos(E))
            cosnu = (np.cos(E) - e) / (1 - e * np.cos(E))

            nu = np.arcsin(sinnu)

        elif e == 1:
            nu = np.arcsin(p*E / r)

        else:
            sinnu = -np.sqrt(e ** 2 - 1) * np.sinh(E) / (e - np.cosh(E))
            cosnu = (np.cosh(E) - e) / (1 - e*np.cosh(E))

            nu = np.arcsin(sinnu)

        return nu

    def keplerCOE(self, R0, V0, dt):
        '''Compute the classical orbital elements from the initial state and time
        A.7
        Parameters
        ----------
        r0 : numpy.array
            Initial position vector in km
        v0 : numpy.array
            Initial velocity vector in km/s
        dt : float
            Time since periapsis in seconds
        Returns
        -------
        r : numpy.array
            Position vector in km
        v : numpy.array
            Velocity vector in km/s
        '''
        r0 = np.linalg.norm(R0)
        v0 = np.linalg.norm(V0)
        orb = self.rv2coe(R0, V0, dt, special=None)

        a = orb[0]
        e = orb[1]
        i = orb[2]
        Omega = orb[3]
        omega = orb[4]
        nu = orb[5]
        u = orb[6]
        lambda_true = orb[7]
        w_true = orb[8]

        if e != 0:
            E = self.nu2anomaly(e, nu)
        else:
            E0 = u
            E0 = lambda_true

        # Mean velocity
        # n = np.sqrt(self.mu / a ** 3)
        n = np.sqrt(self.mu / a ** 3)

        if e < 1.0:
            M = E - e * np.sin(E)
            M = M + n*dt

            E = self.KepEqtnE(M, e)
        elif e == 1:
            H = np.cross(R0, V0)
            p = np.linalg.norm(H) ** 2 / self.mu
            M = E + E**3 / 3

            E = self.KepEqtnP(dt, p)
        else:
            M = e*np.sinh(E) - E
            M = M + n*dt

            E = self.KepEqtnH(M, e)

        if e != 0:
            self.nu = self.anomaly2nu(e, E)
        else:
            u = E
            lambda_true = E

        # Compute the position and velocity vectors
        self.coe2rv(a, e, i, Omega, omega, nu)

    def kepler(self, R, V, dt):
        '''Compute the classical orbital elements from the initial state and time
        universal variable
        Parameters
        ----------
        r : numpy.array
            Position vector in km
        v : numpy.array
            Velocity vector in km/s
        dt : float
            Time since periapsis in seconds
        Returns
        -------
        r : numpy.array
            Position vector in km
        v : numpy.array
            Velocity vector in km/s
        '''
        if dt == 0:
            return R, V

        r = np.linalg.norm(R)
        v = np.linalg.norm(V)

        eta = v**2 / 2 - self.mu / r
        alpha = -v**2 / self.mu + 2 / r

        if alpha > 0.000001:
            chi = np.sqrt(self.mu) * dt * alpha
        elif abs(alpha) < 0.000001:
            H = np.cross(R, V)
            h = np.linalg.norm(H)

            p = h**2 / self.mu
            s = 0.5 * np.arctan(1 /(3 * np.sqrt(self.mu/p**3) *dt))
            w = np.arctan(np.cbrt(np.tan(s)))

            chi = np.sqrt(p) * 2 * 1 / np.tan(2*w)
        elif alpha < -0.000001:
            a = 1 / alpha
            chi = np.sign(dt) * np.sqrt(-a) * np.log(-2 * self.mu * alpha * dt / (np.dot(R, V) + np.sign(dt) * np.sqrt(-self.mu * a) * (1 - r * alpha)))

        ratio = 1
        max_iter = 1000
        iter = 0

        while np.abs(ratio) > 1e-10 and max_iter > iter:
            phi = chi**2 * alpha

            c2, c3 = self.c2c3fromchi(phi)

            Q1 = np.sqrt(self.mu)*dt-chi**3*c3-np.dot(R,V)/np.sqrt(self.mu)*chi**2*c2-r*chi*(1-phi*c3)
            Q2 = chi**2*c2+np.dot(R,V)/np.sqrt(self.mu)*chi*(1-phi*c3)+r*(1-phi*c2)

            ratio = Q1/Q2
            chi = chi + ratio
            iter += 1

        f = 1 - chi**2 * c2 / r
        g = dt - chi**3 * c3 / np.sqrt(self.mu)

        dg = 1 - chi**2 * c2 /Q2
        df = np.sqrt(self.mu) / (Q2 * r) * chi * (phi * c3 - 1)

        R_final = f * R + g * V
        V_final = df * R + dg * V

        return R_final, V_final

    def rv2coe(self, R, V, special=None):
        '''Compute the classical orbital elements from the position and velocity
        A.8
        Parameters
        ----------
        R : numpy.array
            Position vector in km
        V : numpy.array
            Velocity vector in km/s
        special : str
            Special case (For parabolic and hyperbolic orbits)
        Returns
        -------
        orb : OrbitalElements
            a : float
                Semi-major axis in km
            e : float
                Eccentricity
            i : float
                Inclination in radians
            Omega : float
                Right ascension of the ascending node in radians
            omega : float
                Argument of periapsis in radians
            nu : float
                True anomaly in radians
        '''
        r = np.linalg.norm(R)
        v = np.linalg.norm(V)

        H = np.cross(R, V)
        h = np.linalg.norm(H)

        N = np.cross([0, 0, 1], H)
        n = np.linalg.norm(N)

        ecc = (np.dot((v**2-self.mu/r),R) - np.dot(np.dot(R, V),V)) / self.mu
        e = np.linalg.norm(ecc)

        eta = v**2/2 - self.mu/r

        if e != 0:
            a = -self.mu / (2 * eta)
            p = a * (1 - e**2)
        else:
            p = h**2 / self.mu
            a = np.inf

        i = np.arccos(H[2] / h)

        Omega = np.arccos(N[0] / n)

        if N[1] < 0:
            Omega = 2 * np.pi - Omega

        omega = np.arccos(np.dot(N, ecc) / (n * e))
        if ecc[2] < 0:
            omega = 2 * np.pi - omega

        nu = np.arccos(np.dot(ecc, R) / (e * r))
        if np.dot(R, V) < 0:
            nu = 2 * np.pi - nu

        orb = [a, e, i, Omega, omega, nu]

        # Special cases
        # Ecliptical orbit
        w_true = np.arccos(ecc[0] / e)
        if ecc[1] < 0:
            w_true = 2 * np.pi - w_true
            orb.append(w_true)
                

        # Circular orbit

        u = np.arccos(np.dot(N, R) / (n * r))
        if R[2] < 0:
            u = 2 * np.pi - u
            orb.append(u)


        # Equatorial orbit
        lambda_true = np.arccos(R[0] / r)
        if R[1] < 0:
            lambda_true = 2 * np.pi - lambda_true
            orb.append(lambda_true)

        if special == 'equatorial':
            print('Equatorial orbit')
        elif special == 'circular':
            print('Circular orbit')
        elif special == 'ecliptic':
            print('Ecliptical orbit')

        return orb

    def coe2rv(self, p, e, i, Omega, omega, nu, extra=None, special=None):
        '''Compute the position and velocity vectors from the classical orbital elements
        A.9
        Parameters
        ----------
        p : float
            Semi-l
        e : float
            Eccentricity
        i : float
            Inclination in radians
        Omega : float
            Right ascension of the ascending node in radians
        omega : float
            Argument of periapsis in radians
        nu : float
            True anomaly in radians
        extra : list
            Special case (For parabolic and hyperbolic orbits)
            - u : float
                True longitude in radians
            - lambda_true : float
                True longitude in radians
            - w_true : float
                True longitude in radians
        special : str
            Special case (For parabolic and hyperbolic orbits)
            circular, equatorial, ecliptic
        Returns
        -------
        r : numpy.array
            Position vector in km
        v : numpy.array
            Velocity vector in km/s
        '''
        if extra is not None:
            if special == 'equatorial':
                omega = 0
                Omega = 0
                nu = extra[1]
            elif special == 'circular':
                omega = 0
                nu = extra[0]
            elif special == 'ecliptic':
                Omega = 0
                omega = extra[2]

        r = np.array([p * np.cos(nu) / (1 + e * np.cos(nu)), p * np.sin(nu) / (1 + e * np.cos(nu)), 0]) 
        v = np.array([-np.sqrt(self.mu / p) * np.sin(nu), np.sqrt(self.mu / p) * (e + np.cos(nu)), 0])

        # Rotate to the correct frame

        R1 = np.array([[np.cos(Omega), np.sin(Omega), 0],
                       [-np.sin(Omega), np.cos(Omega), 0],
                       [0, 0, 1]])

        R2 = np.array([[1, 0, 0],
                       [0, np.cos(i), np.sin(i)],
                       [0, -np.sin(i), np.cos(i)]])

        R3 = np.array([[np.cos(omega), np.sin(omega), 0],
                       [-np.sin(omega), np.cos(omega), 0],
                       [0, 0, 1]])

        rotation_matrix = np.transpose(R1) @ np.transpose(R2) @ np.transpose(R3)

        r = rotation_matrix @ r
        v = rotation_matrix @ v

        return r, v


if __name__ == "__main__":

    from CelestialBodies import CelestialBodies

    # Earth
    earth = CelestialBodies()
    earth.earth()

    kepler = KeplerPropagator(earth)

    orb = kepler.rv2coe([6524.834, 6862.875, 6448.296], [4.901327, 5.533756, -1.976341])
    
    print(orb)

    r, v = kepler.coe2rv(11067.79, 0.83285, 87.87 * np.pi/180, 227.89 * np.pi/180, 53.38 * np.pi/180, 92.335 * np.pi/180)
    print(r)
    print(v)

    E = kepler.KepEqtnE(235.4*np.pi/180, 0.4)
    print(E)
    B = kepler.KepEqtnP(53.7874*60, 25512)
    print(B)
    H = kepler.KepEqtnH(235.4*np.pi/180, 2.4)
    print(H)


    R = np.array([1131.340, -2282.343, 6672.423])
    V = np.array([-5.64305, 4.30333, 2.42879])
    dt = 40 * 60

    r, v = kepler.kepler(R, V, dt)
    print(r)
    print(v)



