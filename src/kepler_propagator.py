import numpy as np
import pandas as pd
import math
from src.utilities import Utilities

class KeplerPropagator:

    def __init__(self, body):
        self.mu = body.mu
        self.radius = body.radius
        self.body = body
        self.small = 1e-10
        self.iterations = 50

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
        M : float
            Mean anomaly in radians
        '''

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

        '''


        E = 999999.9
        m = 999999.9
        small = 0.00000001

        # circular
        if abs(e) < small:
            m = nu
            E = nu
        else:
            # elliptical
            if e < 1.0 - small:
                sine = (np.sqrt(1.0 - e**2) * np.sin(nu)) / (1.0 + e*np.cos(nu))
                cose = (e + np.cos(nu)) / (1.0 + e*np.cos(nu))
                E = np.arctan2(sine, cose)
                m = E - e*np.sin(E)
            else:
                # hyperbolic
                if e > 1.0 + small:
                    if e > 1.0 and abs(nu) + 0.00001 < math.pi - np.arccos(1.0 / e):
                        sine = (np.sqrt(e**2 - 1.0) * np.sin(nu)) / (1.0 + e*np.cos(nu))
                        E = math.asinh(sine)
                        m = e*np.sinh(E) - E
                else:
                    # parabolic
                    if abs(nu) < 168.0*math.pi/180.0:
                        E = np.tan(nu*0.5)
                        m = E + (E**3)/3.0

        if e < 1.0:
            m = m % (2.0*np.pi)
            if m < 0.0:
                m = m + 2.0*np.pi
            E = E % (2.0*np.pi)

        return E, m

    def anomaly2nu2(self, e, E, extra=None):
        '''Compute true anomaly from eccentric anomaly

        DEPRECATED

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
    
    def anomaly2nu(self, e, M):
        '''Compute true anomaly from eccentric anomaly
        Vallado
        Parameters
        ----------
        e : float
            Eccentricity
        M : float
            Mean anomaly in radians
        Returns
        -------
        nu : float
            True anomaly in radians
        E : float
            Eccentric anomaly in radians
        '''

        small = self.small
        numiter = self.iterations

        # Hyperbolic
        if (e - 1) > small:
            if e < 1.6:
                if (M < 0.0 and M > -np.pi) or M > np.pi:
                    E = M - e
                else:
                    E = M + e
            else:
                if e < 3.6 and np.abs(M) > np.pi:
                    E = M - np.sign(M) * e
                else:
                    E = M / (e - 1)
            
            cont = 0
            E1 = E + (M-e*np.sinh(E)-E) / (e*np.cosh(E)-1.0)

            while np.abs(E1-E) > small and cont < numiter:
                E = E1
                E1 = E + (M-e*np.sinh(E)-E) / (e*np.cosh(E)-1.0)
                cont += 1

            sinnu = -np.sqrt(e ** 2 - 1) * np.sinh(E) / (e - np.cosh(E))
            cosnu = (np.cosh(E) - e) / (1 - e*np.cosh(E))

            nu = np.arctan2(sinnu, cosnu)

        else:
            # Parabolic
            if np.abs(e -1) < small:
                s = 0.5 * (np.pi/2 - np.arctan(1.5*M))
                w = np.arctan(np.tan(s)**(1/3))
                E0 = 2.0* 1/np.tan(2.0*w)
                nu = 2.0*np.arctan(E0)
            else:
                # Elliptical
                if e > small:
                    if (M<0.0 and M>-np.pi) or M>np.pi:
                        E = M - e
                    else:
                        E = M + e

                    cont = 0
                    E1 = E + (M+e*np.sin(E)-E) / (1.0-e*np.cos(E))

                    while np.abs(E1-E) > small and cont < numiter:
                        E = E1
                        E1 = E + (M+e*np.sin(E)-E) / (1.0-e*np.cos(E))
                        cont += 1
                    
                    E = E1

                    sinnu = np.sqrt(1.0-e**2) * np.sin(E) / (1.0-e*np.cos(E))
                    cosnu = (np.cos(E)-e) / (1.0-e*np.cos(E))
                    nu = np.arctan2(sinnu, cosnu)

                else:
                    cont = 0
                    nu = M
                    E = M

        return nu, E

    def kepler(self, R, V, dt, savedata=False):
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

        print('r =', R_final, 'km')
        print('v =', V_final, 'km/s')

        if savedata:
            # Save data to file .csv
            data = np.array([R_final, V_final])
            np.savetxt('results/kepler.csv', data, delimiter=',')
            print('Datos guardados en results/kepler.csv')
            print('\n')

        return R_final, V_final

    def rv2coe2(self, R, V, special=None):
        '''Compute the classical orbital elements from the position and velocity
        
        DEPRECATED

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

        print('a =', a, 'km')
        print('e =', e)
        print('i =', i * 180 / np.pi, 'º')
        print('Omega =', Omega * 180 / np.pi, 'º')
        print('omega =', omega * 180 / np.pi, 'º')
        print('nu =', nu * 180 / np.pi, 'º')

        # Special cases

        w_true = np.arccos(ecc[0] / e)
        u = np.arccos(np.dot(N, R) / (n * r))
        lambda_true = np.arccos(np.dot([1, 0, 0], R) / r)


        if ecc[1] < 0: # Ecliptical orbit
            w_true = 2 * np.pi - w_true
            orb.append(['Elliptical equatorial', w_true])

        elif R[2] < 0: # Circular orbit
            u = 2 * np.pi - u
            orb.append(['Circular Inclined', u])

        elif R[1] < 0: # Equatorial orbit
            lambda_true = 2 * np.pi - lambda_true
            orb.append(['Circular equatorial', lambda_true])

        else: # Non special case
            orb.append(['Elliptical', None])

        if special == 'Circular equatorial':
            print('Circular equatorial orbit')
        elif special == 'Circular Inclined':
            print('Circular Inclined orbit')
        elif special == 'Elliptical equatorial':
            print('Elliptical equatorial orbit')

        return orb
    
    def rv2coe(self, R, V):
        '''Compute the classical orbital elements from the position and velocity
        A.8
        Parameters
        ----------
        R : numpy.array
            Position vector in km
        V : numpy.array
            Velocity vector in km/s
        Returns
        -------
        p : float
            Semi-l
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
        arglat : float
            Argument of latitude in radians
        lonper : float
            Longitude of periapsis in radians
        truelon : float
            True longitude in radians
        M : float
            Mean anomaly in radians
        '''
        

        small = self.small
        
        r = np.linalg.norm(R)
        v = np.linalg.norm(V)

        H = np.cross(R, V)
        h = np.linalg.norm(H)

        if h >= small:
            N = np.array([-H[1], H[0], 0.0])
            n = np.linalg.norm(N)

            ecc = np.dot((v**2 - self.mu / r) , R) / self.mu - np.dot(np.dot(R, V) , V) / self.mu
            e = np.linalg.norm(ecc)

            eta = v**2/2 - self.mu/r

            if abs(eta) > small:
                a = -self.mu / (2 * eta)
            else:
                a = np.inf

            p = h**2 / self.mu

            i = np.arccos(H[2] / h)

            case = 'ei'
            if e < small:
                if i < small or np.abs(np.pi - i) < small:
                    case = 'ce'
                else:
                    case = 'ci'
            else:
                if i < small or np.abs(np.pi - i) < small:
                    case = 'ee'

            # right ascension of ascending node
            if n > small:
                temp = N[0] / n
                if abs(temp) > 1:
                    temp = np.sign(temp)
                Omega = np.arccos(temp)

                if N[1] < 0:
                    Omega = 2 * np.pi - Omega

            else:
                Omega = np.nan

            # argument of perigee
            if case == 'ei':
                omega = np.arccos(np.dot(N, ecc) / (n * e)) 

                if ecc[2] < 0:
                    omega = 2 * np.pi - omega

            else :
                omega = np.nan

            # true anomaly
            if case[0] == 'e':
                nu = np.arccos(np.dot(ecc, R) / (e * r))

                if np.dot(R, V) < 0:
                    nu = 2 * np.pi - nu

            else:
                nu = np.nan

            # argument of latitude - circular inclined
            if case == 'ci' or case == 'ei':
                arglat = np.arccos(np.dot(N, R) / (n * r)) 
                
                if R[2] < 0:
                    arglat = 2 * np.pi - arglat
                m = arglat
            else:
                arglat = np.nan

            # longitude of perigee - elliptical equatorial
  
            if e > small and case == 'ee':
                temp = ecc[0] / e

                if abs(temp) > 1:
                    temp = np.sign(temp)
                lonper = np.arccos(temp)

                if ecc[1] < 0:
                    lonper = 2 * np.pi - lonper

                if i > np.pi/2:
                    lonper = 2 * np.pi - lonper

            else:
                lonper = np.nan

            # true longitude - circular equatorial

            if r > small and case == 'ce':
                temp = R[0] / r

                if np.abs(temp) > 1:
                    temp = np.sign(temp)
                truelon = np.arccos(temp)

                if R[1] < 0:
                    truelon = 2 * np.pi - truelon

                if i > np.pi/2:
                    truelon = 2 * np.pi - truelon
            else:
                truelon = np.nan

            E, M = self.nu2anomaly(e, nu)

        else:
            p = np.nan
            a = np.nan
            e = np.nan
            i = np.nan
            Omega = np.nan
            omega = np.nan
            nu = np.nan
            arglat = np.nan
            lonper = np.nan
            truelon = np.nan
            M = np.nan


        return [p, a, e, i, Omega, omega, nu, arglat, lonper, truelon, M]

    def coe2rv(self, p, e, i, Omega, omega, nu, arglat, lonper, truelon):
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
        arglat : float
            Argument of latitude in radians
        lonper : float
            Longitude of periapsis in radians
        truelon : float
            True longitude in radians
        Returns
        -------
        r : numpy.array
            Position vector in km
        v : numpy.array
            Velocity vector in km/s
        '''
        small = self.small

        if e < small:
            # circular equatorial
            omega = 0
            if i < small or np.abs(np.pi - i) < small:
                Omega = 0
                nu = truelon
            else:
                # circular inclined
                nu = arglat
        else:
            # elliptical equatorial
            if i < small or np.abs(np.pi - i) < small:

                Omega = 0
                omega = lonper

        
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

    def coe2rv2(self, p, e, i, Omega, omega, nu, extra=None):
        '''Compute the position and velocity vectors from the classical orbital elements
        
        DEPRECATED

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
            *special : str
                Special case (For parabolic and hyperbolic orbits)
                circular, equatorial, ecliptic
            *Special case (For parabolic and hyperbolic orbits)
            - u : float
                True longitude in radians
            - lambda_true : float
                True longitude in radians
            - w_true : float
                True longitude in radians
        Returns
        -------
        r : numpy.array
            Position vector in km
        v : numpy.array
            Velocity vector in km/s
        '''
        
        if extra is not None:
            if extra[0] == 'Circular equatorial':
                omega = 0
                Omega = 0
                nu = extra[1]
            elif extra[0] == 'Circular Inclined':
                omega = 0
                nu = extra[1]
            elif extra[0] == 'Elliptical equatorial':
                Omega = 0
                omega = extra[1]

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


    def Pkepler(self, ro, vo, dtsec, ndot=0, nddot=0, savedata=False):
        re = self.radius         # km
        mu = self.mu      # km3/s2
        j2 = self.body.J2      # km5/s2

        [p, a, e, i, Omega, omega, nu, arglat, truelon, lonper, m] = self.rv2coe(ro, vo)
        n = np.sqrt(mu/(a**3))

        # ------------- find the value of j2 perturbations -------------
        j2op2 = (n*1.5*re**2*j2) / (p**2)
        Omegadot = -j2op2 * np.cos(i)
        omegadot = j2op2 * (2.0-2.5*np.sin(i)**2)
        mdot = n

        a = a - 2.0*ndot*dtsec * a / (3.0*n)
        e = e - 2.0*(1.0 - e)*ndot*dtsec / (3.0*n)
        p = a*(1.0 - e**2)

        # ----- update the orbital elements for each orbit type --------
        small = 1e-10
        if e < small:
            # circular equatorial 
            if (i < small) or (abs(i-np.pi) < small):
                truelondot = Omegadot + omegadot + mdot
                truelon = truelon + truelondot * dtsec
                truelon = np.remainder(truelon, np.pi*2)
            else:
                #circular inclined  
                Omega = Omega + Omegadot * dtsec
                Omega = np.remainder(Omega, np.pi*2)
                arglatdot = omegadot + mdot
                arglat = arglat + arglatdot * dtsec
                arglat = np.remainder(arglat, np.pi*2)
        else:
            # elliptical, parabolic, hyperbolic equatorial 
            if (i < small) or (abs(i-np.pi) < small):
                lonperdot = Omegadot + omegadot
                lonper = lonper + lonperdot * dtsec
                lonper = np.remainder(lonper, np.pi*2)
                m = m + mdot*dtsec + ndot*dtsec**2 + nddot*dtsec**3
                m = np.remainder(m, np.pi*2)
                [e0, nu] = self.anomaly2nu(e, m)
            else:
                # elliptical, parabolic, hyperbolic inclined 
                Omega = Omega + Omegadot * dtsec
                Omega = np.remainder(Omega, np.pi*2)
                omega = omega + omegadot * dtsec
                omega = np.remainder(omega, np.pi*2)
                m = m + mdot*dtsec + ndot*dtsec**2 + nddot*dtsec**3
                m = np.remainder(m, np.pi*2)
                [e0, nu] = self.anomaly2nu(e, m)

        # Use coe2rv to find new

        r, v = self.coe2rv(p, e, i, Omega, omega, nu, arglat, lonper, truelon)

        if savedata:
            # Save data to file .csv
            data = np.array([r, v])
            np.savetxt('results/Pkepler.csv', data, delimiter=',')
            print('Datos guardados en results/Pkepler.csv')
            print('\n')

        return 

    def propagate(self, r0, v0, tf, dt, dn=0, ddn=0, type='Pkepler', savedata=True):
        '''Propagate the orbit forward in time

        Parameters
        ----------
        r0 : float
            Initial position vector in km
        v0 : float
            Initial velocity vector in km/s
        tf : float
            Final time in seconds
        dt : float
            Time since periapsis in seconds
        dn : float
            First derivative of mean motion vector
        ddn : float
            Second derivative of mean motion vector 
        type : str
            Type of propagator
            'Pkepler' : Perturbed Keplerian
            'kepler' : Keplerian
        Returns
        -------
        r : float
            Final position vector in km
        v : float
            Final velocity vector in km/s
        '''

        r0 = np.array(r0)
        v0 = np.array(v0)
        
        time = np.arange(0, tf, dt, dtype=float)

        # vector of position and velocity

        r = np.zeros((len(time)+1, 3))
        v = np.zeros((len(time)+1, 3))

        r[0, :] = r0
        v[0, :] = v0

        for i, dt in enumerate(time):
            if type == 'Pkepler':
                r0, v0 = self.Pkepler(r0, v0, dt, dn, ddn)
            elif type == 'kepler':
                r0, v0 = self.kepler(r0, v0, dt)
            r[i+1 , :] = r0
            v[i+1, :] = v0

        if savedata:
            # save data to file .csv with a pandas dataframe
            data = pd.DataFrame({'x': r[:, 0], 'y': r[:, 1], 'z': r[:, 2], 'vx': v[:, 0], 'vy': v[:, 1], 'vz': v[:, 2]})
            data.to_csv('results/propagate.csv', index=False)
            print('Datos guardados en results/propagate.csv')

        return r, v
    
    def RVFromTLE(self, tle_file):
        '''Calculate the position and velocity vectors from TLE data.
        Parameters
        ----------
        tle_data : dict
            Dictionary containing the TLE data.
        Returns
        -------
        R : numpy array
            Position vector in km.
        V : numpy array
            Velocity vector in km/s.
        '''

        # Decode the TLE string
        util = Utilities(self.body)

        tle_data = util.decode_tle(tle_file)

        print(f'Name: {tle_data["name"]}')

        deg = np.pi/180

        # Extract the TLE data
        inclination = tle_data['inclination'] * deg  # inclination in rad
        raan = tle_data['raan'] * deg # right ascension of the ascending node in rad
        eccentricity = tle_data['eccentricity'] # eccentricity
        arg_of_perigee = tle_data['arg_of_perigee'] * deg # argument of perigee in rad
        mean_anomaly = tle_data['mean_anomaly'] * deg  # mean anomaly in rad
        mean_motion = tle_data['mean_motion'] /24/3600 # mean motion in rad/s

        # Calculate de p and a
        a = (self.mu / mean_motion**2)**(1/3) # semi-major axis in km
        p = (1 - eccentricity**2) * a # Semi-l in km

        nu, E = self.anomaly2nu(eccentricity, mean_anomaly)

        r, v = self.coe2rv(p=p,
                            e=eccentricity,
                            i=inclination,
                            Omega=raan,
                            omega=arg_of_perigee,
                            nu=nu,
                            arglat=0,
                            lonper=0,
                            truelon=0)
        
        return r, v


if __name__ == "__main__":

    from src.CelestialBodies import CelestialBodies

    # Earth
    earth = CelestialBodies()
    earth.earth()

    kepler = KeplerPropagator(earth)

    # Test
    R = np.array([1131.340, -2282.343, 6672.423])
    V = np.array([-5.64305, 4.30333, 2.42879])
    dt = 40 * 60

    r, v = kepler.kepler(R, V, dt)
    print(r)
    print(v)

    util = Utilities(earth)

    # Decode the TLE string
    
    r0, v0 = kepler.RVFromTLE('data/paz.tle')
    print(r0)
    print(v0)

    tle_data = util.decode_tle('data/paz.tle')

    r, v = kepler.Pkepler(r0, v0, dt, .00000107, tle_data['dd_mean_motion'])
    r, v = kepler.Pkepler(R, V, dt, 0, 0)

    print(r)
    print(v)

    r, v = kepler.propagate(r0, v0, 60, .1, .00000107, 0)
    
    
    # 3d plot
    import matplotlib.pyplot as plt
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(r[:, 0], r[:, 1], r[:, 2])
    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')
    ax.set_zlabel('Z (km)')
    plt.show()
