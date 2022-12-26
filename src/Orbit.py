import numpy as np
from matplotlib import pyplot as plt

class Orbit():
    def __init__(self, body):
        self.mu = body.mu
        self.radius = body.radius

    def kepler_elliptic(self, M, e):
        '''Kepler equation for elliptic orbits by Newton-Raphson method
        Parameters
        ----------
        M : float 
            Mean anomaly rad 
        e : float
            Eccentricity 
        Returns
        -------
        E : float
            Eccentric anomaly rad
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
            Time since periapsis in seconds
        p : float

        e : float
            Eccentricity
        Returns
        -------
        B : float
            Parabolic anomaly in radians
        '''

        n_p = 2* np.sqrt(self.mu / p**3)

        cot2s = 3/2 * dt * n_p

        s = 0.5*np.arctan(1/cot2s)

        w = np.arctan(np.cbrt(np.tan(s)))

        return 2 * 1 / np.tan(2 * w)

    def kepler_hyperbolic(self, M, e):
        '''Kepler equation for hyperbolic orbits by Newton-Raphson method

        Parameters
        ----------
        M : float
            Mean anomaly in radians
        e : float
            Eccentricity
        Returns 
        -------
        F : float
            Hyperbolic anomaly in radians
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

        while abs(F - e * np.sinh(F) + M) > 1e-10:
            F = F + (F - e * np.sinh(F) + M) / (e * np.cosh(F) - 1)

        return F

    def stumpff_s(self, z):
        '''Stumpff function S(z)

        Parameters
        ----------
        z : float
            z = x / a
        Returns
        -------
        S : float
            Stumpff function
        '''
        if z > 0:
            S = (np.sqrt(z) - np.sin(np.sqrt(z))) / np.sqrt(z)**3
        elif z < 0:
            S = (np.sinh(np.sqrt(-z)) - np.sqrt(-z)) / np.sqrt(-z)**3
        else:
            S = 1 / 6

        return S

    def stumpff_c(self, z):
        '''Stumpff function C(z)

        Parameters
        ----------
        z : float
            z = x / a
        Returns
        -------
        C : float
            Stumpff function
        '''
        if z > 0:
            C = (1 - np.cos(np.sqrt(z))) / z
        elif z < 0:
            C = (np.cosh(np.sqrt(-z)) - 1) / (-z)
        else:
            C = 1 / 2

        return C

    def Kepler_universal(self, r0, v0 , dt, alpha):
        '''Universal Kepler equation for elliptic, parabolic and hyperbolic orbits

        Parameters
        ----------
        r0 : float
            Initial position vector in km
        v0 : float
            Initial velocity vector in km/s
        dt : float
            Time since periapsis in seconds
        alpha : float
            alpha = 1/a  (a = semi-major axis) 1/km
        Returns
        -------
        chi : float
            Universal anomaly in Km^0.5
        '''

        # Initial guess

        chi = np.sqrt(self.mu) * np.abs(alpha) * dt

        # Newton-Raphson method

        ratio =  1
        max_iter = 1000
        iter = 0

        while np.abs(ratio) > 1e-8 and max_iter > iter:

            C = self.stumpff_c(chi**2 * alpha)
            S = self.stumpff_s(chi**2 * alpha)

            P1 = (r0 * v0) / np.sqrt(self.mu) * chi**2 * C
            P2 = (1 - r0 * alpha) * chi**3 * S + r0 * chi - np.sqrt(self.mu) * dt
            f_chi =  P1 + P2

            Q1 =(r0 * v0) / np.sqrt(self.mu) * chi * (1 - chi**2 * alpha * S) 
            Q2 = (1 - r0 * alpha) * chi**2 * C + r0
            f_prime_chi = Q1 + Q2 

            ratio =  f_chi / f_prime_chi

            chi = chi - ratio

            iter += 1
        return chi

    def lagrange_coeff(self, chi, t, r0, alpha):
        '''Lagrange coefficients for elliptic, parabolic and hyperbolic orbits

        Parameters
        ----------
        chi : float
            Universal anomaly in Km^0.5
        t : float
            Time since periapsis in seconds
        ro : float
            Initial position vector in km
        alpha : float
            alpha = 1/a  (a = semi-major axis) 1/km
        Returns
        -------
        f : float
            Lagrange coefficient
        g : float
            Lagrange coefficient

        '''

        C = self.stumpff_c(chi**2 * alpha)
        S = self.stumpff_s(chi**2 * alpha)

        f = 1 - chi**2 / r0 * C
        g = t - chi**3 * S / np.sqrt(self.mu)

        return f, g

    def d_lagrange_coeff(self, chi, r, r0, alpha):
        '''Derivative of Lagrange coefficients for elliptic, parabolic and hyperbolic orbits

        Parameters
        ----------
        chi : float
            Universal anomaly in Km^0.5
        r : float
            Radial position in km
        ro : float
            Initial position vector in km
        alpha : float
            alpha = 1/a  (a = semi-major axis) 1/km
        Returns
        -------
        f : float
            Lagrange coefficient
        g : float
            Lagrange coefficient

        '''

        C = self.stumpff_c(chi**2 * alpha)
        S = self.stumpff_s(chi**2 * alpha)

        df = np.sqrt(self.mu) / r0 / r *(chi**2 * alpha * S - 1) * chi
        dg = 1 - chi**2 / r * C

        return df, dg

    def r0v02rv(self, r0, v0, dt):
        '''Final position and velocity for elliptic, parabolic and hyperbolic orbits 
        from initial position, velocity and time since periapsis
        
        Parameters
        ----------
        r0 : float
            Initial position vector in km
        v0 : float
            Initial velocity vector in km/s
        dt : float
            Elapsed time in senconds
        Returns
        -------
        r : float
            Final position vector in km
        v : float
            Final velocity vector in km/s
        '''

        # Normalization
        r0n = np.linalg.norm(np.array(r0))
        v0n = np.linalg.norm(np.array(v0))

        # Radial component of initial velocity
        vr0 = np.dot(r0, v0) / r0n

        # alpha
        alpha = (2 / r0n - v0n**2 / self.mu)

        # The sign of alpha determines the type of orbit
        # alpha > 0 : elliptic orbit
        # alpha = 0 : parabolic orbit
        # alpha < 0 : hyperbolic orbit

        chi = self.Kepler_universal(r0n, vr0, dt, alpha)

        f, g = self.lagrange_coeff(chi, dt, r0n, alpha)

        r = f * r0 + g * v0

        rn = np.linalg.norm(np.array(r))

        df, dg = self.d_lagrange_coeff(chi, rn, r0n, alpha)

        v = df * r0 + dg * v0

        return r, v

    def propagate(self, r0, v0, tf, dt):
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

        for i, t in enumerate(time):
            r0, v0 = self.r0v02rv(r0, v0, dt)
            r[i+1 , :] = r0
            v[i+1, :] = v0

        plt.style.use('dark_background')
        # plt.rcParams['grid.color'] = "black"

        fig = plt.figure()

        ax = fig.add_subplot(111, projection='3d')
        # ax.w_xaxis.set_pane_color((0.0, 0.0, 0.0, 0.8))
        # ax.w_yaxis.set_pane_color((0.0, 0.0, 0.0, 0.6))
        # ax.w_zaxis.set_pane_color((0.0, 0.0, 0.0, 1.0))

        ax.plot(r[:, 0], r[:, 1], r[:, 2], label='orbit')

        ax.plot([0], [0], [0], 'o', label='Earth', color='green', markersize=self.radius/1000)


        ax.set_xlabel('X [km]')

        ax.set_ylabel('Y [km]')

        ax.set_zlabel('Z [km]')

        ax.legend()

        plt.show()


        return r, v, ax




    



        

    

    