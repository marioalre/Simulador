import numpy as np
from src.Orbit import Orbit 
from src.orb2rv import Orb2rv
import matplotlib.pyplot as plt
from src.CelestialBodies import CelestialBodies


class Propagator(Orbit):
    '''Class for perturbations propagation and computation''' 

    def __init__(self, body, R0, V0, R3, V3):
        self.mu = body.mu
        self.radius = body.radius
        self.R0 = np.array(R0)
        self.V0 = np.array(V0)
        self.R3 = np.array(R3)
        self.V3 = np.array(V3)
        
    def Encke(self, t0, tf, dt, bodies, r):
        '''Encke propagation method
        Parameters
        ----------
        t0 : float
            Initial time in seconds
        tf : float
            Final time in seconds
        dt : float
            Time step in seconds
        bodies : list
            List of bodies to compute perturbations
        r : list
            List of position vectors in km of the bodies to compute perturbations
        Returns
        -------
        r : ndarray
            Position vector in km
        v : ndarray
            Velocity vector in km/s
        '''
        # Initial conditions
        r0 = np.linalg.norm(self.R0)
        v0 = np.linalg.norm(self.V0)

        dr = np.array([0, 0, 0])
        epsilon = 0
        f = 0
        Rp = self.R0
        Vp = self.V0

        t = t0 + dt

        # List of position and velocity vectors
        N = int((tf - t0)/dt)
        r_osc_l = self.R0
        v_osc_l = self.V0
        r_l = self.R0
        v_l = self.V0
        r3_l = self.R3
        v3_l = self.V3


        while t < tf:
            R_osc, V_osc = self.r0v02rv(Rp, Vp, t-t0)

            r_osc = np.linalg.norm(R_osc)
            v_osc = np.linalg.norm(V_osc)
            
            dr = R_osc - Rp # Vector perturbation of position
            dv = V_osc - Vp # Vector perturbation of velocity

            # Compute the perturbation
            epsilon = np.dot(R_osc, dr) / r_osc**2

            f = 1/epsilon * (1 - 1/(1-2*epsilon)**(3/2))
            if t == 24451200:
                print(f)

            r3, v3 = self.r0v02rv(self.R3, self.V3, t-t0)
            ad = self.acc_third_body(bodies, r)

            d2dr = ad + self.mu / r_osc**3 * (f/epsilon*Rp - dr)


            r3_l = np.concatenate((r3_l, r3))
            v3_l = np.concatenate((v3_l, v3))

            r_l = np.concatenate((r_l, R_osc))
            v_l = np.concatenate((v_l, V_osc))
  
            if (np.linalg.norm(dr)/np.linalg.norm(Rp)) > 0.01:
                R_osc = Rp
                V_osc = Vp
            else:
                R_osc = Rp + dr
                V_osc = Vp + dv

            r_osc_l = np.concatenate((r_osc_l, R_osc))
            v_osc_l = np.concatenate((v_osc_l, V_osc))
            t += dt

        r_osc_l = r_osc_l.reshape(N, 3)
        v_osc_l = v_osc_l.reshape(N, 3)
        r_l = r_l.reshape(N, 3)
        v_l = v_l.reshape(N, 3)
        r3_l = r3_l.reshape(N, 3)
        v3_l = v3_l.reshape(N, 3)

        return r_osc_l, v_osc_l, r_l, v_l, r3_l, v3_l


    def acc_third_body(self, bodies, r):
        '''Compute the acceleration due to a third body
        Parameters
        ----------
        bodies : CelestialBody
            Celestial body
        r : ndarray
            Position vector in km
        Returns
        ------- 
        acc : ndarray
            Acceleration due to third body'''

        # Mu of the bodies
        mu = np.array([body.mu for body in bodies])

        # Radius of the bodies normalized (Unpack)
        r_norm = np.array([np.linalg.norm(i) for i in r])

        # Acceleration due to third body
        a_d = np.zeros((3, 1))

        r13 = r[0] - r[2] # Radial distance between Jupiter and the Sun
        r23 = r[1] - r[2] # Radial distance between Jupiter and Mars

        r13_norm = np.linalg.norm(r13) # Normalized radius between Jupiter and the Sun
        r23_norm = np.linalg.norm(r23) # Normalized radius between Jupiter and Mars

        a_d = (((1 / mu[1]) * ((mu[1] * mu[2])/ (r23_norm ** 3))) * r23) - (((1 / mu[0]) * ((mu[0] * mu[2]) / (r13_norm ** 3))) * r13)

        return a_d
        

    def plot_orbits(self, r_osc_l, r_l, r3_l):
        '''Plot the orbits of the bodies in 3d
        Parameters
        ----------
        r_osc_l : ndarray
            Position vector in km
        r_l : ndarray
            Position vector in km
        r3_l : ndarray
            Position vector in km
        '''

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        ax.plot(r_osc_l[:, 0], r_osc_l[:, 1], r_osc_l[:, 2], label='Encke Osc')
        ax.plot(r_l[:, 0], r_l[:, 1], r_l[:, 2], label='Keplerian')
        #ax.plot(r3_l[:, 0], r3_l[:, 1], r3_l[:, 2], label='Third body')

        ax.set_xlabel('X (km)')
        ax.set_ylabel('Y (km)')
        ax.set_zlabel('Z (km)')

        # Title and legend in Spanish
        ax.set_title('Órbitas de los cuerpos')
        ax.legend()

        plt.show()


if __name__ == "__main__":
    # Bodies
    sun = CelestialBodies()
    sun.sun()
    mars = CelestialBodies()
    mars.mars()
    jupiter = CelestialBodies()
    jupiter.jupiter()

    bodies = [sun, mars, jupiter]
    # Initial conditions
    marte = Orb2rv(a = 228000000, e = 0.0934, i = 1.85, Omega = 49.562, omega = 286.537, nu = 23.33, body= mars)
    r0 = marte.r
    v0 = marte.v

    Jupiter = Orb2rv(a = 778000000, e = 0.0489, i = 1.305, Omega = 100.556, omega = 273.867, nu = 20.02, body= jupiter)
    r3 = Jupiter.r
    v3 = Jupiter.v

    # Propagation method
    encke = Propagator(sun, r0, v0, r3, v3)

    rs = np.array([0., 0., 0])           # Initial radius of the Sun
    vs = np.array([0., 0., 0])           # Initial velocity of the Sun

    r = np.array([r0, r3, rs])         # Initial radius of the bodies

    # Propagation
    r_osc_l, v_osc_l, r_l, v_l, r3_l, v3_l = encke.Encke(0, 24*3600*4000, 6*3600, bodies, r)

    # Plot orbits
    encke.plot_orbits(r_osc_l, r_l, r3_l)

    v3_0 = np.array([-13.04, -.713])  # Initial velocity of Jupiter


