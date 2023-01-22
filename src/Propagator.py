# My libraries
from src.Orbit import Orbit 
from src.orb2rv import Orb2rv
from src.perturbations import third_body, atmopheric_drag, solar_pressure
from src.CelestialBodies import CelestialBodies
from src.interplanetary import Interplanetary

# Python libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp



class Propagator(Orbit):
    '''Class for perturbations propagation and computation''' 

    def __init__(self, body):
        '''Initialize the class with the central body and the initial conditions
        Parameters
        ----------
        body : CelestialBody
            Central body
        R0 : ndarray
            Initial position vector in km
        V0 : ndarray
            Initial velocity vector in km/s
        '''

        self.body = body
        self.mu = body.mu
        self.radius = body.radius


    def perturbations(self, third_body_p=None, solar_pressure=None, drag=None):
        '''Perturbations
        Parameters
        ----------
        third_body_p : list
            - R_tb : ndarray
                Position vector of the third body
            - V_tb : ndarray
                Velocity vector of the third body
            - body_tb : CelestialBody
                Third body
        solar_pressure : ndarray
            - A : float
                Reference area in m^2
            - C_S : float
                Solar pressure coefficient
            - body : CelestialBody
                Central body
        drag : ndarray
            - A : float
                Reference area in m^2
            - C_D : float
                Drag coefficient
            - body : CelestialBody
                Central body
        Returns
        -------
        ad : ndarray
            Acceleration due to perturbations in km/s^2
            '''

        ad = 0
        # Third body perturbation
        if third_body_p is not None:
            R_tb = third_body_p[0]
            R_sc = third_body_p[1]  # spacecraft position
            body_tb = third_body_p[2]

            ad += third_body(R_tb, R_sc, body_tb)   
        else:
            ad += 0

        if solar_pressure is not None:
            ad += solar_pressure()
        else:
            ad += 0

        # Drag perturbation
        if drag is not None:
            ad += atmopheric_drag()
        else:
            ad += 0

        return ad

        
    def Encke(self, R0, V0, dt, trange, perturbations=None):
        '''Encke propagation method
        Parameters
        ----------
        R0 : ndarray
            Initial position vector in km
        V0 : ndarray
            Initial velocity vector in km/s
        dt : float
            Time step in s
        trange : list
            Time range in s [t0, tf]
        perturbations : list
            [Third body, solar pressure and drag perturbations]
            (See the perturbations function for more details)
            If there are no perturbations, set perturbations = [Third body, None, None]
            That means that the third body perturbation is computed, but the solar pressure
            and drag perturbations are not
        Returns
        -------
        rrs : list
            List of position vectors in km
        vvs : list
            List of velocity vectors in km/s
        '''

        # Compuete the perturbations
        if perturbations is not None:
            tb_params = perturbations[0]
            sp_params = perturbations[1]
            drag_params = perturbations[2]
        else:
            tb_params = None
            sp_params = None
            drag_params = None

        # Initial conditions
        t0 = trange[0]
        tf = trange[1]

        t=t0
        R = R0
        V = V0

        n = int((tf-t0)/dt)
        Rs = np.zeros((n+2, 3))
        Vs = np.zeros((n+2, 3))
        Ts = np.zeros((n+2, 1))

        Rs[0, :] = R
        Vs[0, :] = V
        Ts[0] = t

        # Third body perturbation
        if tb_params is not None:
            R3 = tb_params[0]
            V3 = tb_params[1]
            body3 = tb_params[2]

        i = 0
        while t <= tf:
            alpha = np.zeros(3)
            beta = np.zeros(3)

            orb = Orbit(self.body) 

            Rb, Vb = orb.r0v02rv(R, V, dt)

            rb = np.linalg.norm(Rb)
            r = np.linalg.norm(R)

            # Third body position and velocity
            R3, V3 = orb.r0v02rv(R3, V3, t)

            tb_params[0] = R3
            tb_params[1] = R

            if perturbations is not None:
                ad = self.perturbations(third_body_p=tb_params, solar_pressure=sp_params, drag=drag_params)
            else:
                ad = 0

            beta  = self.mu * dt * (1 - (rb/r)**3)/rb**3 + ad*dt
            alpha = beta * dt

            R = Rb + alpha
            V = Vb + beta

            Rs[i+1, :] = R
            Vs[i+1, :] = V
            Ts[i+1] = t

            i += 1
            t += dt 

        return Rs, Vs, Ts
    
    def cowell(self, t, k):
        '''Cowell method to propagate the orbit from poliastro'''

        # Initial conditions
        x, y, z = self.R0
        vx, vy, vz = self.V0

        u0 = np.array([x, y, z, vx, vy, vz])

        # Propagation with scipy solve_ivp
        sol = solve_ivp(self.func_twobody, (t[0], t[-1]), u0, args=(k,), t_eval=t)

        rrs = []
        vvs = []
        for i in range(len(t)):
            tt = t[i]
            y = sol.sol(tt)
            rrs.append(y[:3])
            vvs.append(y[3:])

        return rrs, vvs

    def func_twobody(self, t0, u_, k):
        """Differential equation for the initial value two body problem.
        From Poliastro
        Parameters
        ----------
        t0 : float
            Time.
        u_ : numpy.ndarray
            Six component state vector [x, y, z, vx, vy, vz] (km, km/s).
        k : float
            Standard gravitational parameter.
        """
        x, y, z, vx, vy, vz = u_
        r3 = (x**2 + y**2 + z**2) ** 1.5

        du = np.array([vx, vy, vz, -k * x / r3, -k * y / r3, -k * z / r3])
        return du

if __name__ == "__main__":
    # Bodies
    sun = CelestialBodies()
    sun.sun()
    tierra = CelestialBodies()
    tierra.earth()

    # spacecraft
    R0 = np.array([-0.27, 1.475, 0.001]) *1e8
    V0 = np.array([-33, -10, 1])  

    # Earth orbital elements
    a = 1.49597870700e8
    e = 0.01671123
    i = 0.0
    raan = 0.0
    argp = 0.0
    tau = -100*24*3600

    n = np.sqrt(tierra.mu/a**3)
    M = -n * tau

    # Initial conditions

    # orb = Orb2rv(a=a, e=e, Omega=raan, omega=argp, i=i, M=M, body= sun)
    # R0_tierra = orb.r
    # V0_tierra = orb.v

    R0_tierra = np.array([-2.719172990721968e+07, 1.475245243392482e+08, 0])
    V0_tierra = np.array([-29.295349119787478, -4.903140836552026, 0])

    # Inicialization
    propagator = Propagator(sun)

    # Propagation
    dt = 24*3600
    trange = [0, 100*dt]
    perturbations = [[R0_tierra, R0, tierra], None, None]

    # Encke
    Rs, Vs, Ts = propagator.Encke(R0, V0, dt, trange, perturbations)

    orb = Orbit(sun)
    r, v = orb.propagate(R0_tierra, V0_tierra, trange[1] ,dt)
    rsat, vsat = orb.propagate(R0, V0, trange[1] ,dt)

    # mirar vectores para que esten en posición heliocentrica
    ##########################################################

    # Plot

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.plot(Rs[:, 0], Rs[:, 1], Rs[:, 2], label='Encke')
    ax.plot(rsat[:, 0], rsat[:, 1], rsat[:, 2], label='Satellite')  
    ax.plot(r[:, 0], r[:, 1], r[:, 2], label='Earth')
    ax.plot([0], [0], [0], 'o', label='Sun')

    ax.set_xlabel('X [km]')
    ax.set_ylabel('Y [km]')
    ax.set_zlabel('Z [km]')


    ax.legend()
    plt.show()
    print('Done')

