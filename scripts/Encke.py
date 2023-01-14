from matplotlib import pyplot
import numpy
from numpy import linalg
from matplotlib import rcParams
rcParams['font.family'] = 'serif'
rcParams['font.size'] = 12

def Kepler_eqn(e, M):
    """Takes the eccentricity and mean anomaly of an orbit to solve Kepler's equation
    
    Parameters:
    ----------
        e : float
            eccentricity of orbit
        M : float
            Mean anomaly of orbit
        
    Returns:
    -------
        E : float
            Eccentric anomaly
    """
    
    E = M + e * numpy.sin(M) # eccentric anomoaly
    fofE = E - e * numpy.sin(E) - M #eccentric anomaly as a function of E
    fdotE = 1 - e * numpy.cos(E) #derivative with respect to E of fofE
    dE = - fofE / fdotE # change in E
    Enew = E + dE
    tolerance = 1e-2
    
    while abs(fofE) > tolerance:
        E = M + e * numpy.sin(Enew)
        fofE = E - numpy.sin(E) - M
        fdotE = 1 - e * numpy.cos(E)
        dE = - fofE / fdotE
        Enew = E + dE
    
    return E
    
    #Based on code from Ashish Tewari

def ellip_orb(a, Period, mu, e, t0, r0, v0, t):
    
    """Calculates the orbital position for an elliptical orbit
    
    Parameters:
    ----------
    
    a  : float
        Semi-major axis
    Period : float
        Period of planetary orbit
    mu : float
        Gravitational parameter
    t0 : float
        Initial time t = 0
    r0 : array of float
        Initial positional array
    v0 : array of float
        Initial velocity array
    t : float
        time
    
    Returns:
    -------
    
    r : array of float
        Array of radius at each time t
    v : array of float
        Array of velocity at each time t
    """
    
    r0_norm = numpy.linalg.norm(r0) # Normalized initial radius
    
    v0_norm = numpy.linalg.norm(v0) # Normalized initial velocity
    
    alpha = r0 * v0 # Constant used for Legrangian coefficients
    
    theta0 = numpy.pi # Initial true anomaly
            
    n = 2 * numpy.pi / (Period) # Mean motion, given the period

    E0 = 2 * numpy.arctan(numpy.sqrt((1 - e) / (1 + e)) * numpy.tan(0.5 * theta0)) # Initial eccentric anomaly
    
    tau = t0 + (- E0 + e * numpy.sin(E0)) / n # t - tau is time since Perigee

    M = n * (t - tau) # Mean anomaly

    E = Kepler_eqn(e, M) # Eccentric anomaly found through Kepler's equation
    
    r_leg = a * (1 - e * numpy.cos(E)) # Radius used for legrangian coefficients
    
    #Legrangian Coefficients
    
    F = 1 + a * (numpy.cos(E - E0) - 1) * r0_norm
    
    G = a * alpha * (1 - numpy.cos(E - E0)) / mu + r0_norm * numpy.sqrt(a / mu) * numpy.sin(E - E0)
    
    F_dot = - numpy.sqrt(mu * a) * (numpy.sin(E - E0)) / (r_leg * r0_norm)
    
    G_dot = 1 + a * (numpy.cos(E - E0) - 1) / r_leg
    
    r = numpy.zeros_like(r0)
    
    v = numpy.zeros_like(v0)
    
    r = F * r0 + G * v0 # Radial value of orbit for specified time
    
    v = F_dot * r0 + G_dot * v0 # Velocity value of orbit for specified time
    
    return r, v

def acceleration_d(m1, m2, m3, r, r3):
    
    """Calculates the acceleration due to the disturbing orbit
    
    Parameters:
    ----------
    m1 : float
        Mass of central body
    m2 : float
        Mass of second body
    m3 : float
        Mass of third (disturbing) body
    r : array of float
        Radial distance between body two and one
    r3: array of float
        Radial distance between body three and one
        
    Returns:
    -------
    a_d : array of float
        Acceleration due to the disturbing orbit
    """
    a_d = numpy.zeros((2, 1))
    
    G = 6.674e-11 # Gravitational constant
    
    r13 = r3 # Radial distance between Jupiter and the Sun
    
    r23 = r - r3 # Radial distance between Jupiter and Mars
    
    r23_norm = numpy.linalg.norm(r23) # Normalized radius between Jupiter and Mars
    
    r13_norm = numpy.linalg.norm(r13) # Normalized radius between Jupiter and the Sun
    
    a_d = (((1 / m2) * ((G* m2 * m3)/ (r23_norm ** 3))) * r23) - (((1 / m1) * ((G * m1 * m3) / (r13_norm ** 3))) * r13)
    
    return a_d

mu3 = 1.2669e17                  # Standard gravitational parameter of Jupiter in m^3 / s^2
m3 = 1.8983e27                   # Mass of Jupiter in kg
e3 = .0489                       # Eccentricity of Jupiter
a3 = 778000000.                  # Semi-major Axis of Jupiter in km
Period3 = 4332.589 * 3600 * 24   # Period of Jupiter Orbit in seconds

mu = 4.2828e13                   # Standard gravitational parameter of Mars in m^3 / s^2
m2 = 6.4174e23                   # Mass of Mars in kg
e = .0934                        # Eccentricity of Mars
a = 228000000.                   # Semi-major Axis of Mars in km
Period = 686.980 * 3600 * 24     # Period of Mars Orbit in seconds

mu1 = 1.3271e20                  # Standard gravitational parameters of the Sun in m^3 / s^2
m1 = 1.989e30                    # Mass of the Sun in kg

dt = 24 * 3600                   # Time step
tfinal = 4400 * dt               # Final time
N = int(tfinal / dt) + 1         # Number of time steps

t = numpy.linspace(0,tfinal,N)       # Time array
r0 = numpy.array([228000000., 0.])   # Initial radius of Mars
v0 = numpy.array([-21.84, -10.27])   # Initial velocity
r3_0 = numpy.array([778000000., 0.]) # Initial radius of Jupiter
v3_0 = numpy.array([-13.04, -.713])  # Initial velocity of Jupiter

# Set arrays for radial and velocity components

r = numpy.empty((N, 2))
v = numpy.empty((N, 2))
gamma =  numpy.empty((N, 2))
delta = numpy.empty((N, 2))
r_n = numpy.empty((N, 2))
v_n = numpy.empty((N, 2))
a_d = numpy.empty((N, 2))
r_osc = numpy.empty((N, 2))
r_osc_n = numpy.empty((N, 2))
v_osc = numpy.empty((N, 2))
v_osc_n = numpy.empty((N, 2))
r3_n = numpy.empty((N, 2))

for i,ts in enumerate(t):
    delta = numpy.zeros_like(r0)
    gamma = numpy.zeros_like(r0)
    r_osc, v_osc = ellip_orb(a, Period, mu1, e, t[0], r0, v0, ts) # Trajectory of the osculating orbit of Mars
    r_osc_norm = numpy.linalg.norm(r_osc) # Normalized osculating orbit of Mars
    r0_norm = numpy.linalg.norm(r0) # Normalized initial orbit of Mars
    r3, v3 = ellip_orb(a3, Period3, mu3, e3, t[0], r3_0, v3_0, ts) # Trajectory of Jupiter
    a_d = acceleration_d(m1, m2, m3, r_osc, r3) # Acceleration due to Jupiter
    gamma = mu3 * (dt) * ((1 - (r_osc_norm / r0_norm) ** 3) / r_osc_norm ** 3) + a_d * (dt) # Difference in velocity between osculating orbit and perturbed
    delta = gamma * (dt) # Difference between osculating orbit and perturbed orbit radius
    r = r_osc + delta # Perturbed orbital radius
    v = v_osc + gamma # Perturbed orbital velocity
    r_osc_n[i,:] = r_osc # Value of osculating orbital radius for every time step
    v_osc_n[i,:] = v_osc # Value of osculating orbital velocity for every time step
    r3_n[i,:] = r3 # Value of Jupiter's radius for every time step
    r_n[i,:] = r # Value of the perturbed orbital radius for every time step
    v_n[i,:] = v # Value of the perturbed orbital velocity for every time step

x = numpy.linspace(t[0], t[-1], N)

pyplot.figure(figsize = (8,6))
pyplot.grid(True)
pyplot.xlabel(r'X Distance (km)', fontsize = 16)
pyplot.ylabel(r'Y Distance (km)', fontsize = 16)
pyplot.title('Trajectory of Osc vs Perturbed Orbit, Flight Time = %.2f days' %(tfinal / dt), fontsize=14)
pyplot.plot(r_n[:,0], r_n[:,1])
pyplot.plot(r_osc_n[:,0], r_osc_n[:,1])
pyplot.legend(['Perturbed Orbit', 'Osculating Orbit'])
pyplot.show()

pyplot.figure(figsize = (8,6))
pyplot.grid(True)
pyplot.xlabel(r'X Distance (km)', fontsize = 16)
pyplot.ylabel(r'Y Distance (km)', fontsize = 16)
pyplot.title('Trajectory of Osc, Perturbed and Jupiter Orbit, Flight Time = %.2f days' %(tfinal / dt), fontsize=14)
pyplot.plot(r_n[:,0], r_n[:,1])
pyplot.plot(r_osc_n[:,0], r_osc_n[:,1])
pyplot.plot(r3_n[:,0], r3_n[:,1])
pyplot.legend(['Perturbed Orbit', 'Osculating Orbit', 'Jupiter\'s Orbit'])
pyplot.show()