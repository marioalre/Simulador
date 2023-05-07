import numpy as np
import datetime
from src.Orbit import Orbit
from src.rv2orb import Rv2orb

class Utilities:
    def __init__(self, body=None):
        if body is not None:
            self.mu = body.mu
            self.radius = body.radius
            self.central_body = body
        else:
            self.mu = 398600.4418 # km^3/s^2
            self.radius = 6378.1363 # km
            print('Warning: No body selected. Using Earth as default.')

    def findTOF(self, R0, R1, p):
        '''
        Find time of flight
        Parameters
        ----------
        r0 : numpy array
            Initial position vector
        r1 : numpy array
            Final position vector
        p : float
            Semi-latus
        Returns
        -------
        TOF : float
            Time of flight


        ATTENTION: This function is not working properly
        '''

        r0 = np.linalg.norm(R0)
        r1 = np.linalg.norm(R1)

        cosdv = np.dot(R0, R1) / (r0 * r1)

        k = r0 * r1 * (1 - cosdv)

        l = r0 + r1

        m = r0 * r1 * (1 + cosdv)

        a = (m * k * p) / ((2*m-l**2)*p**2 + 2*k*l*p - k**2)

        f = 1- r1 / p * (1 - cosdv)

        g = r0 * r1 * np.sin(np.arccos(cosdv)) / np.sqrt(p * self.mu)

        if a > 0:
            # Elliptic
            df = np.sqrt(self.mu / p) * np.tan(np.arccos(cosdv)/2) *((1 - cosdv) / p - 1 / r0 - 1 / r1) 

            cosde = 1 - r0 / a * (1 - f)

            sinde = - r0 * r1 * df / (np.sqrt(a * self.mu))

            TOF = g + np.sqrt(a**3 / self.mu) * (np.arccos(cosde) - sinde)

        elif a == np.inf:
            # Parabolic
            c = np.sqrt(r0**2 + r1**2 - 2 * r0 * r1 * cosdv)

            s = (r0 + r1 + c) / 2

            TOF = 2/3 * np.sqrt(s**3 / (2*self.mu)) * (1 - ((s - c) / s)**(3/2))

        else:
            # Hyperbolic
            cosdh = 1 + r0 / a * (f - 1)

            TOF = g + np.sqrt((-a)**3 / self.mu) * (np.sinh(np.arccosh(cosdh)) - np.arccosh(cosdh))

        return TOF


    def Gibbs(self, R0, R1, R2):
        '''
        Gibs method to find v2. the position vectors are given in the geometric frame.
        This method performs well when the angles between the vectors are larger than 5º.
        Parameters
        ----------
        r0 : numpy array
            Initial position vector
        r1 : numpy array
            Intermediate position vector
        r2 : numpy array
            Final position vector
        Returns
        ----------
        v1 : numpy array
            Intermediate velocity vector'''

        # Parameters
        irr = 0
        tol = 1e-10

        R0 = np.array(R0)
        R1 = np.array(R1)
        R2 = np.array(R2)
        
        # Normalize position vectors
        r0 = np.linalg.norm(R0)
        r1 = np.linalg.norm(R1)
        r2 = np.linalg.norm(R2)

        # Angles between position vectors
        alpha01 = np.arccos(np.dot(R0, R1) / r0 / r1)
        alpha12 = np.arccos(np.dot(R1, R2) / r1 / r2)

        # Show yhe angles in degrees
        print('The angles between the vectors are: ', np.rad2deg(alpha01), np.rad2deg(alpha12))

        # Cross product
        r0xr1 = np.cross(R0, R1)
        r1xr2 = np.cross(R1, R2)
        r2xr0 = np.cross(R2, R0)

        # Coplanarity
        if np.abs(np.dot(R0, r1xr2) / r0 / np.linalg.norm(r1xr2)) < 1e-10:
            print('Coplanarity condition not satisfied')
            irr = 1

        N = r0 * r1xr2 + r1 * r2xr0 + r2 * r0xr1

        D = r0xr1 + r1xr2 + r2xr0

        S = R0 * (r1 - r2) + R1 * (r2 - r0) + R2 * (r0 - r1)

        return np.sqrt(self.mu / np.linalg.norm(N) / np.linalg.norm(D)) * (np.cross(D, R1) / r1 + S), irr

    def HERRICK_GIBBS(self, R0, R1, R2, t0, t1, t2):
        '''This method is used to find the velocity vectors at the intermediate position like
        the Gibbs method, but it uses the time of flight to find the velocity vectors. It performs
        better than the Gibbs method when the position vectors are close to each other (<1º).
        
        Parameters
        ----------
        R0 : numpy array
            Initial position vector
        R1 : numpy array
            Intermediate position vector
        R2 : numpy array
            Final position vector
        t0 : float
            Initial time
        t1 : float
            Intermediate time
        t2 : float
            Final time
        Returns
        ----------
        V1 : numpy array
            Intermediate velocity vector
        '''

        # Parameters
        dt01 = np.abs(t1 - t0)
        dt12 = np.abs(t2 - t1)
        dt02 = np.abs(t2 - t0)

        r0 = np.linalg.norm(R0)
        r1 = np.linalg.norm(R1)
        r2 = np.linalg.norm(R2)

        Z12 = np.cross(R1, R2)

        alpha_cop = np.pi / 2 - np.arccos(np.dot(Z12, R0) / np.linalg.norm(Z12) / np.linalg.norm(R0))

        # Angles between the vectors
        alpha01 = np.arccos(np.dot(R0, R1) / r0 / r1)
        alpha12 = np.arccos(np.dot(R1, R2) / r1 / r2)

        # Show yhe angles in degrees
        print('The angles between the vectors are: ', np.rad2deg(alpha01), np.rad2deg(alpha12))

        A = np.dot(-dt12 * (1/(dt01 * dt02) + self.mu/(12*r0**3)), R0)
        B = np.dot((dt12 - dt01) * (1/(dt01 * dt12) + self.mu/(12*r2**3)), R1)
        C = np.dot(dt01 * (1/(dt12 * dt02) + self.mu/(12*r2**3)), R2)

        v2 = A + B  + C 

        return v2


##############################################################################################################
#  OBSERVATION FUNCTIONS to calculate the position and velocity vectors from the observation of a satellite  #
##############################################################################################################

    def rv_from_observation(self, rho, drho, A, dA, a, da, theta, phi, H):
        '''This function calculates the position and velocity vectors from the observation of a satellite.
        Parameters
        ----------
        rho : float
            Range in km
        drho : float
            Range rate in km/s
        A : float
            Azimuth in degrees
        dA : float
            Azimuth rate in degrees/s
        a : float
            Elevation in degrees
        da : float
            Elevation rate in degrees/s
        theta : float
            Local sideral time in degrees
        phi : float
            Latitude in degrees
        H : float
            Height in km
        Returns
        ----------
        R : numpy array
            Position vector in km
        V : numpy array
            Velocity vector in km/s
        '''

        deg2rad = np.pi / 180
        w = 7.292115e-5
        f = 1 / 298.257223563 # WGS84 flattening
        omega = np.array([0, 0, w])

        # Convert to radians
        A = A * deg2rad
        dA = dA * deg2rad
        a = a * deg2rad
        da = da * deg2rad
        phi = phi * deg2rad
        theta = theta * deg2rad

        # Calculate the position vector
        R = (self.radius/np.sqrt(1 - (2 * f - f**2)*np.sin(phi)**2) + H) * np.array([np.cos(phi) * np.cos(theta), np.cos(phi) * np.sin(theta), 0]) + (self.radius*(1-f)**2/np.sqrt(1 - (2 * f - f**2)*np.sin(phi)**2) + H)* np.array([0, 0, np.sin(phi)])

        # Calculate the topocentric declination
        delta = np.arcsin(np.sin(phi) * np.sin(a) + np.cos(phi) * np.cos(a) * np.cos(A))

        # topocentric right ascension
        if 0 < A < 2*np.pi:
            h = 2*np.pi -  np.arccos((np.sin(a) * np.cos(phi) - np.sin(phi) * np.cos(a) * np.cos(A)) / (np.cos(delta)))
        elif 180 <= A <= 360 or A == 0:
            h = np.arccos((np.sin(a) * np.cos(phi) - np.sin(phi) * np.cos(a) * np.cos(A)) / (np.cos(delta)))
        else:
            Warning('The azimuth angle is not in the range [0, 360)')

        alpha = theta - h

        # direction cosine unit vector
        u = np.array([np.cos(delta) * np.cos(alpha), np.cos(delta) * np.sin(alpha), np.sin(delta)])

        # goecentric position vector
        r = R + rho * u

        # inertial velocity vector
        dR = np.cross(omega, R)

        # declination rate
        ddelta = (-np.cos(a) * np.sin(A) * np.cos(phi) * dA - np.cos(phi) * np.sin(a) * np.cos(A) * da + np.sin(phi) * np.cos(a) * da) / (np.cos(delta))

        # Right ascension rate
        dalpha = w + (np.cos(a) * np.cos(A) * dA - np.sin(a) * np.sin(A) * da + np.sin(A) * np.cos(a)* np.tan(delta) * ddelta) / (np.cos(phi) * np.sin(a) - np.sin(phi) * np.cos(a) * np.cos(A))

        # direction cosine rate vector
        du = np.array([-np.sin(alpha) * np.cos(delta) * dalpha - np.cos(alpha) * np.sin(delta) * ddelta, -np.sin(delta) * np.sin(alpha) * ddelta + np.cos(delta) * np.cos(alpha) * dalpha, np.cos(delta) * ddelta])

        # geocentric velocity vector
        v = dR + drho * u + rho * du

        return r, v

###############################################################################################################
# Gauss' method with iterative improvement
###############################################################################################################
    def Gauss_POD(self, q1, q2, q3, R1, R2, R3, t1, t2, t3):
        '''This function calculates the position and velocity vectors from the observation of a satellite.
        Gauss' method is used to calculate the position and velocity vectors.
        Parameters
        ----------
        q1 : numpy array
            Position vector in km (satellite from earth surface)
        q2 : numpy array
            Position vector in km (satellite from earth surface)
        q3 : numpy array
            Position vector in km (satellite from earth surface)
        R1 : numpy array
            Position vector in km (observer)
        R2 : numpy array
            Position vector in km (observer)
        R3 : numpy array
            Position vector in km (observer)
        t1 : float
            Time in seconds first observation
        t2 : float
            Time in seconds second observation
        t3 : float
            Time in seconds third observation
        Returns
        ----------
        R : numpy array
            Position vector in km
        V : numpy array
            Velocity vector in km/s
        '''
        self.q1 = q1
        self.q2 = q2
        self.q3 = q3
        self.R1 = R1
        self.R2 = R2
        self.R3 = R3
        self.t1 = t1
        self.t2 = t2
        self.t3 = t3

        # time differences
        self.tau1 = t1 - t2
        self.tau3 = t3 - t2
        self.tau = self.tau3 - self.tau1

        # cross products
        p1 = np.cross(q2, q3)
        p2 = np.cross(q1, q3)
        p3 = np.cross(q1, q2)

        D0 = np.dot(q1, p1)
        self.D0 = D0

        D11 = np.dot(R1, p1)
        D21 = np.dot(R2, p1)
        D31 = np.dot(R3, p1)
        self.D11 = D11
        self.D21 = D21
        self.D31 = D31

        D12 = np.dot(R1, p2)
        D22 = np.dot(R2, p2)
        D32 = np.dot(R3, p2)
        self.D12 = D12
        self.D22 = D22
        self.D32 = D32

        D13 = np.dot(R1, p3)
        D23 = np.dot(R2, p3)
        D33 = np.dot(R3, p3)
        self.D13 = D13
        self.D23 = D23
        self.D33 = D33

        A = 1/D0 * (-D12*(self.tau3/self.tau) + D22 + D32 * (self.tau1/self.tau))
        B = 1/(D0 * 6) * (D12 * (self.tau3**2 - self.tau**2)*(self.tau3/self.tau) + D32 * (self.tau**2 - self.tau1**2)*(self.tau1/self.tau))

        E = np.dot(R2, q2)
        r22 = np.dot(R2, R2)

        a = -(A**2 + 2* A * E + r22)
        b = -2 * B * (A + E) * self.mu
        c = -B**2 * self.mu**2

        roots = np.roots([1, 0, a, 0, 0, b, 0, 0, c])

        for root in roots:
            if np.isreal(root) and root > 0:
                r = root.real
                print(f'Radius: {r} km')
                if r < self.radius:
                    Warning('The radius is smaller than the radius of the earth')
                break
        
        Q2 = A + (self.mu * B) / r**3
        self.Q2 = Q2

        num = 6 *(D31 * self.tau1/self.tau3 + D21 * self.tau/self.tau3) * r**3 + self.mu * D31 * (self.tau**2 - self.tau1**2) * self.tau1/self.tau3
        den = 6*r**3 + self.mu *(self.tau**2 - self.tau3**2)
        Q1 = 1/D0 * ((num)/(den) - D11)
        self.Q1 = Q1

        num1 = 6 *(D13 * self.tau3/self.tau1 - D23 * self.tau/self.tau1) * r**3 + self.mu * D13 * (self.tau**2 - self.tau3**2) * self.tau3/self.tau1
        Q3 = 1/D0 * ((num1)/(den) - D33)
        self.Q3 = Q3

        # position vectors
        r1 = R1 + Q1 * q1
        r2 = R2 + Q2 * q2
        r3 = R3 + Q3 * q3

        # velocity vectors
        self.f1 = 1 - 0.5*self.mu * (self.tau1**2) / r**3
        self.g1 = self.tau1 - (1/6) * self.mu * self.tau1**3 / r**3

        self.f3 = 1 - 0.5*self.mu * (self.tau3**2) / r**3
        self.g3 = self.tau3 - (1/6) * self.mu * self.tau3**3 / r**3

        v2 = (-self.f3 * r1 + self.f1 * r3) / (self.g3 * self.f1 - self.g1 * self.f3)

        print('Now you can use the function "rv2coe" to calculate the orbital elements.')

        self.r2 = r2
        self.v2 = v2

        return r2, v2

    def rv2coe(self):
        '''This function calculates the orbital elements from the position and velocity vectors.'''
        orb = Rv2orb(self.r2, self.v2, self.central_body, 0)
        deg = 180 / np.pi
        # Print the orbital elements
        print(f'Orbital elements of the satellite:')
        print(f'Inclination: {orb.i * deg} deg')
        print(f'Right ascension of the ascending node: {orb.Omega* deg} deg')
        print(f'Eccentricity: {orb.ecc}')
        print(f'Argument of perigee: {orb.omega * deg} deg')
        print(f'Semi-major axis: {orb.a} km')
        print(f'true anomaly: {orb.nu * deg} deg')

        return orb

    def iterative(self):
        '''Using the Gauss method to calculate the position and velocity vectors by iterating over the observations.
        '''

        r, v = self.Gauss_POD(self.q1, self.q2, self.q3, self.R1, self.R2, self.R3, self.t1, self.t2, self.t3)

        diff1 = 1
        diff2 = 1
        diff3 = 1

        tol = 1e-9

        Q1_old = self.Q1
        Q2_old = self.Q2
        Q3_old = self.Q3 

        nmax = 1000
        n = 0

        while (diff1 > tol and diff2 > tol and diff3 > tol) and n < nmax:
            
            n += 1

            r_n = np.linalg.norm(r)
            v_n = np.linalg.norm(v)

            alpha = 2/r_n - v_n**2/self.mu

            vr2 = np.dot(v, r) / r_n

            prop = Orbit(self.central_body)

            chi1 = prop.Kepler_universal(r_n, vr2, self.tau1, alpha)
            chi3 = prop.Kepler_universal(r_n, vr2, self.tau3, alpha)

            # Lagrange coefficients
            ff1, gg1 = prop.lagrange_coeff(chi1, self.tau1, r_n, alpha)
            ff3, gg3 = prop.lagrange_coeff(chi3, self.tau3, r_n, alpha)

            # Update f and g
            self.f1 = (ff1 + self.f1) / 2
            self.g1 = (gg1 + self.g1) / 2
            self.f3 = (ff3 + self.f3) / 2
            self.g3 = (gg3 + self.g3) / 2

            c1 = self.g3 / (self.f1 * self.g3 - self.f3 * self.g1)
            c3 = -self.g1 / (self.f1 * self.g3 - self.f3 * self.g1)

            self.Q1 = 1/self.D0 * (-self.D11 + 1/c1 * self.D21 - c3/c1 * self.D31)
            self.Q2 = 1/self.D0 * (-c1 * self.D12 + self.D22 - c3 * self.D32)
            self.Q3 = 1/self.D0 * (-c1/c3 * self.D13 + 1/c3*self.D23 - self.D33)

            self.r1 = self.R1 + self.Q1 * self.q1
            self.r2 = self.R2 + self.Q2 * self.q2
            self.r3 = self.R3 + self.Q3 * self.q3

            self.v2 = (-self.f3 * self.r1 + self.f1 * self.r3) / (self.g3 * self.f1 - self.g1 * self.f3)

            # difference between old and new values
            diff1 = np.linalg.norm(self.Q1 - Q1_old)
            diff2 = np.linalg.norm(self.Q2 - Q2_old)
            diff3 = np.linalg.norm(self.Q3 - Q3_old)

            # Update old values
            Q1_old = self.Q1
            Q2_old = self.Q2
            Q3_old = self.Q3

            # print iiterations and values 
            print('Iteration: ', n)
            print(f'Q1: {self.Q1}   Q2: {self.Q2}   Q3: {self.Q3}')
            print(f'f1: {self.f1}   g1: {self.g1}   f3: {self.f3}   g3: {self.g3}')
            print(f'chi1: {chi1}   chi3: {chi3}')


        # Display de number of iterations
        print('The number of iterations is: ', n)

        return self.r2, self.v2

    def r_with_ralst(self, phi, theta, f, H):
        '''This function calculates the position vector of the satellite in the ECI frame.
        
        Parameters
        ----------
        phi : float
            Right ascension in degrees.
        theta : float
            Local sidereal time of the satellite in degrees.
        f : float
            Flattening of the Earth.
        H : float
            Height of the satellite above the Earth's surface in km.
        '''
        # Convert to radians
        phi = phi * np.pi / 180
        theta = theta * np.pi / 180
        
        R = (self.radius/np.sqrt(1 - (2 * f - f**2)*np.sin(phi)**2) + H) * np.array([np.cos(phi) * np.cos(theta), np.cos(phi) * np.sin(theta), 0]) + (self.radius*(1-f)**2/np.sqrt(1 - (2 * f - f**2)*np.sin(phi)**2) + H)* np.array([0, 0, np.sin(phi)])

        return R

    def q_with_rade(self, alpha, delta):
        '''This function calculates the position vector of the satellite in the ECI frame.
        
        Parameters
        ----------
        alpha : float
            Latitude of the satellite in degrees.
        delta : float
            declination of the satellite in degrees.
        '''
        # Convert to radians
        alpha = alpha * np.pi / 180
        delta = delta * np.pi / 180

        Q = np.array([np.cos(alpha) * np.cos(delta), np.sin(alpha) * np.cos(delta), np.sin(delta)])
        
        return Q

##############################################################################################################
#  TLE files                                                                                                 #
##############################################################################################################
    
    def decode_tle(self, tle_path):
        '''Decode a two-line element set (TLE) into orbital elements.'''

        # Read the TLE file
        with open(tle_path, 'r') as tle_file:
            tle = tle_file.read()


        # Split the TLE into lines
        lines = tle.strip().split('\n')

        # Extract the satellite name and TLE data from the lines
        name = lines[0].strip()
        line1 = lines[1].strip()
        line2 = lines[2].strip()

        # Extract the orbital elements from the TLE data
        inclination = float(line2[8:16])
        raan = float(line2[17:25])
        eccentricity = float("0." + line2[26:33])
        arg_of_perigee = float(line2[34:42])
        mean_anomaly = float(line2[43:51])
        mean_motion = float(line2[52:63])
        epoch_year = int(line1[18:20])
        epoch_day = float(line1[20:32])

        # derivate of mean motion
        d_mean_motion = float(line1[34:43])
        # Second derivate of mean motion
        dd_mean_motion = float(line1[45:50]) * 10**(float(line1[50:52]) + 5.)

        # Calculate the epoch time in UTC
        epoch_year += 1900 if epoch_year >= 57 else 2000
        epoch_day = int(epoch_day)
        epoch_time = datetime.datetime(year=epoch_year, month=1, day=1) + datetime.timedelta(days=epoch_day - 1)

        # Return the extracted values as a dictionary
        tle_data = {
            'name': name,
            'inclination': inclination,
            'raan': raan,
            'eccentricity': eccentricity,
            'arg_of_perigee': arg_of_perigee,
            'mean_anomaly': mean_anomaly,
            'mean_motion': mean_motion,
            'epoch': epoch_time,
            'd_mean_motion': d_mean_motion,
            'dd_mean_motion': dd_mean_motion
        }
        return tle_data
    

if __name__ == '__main__':

    from CelestialBodies import CelestialBodies

    Tierra = CelestialBodies()
    Tierra.earth()

    util = Utilities(Tierra)

    alpha = np.array([43.537, 54.420, 64.318])
    delta = np.array([-8.7833, -12.074, -15.105])
    lst = np.array([44.506, 45.000,  45.499])

    q = np.zeros((3, 3))
    r = np.zeros((3, 3))

    for i in range(len(alpha)):
        q[i] = util.q_with_rade(alpha[i], delta[i])
        r[i] = util.r_with_ralst(40, lst[i], Tierra.f, 1)
    print(f'q{i}: {q} \nr{i}: {r}')

    t1 = 0
    t2 = 118.10
    t3 = 237.58

    r, v = util.Gauss_POD(q[0], q[1], q[2], r[0], r[1], r[2], t1, t2, t3)
    print('The position vector is: ', r)
    print('The velocity vector is: ', v)


    r, v = util.iterative()
    print('The position vector is: ', r)
    print('The velocity vector is: ', v)

    orb = util.rv2coe()

    data = util.decode_tle('data/paz.tle')
    print(f'Name: {data["name"]}')
    print(f'Inclination: {data["inclination"]}')
    print(f'RAAN: {data["raan"]}')
    print(f'Eccentricity: {data["eccentricity"]}')
    print(f'Argument of perigee: {data["arg_of_perigee"]}')
    print(f'Mean anomaly: {data["mean_anomaly"]}')
    print(f'Mean motion: {data["mean_motion"]}')
    print(f'Epoch: {data["epoch"]}')
    print(f'd_mean_motion: {data["d_mean_motion"]}')
    print(f'dd_mean_motion: {data["dd_mean_motion"]}')
   

        

