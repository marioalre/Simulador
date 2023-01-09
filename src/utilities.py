import numpy as np

class Utilities:
    def __init__(self, body=None):
        if body is not None:
            self.mu = body.mu
            self.radius = body.radius
        else:
            self.mu = 398600.4418 # km^3/s^2

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

    def J0(self, year, month, day):
        '''This function calculates the Julian date of a given date at 0h UT. 
        Between 1900 and 2100.
        Parameters
        ----------
        year : int
            Year
        month : int
            Month
        day : int
            Day
        Returns
        ----------
        J0 : float
            Julian date
        '''

        return 367 * year - np.floor(7 * (year + np.floor((month + 9) / 12)) / 4) + np.floor(275 * month / 9) + day + 1721013.5

    
    def localSideralTime(self, year, month, day, ut, EL):
        '''This function calculates the local sideral time at a given date and time.
        Parameters
        ----------
        year : int
            Year
        month : int
            Month
        day : int
            Day
        ut : float
            Universal time in hours
        EL : float
            East longitude in degrees
        Returns
        ----------
        theta : float
            Local sideral time in degrees
        '''
        
        j0 = self.J0(year, month, day)
        j = (j0 - 2451545) / 36525

        g0 = 100.4606184 + 36000.77004 * j + 0.000387933 * j**2 - 2.583e-8 * j**3

        theta = g0 + 360.98564724 * ut + EL

        # Convert to the range [0, 360)
        g0 = self.zero2360(g0)

        lst = g0 + 360.98564724 * ut/24 + EL

        # Convert to the range [0, 360)
        lst = self.zero2360(lst)
        
        return lst

    def zero2360(self, x):
        '''Converts an angle in degrees to the range [0, 360)'''

        if x >= 360:
            x = x - 360 * np.floor(x / 360)
        elif x < 0:
            x = x + 360 * np.ceil(np.abs(x) / 360)
        
        return x

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
