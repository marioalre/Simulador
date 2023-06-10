# This class allows to convert ECI coordinates to ECEF coordinates following the algorithm described in the
# Fundamentals of Astrodynamics and Applications book by David A. Vallado.

# Import the necessary libraries
import numpy as np
import math

class ECI2ECEF:
    ''' This class allows to convert ECI coordinates to ECEF coordinates following the algorithm described in the
        Fundamentals of Astrodynamics and Applications book by David A. Vallado.
    '''

    # Define the constructor that takes the input parameters
    def __init__(self, reci, veci, aeci, ttt, jdut1, lod, xp, yp, eqeterms, ddpsi, ddeps):
        ''' 
        This function initializes the class with the input parameters.
        Inputs:
            reci: ECI position vector (km)
            veci: ECI velocity vector (km/s)
            ttt: Julian centuries of terrestrial time
            jdut1: Julian date of UT1
            lod: Length of day (sec)
            xp: Polar motion coefficient (arcsec)
            yp: Polar motion coefficient (arcsec)
            eqeterms: Boolean to select if the extra terms for nutation are used (0 or 2)
            ddpsi: Correction to delta psi (arcsec)
            ddeps: Correction to delta eps (arcsec)
        '''
        self.reci = reci
        self.veci = veci
        self.aeci = aeci
        self.ttt = ttt
        self.jdut1 = jdut1
        self.lod = lod
        self.xp = xp
        self.yp = yp
        self.eqeterms = eqeterms
        self.ddpsi = ddpsi
        self.ddeps = ddeps

    # Define the function that converts ECI coordinates to ECEF
    def eci2ecef(self):
        # Call the auxiliary functions to get the precession, nutation and rotation parameters
        prec, psia, wa, ea, xa = self.precess(self.ttt, '80')
        deltapsi, trueeps, meaneps, omega, nut = self.nutation(self.ttt, self.ddpsi, self.ddeps)
        st, stdot = self.sidereal(self.jdut1, deltapsi, meaneps, omega, self.lod, self.eqeterms)
        pm = self.polarm(self.xp, self.yp, self.ttt, '80')

        # Calculate the angular velocity of the Earth
        thetasa = 7.29211514670698e-05 * (1.0 - self.lod / 86400.0)
        omegaearth = np.array([0, 0, thetasa])

        # Convert ECI coordinates to ECEF using transformation matrices
        rpef = st.T @ nut.T @ prec.T @ self.reci
        recef = pm.T @ rpef

        vpef = st.T @ nut.T @ prec.T @ self.veci - np.cross(omegaearth, rpef)
        vecef = pm.T @ vpef

        temp = np.cross(omegaearth, rpef)

        # Two additional terms that are not necessary if the satellite is not on the surface of the Earth
        aecef = pm.T @ (st.T @ nut.T @ prec.T @ self.aeci) - np.cross(omegaearth,temp) - 2.0 * np.cross(omegaearth,vpef)

        # Return the results
        return recef, vecef, aecef

    # Define the auxiliary functions for precession, nutation and rotation

    def get_ttt(self, JD_tt):
        '''This function calculates the Julian centuries of terrestrial time
        Parameters
        ----------
        JD_tt : float
            Julian date of TT
        Returns
        ----------
        ttt : float
            Julian centuries of terrestrial time
        '''

        ttt = (JD_tt - 2451545.0) / 36525.0


    def precess(self, ttt, opt):
        '''Based on Vallado's book, Fundamentals of Astrodynamics and Applications, 4th ed., Algorithm 23, p. 211
        and the Matlab code aviable at celestrak.com

        VALIDATED

        Parameters
        ----------
        ttt : float
            Julian centuries of terrestrial time
        opt : str
            '50' for fk4 b1950 precession angles
            '80' for iau 76 precession angles
            '00' for iau 06 precession angles
        Returns
        ----------
        prec : array
            Precession matrix
        '''
        # Code for precess function goes here
        convrt = np.pi / (180.0*3600.0)
        ttt2 = ttt * ttt
        ttt3 = ttt2 * ttt

        prec = np.zeros((3, 3))

        # fk4 b1950 precession angles
        if opt == '50':
            t1 = 0.0  # (ttt - 2433282.42345905) / 365242.198782
            t2 = (ttt - 2433282.42345905) / 36525
            psia = 50.3708 + 0.0050 * ttt
            wa = 0.0
            ea = 84428.26 - 46.845 * ttt - 0.00059 * ttt2 + 0.00181 * ttt3
            xa = 0.1247 - 0.0188 * ttt
            zeta = (23035.545 + 139.720 * t1 + 0.060 * t1 * t1) * ttt + (30.240 - 0.270 * t1) * ttt2 + 17.995 * ttt3
            theta = (20051.12 - 85.29 * t1 - 0.37 * t1 * t1) * ttt + (-42.65 - 0.37 * t1) * ttt2 - 41.80 * ttt3
            z = (23035.545 + 139.720 * t1 + 0.060 * t1 * t1) * ttt + (109.480 + 0.390 * t1) * ttt2 + 18.325 * ttt3

            prec[0, 0] = 1.0 - 2.9696e-4 * ttt2 - 1.3e-7 * ttt3
            prec[0, 1] = 2.234941e-2 * ttt + 6.76e-6 * ttt2 - 2.21e-6 * ttt3
            prec[0, 2] = 9.7169e-3 * ttt - 2.07e-6 * ttt2 - 9.6e-7 * ttt3
            prec[1, 0] = -prec[0, 1]
            prec[1, 1] = 1.0 - 2.4975e-4 * ttt2 - 1.5e-7 * ttt3
            prec[1, 2] = -1.0858e-4 * ttt2
            prec[2, 0] = -prec[0, 2]
            prec[2, 1] = prec[1, 2]
            prec[2, 2] = 1.0 - 4.721e-5 * ttt2

            # pass these back out for testing
            psia = zeta
            wa = theta
            ea = z

            return prec

        # iau 76 precession angles
        elif opt == '80':

            psia = 5038.7784 * ttt - 1.07259 * ttt2 - 0.001147 * ttt3
            wa = 84381.448 + 0.05127 * ttt2 - 0.007726 * ttt3
            ea = 84381.448 - 46.8150 * ttt - 0.00059 * ttt2 + 0.001813 * ttt3
            xa = 10.5526 * ttt - 2.38064 * ttt2 - 0.001125 * ttt3

            zeta = 2306.2181 * ttt + 0.30188 * ttt2 + 0.017998 * ttt3
            theta = 2004.3109 * ttt - 0.42665 * ttt2 - 0.041833 * ttt3
            z = 2306.2181 * ttt + 1.09468 * ttt2 + 0.018203 * ttt3


        # iau 06 precession angles
        else:
            oblo = 84381.406
            psia = (((( -0.0000000951 * ttt + 0.000132851 ) * ttt - 0.00114045 ) * ttt - 1.0790069 ) * ttt + 5038.481507 ) * ttt
            wa = ((((  0.0000003337 * ttt - 0.000000467 ) * ttt - 0.00772503 ) * ttt + 0.0512623 ) * ttt - 0.025754 ) * ttt + oblo
            ea = (((( -0.0000000434 * ttt - 0.000000576 ) * ttt + 0.00200340 ) * ttt - 0.0001831 ) * ttt - 46.836769 ) * ttt + oblo
            xa = (((( -0.0000000560 * ttt + 0.000170663 ) * ttt - 0.00121197 ) * ttt - 2.3814292 ) * ttt + 10.556403 ) * ttt

            zeta = (((( -0.0000003173 * ttt - 0.000005971 ) * ttt + 0.01801828 ) * ttt + 0.2988499 ) * ttt + 2306.083227 ) * ttt + 2.650545
            theta = (((( -0.0000001274 * ttt - 0.000007089 ) * ttt - 0.04182264 ) * ttt - 0.4294934 ) * ttt + 2004.191903 ) * ttt
            z = ((((  0.0000002904 * ttt - 0.000028596 ) * ttt + 0.01826837 ) * ttt + 1.0927348 ) * ttt + 2306.077181 ) * ttt - 2.650545


        # convert from arcseconds to radians
        psia *= convrt
        wa *= convrt
        ea *= convrt
        xa *= convrt
        zeta *= convrt
        theta *= convrt
        z *= convrt

        # calculate the precession matrix elements
        cz = np.cos(z)
        sz = np.sin(z)
        coszeta = math.cos(zeta)
        sinzeta = math.sin(zeta)
        costheta = math.cos(theta)
        sintheta = math.sin(theta)
        cosz = math.cos(z)
        sinz = math.sin(z)

        # Formar la matriz mod to j2000
        prec = np.zeros((3, 3))
        prec[0, 0] = coszeta * costheta * cosz - sinzeta * sinz
        prec[0, 1] = coszeta * costheta * sinz + sinzeta * cosz
        prec[0, 2] = coszeta * sintheta
        prec[1, 0] = -sinzeta * costheta * cosz - coszeta * sinz
        prec[1, 1] = -sinzeta * costheta * sinz + coszeta * cosz
        prec[1, 2] = -sinzeta * sintheta
        prec[2, 0] = -sintheta * cosz
        prec[2, 1] = -sintheta * sinz
        prec[2, 2] = costheta


        return prec, psia, wa, ea, xa


    def nutation(self, ttt, ddpsi, ddeps):
        '''Based on Vallado's book, Fundamentals of Astrodynamics and Applications, 4th ed., Algorithm 22, p. 205
        and the Matlab code aviable at celestrak.com
        '''

        deg2rad = math.pi / 180.0

        iar80, rar80 = self.iau80in()  # coeff in deg

        # ---- determine coefficients for iau 1980 nutation theory ----
        ttt2 = ttt * ttt
        ttt3 = ttt2 * ttt

        meaneps = -46.8150 * ttt - 0.00059 * ttt2 + 0.001813 * ttt3 + 84381.448
        meaneps = (meaneps / 3600.0) % 360.0
        meaneps = meaneps * deg2rad

        l, l1, f, d, omega, lonmer, lonven, lonear, lonmar, lonjup, lonsat, lonurn, lonnep, precrate = self.fundarg(ttt, '80')

        deltapsi = 0.0
        deltaeps = 0.0
        for i in range(105, -1, -1):
            tempval = iar80[i - 1, 0] * l + iar80[i - 1, 1] * l1 + iar80[i - 1, 2] * f + iar80[i - 1, 3] * d + iar80[i - 1, 4] * omega
            deltapsi = deltapsi + (rar80[i - 1, 0] + rar80[i - 1, 1] * ttt) * math.sin(tempval)
            deltaeps = deltaeps + (rar80[i - 1, 2] + rar80[i - 1, 3] * ttt) * math.cos(tempval)

        # --------------- find nutation parameters --------------------
        a1 = deltapsi + ddpsi
        a2 = deltaeps + ddeps

        if a1 < 0.0:
            deltapsi = np.remainder(-a1,  2.0 * math.pi)
            deltapsi = -deltapsi
        else:
            deltapsi = np.remainder(a1,  2.0 * math.pi)

        if a2 < 0.0:
            deltaeps = np.remainder(-a2, 2.0 * math.pi)
            deltaeps = -deltaeps
        else:
            deltaeps = np.remainder(a2, 2.0 * math.pi)

        trueeps = meaneps + deltaeps

        cospsi = math.cos(deltapsi)
        sinpsi = math.sin(deltapsi)
        coseps = math.cos(meaneps)
        sineps = math.sin(meaneps)
        costrueeps = math.cos(trueeps)
        sintrueeps = math.sin(trueeps)

        nut = np.zeros((3, 3))
        nut[0, 0] = cospsi
        nut[0, 1] = costrueeps * sinpsi
        nut[0, 2] = sintrueeps * sinpsi
        nut[1, 0] = -coseps * sinpsi
        nut[1, 1] = costrueeps * coseps * cospsi + sintrueeps * sineps
        nut[1, 2] = sintrueeps * coseps * cospsi - sineps * costrueeps
        nut[2, 0] = -sineps * sinpsi
        nut[2, 1] = costrueeps * sineps * cospsi - sintrueeps * coseps
        nut[2, 2] = sintrueeps * sineps * cospsi + costrueeps * coseps

        return deltapsi, trueeps, meaneps, omega, nut
    
    def iau80in(self):
        # 0.0001" to rad
        convrt = 0.0001 * np.pi / (180 * 3600.0)

        nut80 = np.loadtxt("./data/nut80.dat")

        iar80 = nut80[:, 0:5]
        rar80 = nut80[:, 5:9]

        for i in range(106):
            for j in range(4):
                rar80[i, j] = rar80[i, j] * convrt

        return iar80, rar80
    

    def fundarg(self, ttt, opt):
        deg2rad = np.pi / 180.0

        # ---- determine coefficients for iau 2000 nutation theory ----
        ttt2 = ttt * ttt
        ttt3 = ttt2 * ttt
        ttt4 = ttt2 * ttt2

        # ---- iau 2006 theory
        if opt == '06':
            l = 134.96340251 + (1717915923.2178 * ttt + 31.8792 * ttt2 + 0.051635 * ttt3 - 0.00024470 * ttt4) / 3600.0
            l1 = 357.52910918 + (129596581.0481 * ttt - 0.5532 * ttt2 - 0.000136 * ttt3 - 0.00001149 * ttt4) / 3600.0
            f = 93.27209062 + (1739527262.8478 * ttt - 12.7512 * ttt2 + 0.001037 * ttt3 + 0.00000417 * ttt4) / 3600.0
            d = 297.85019547 + (1602961601.2090 * ttt - 6.3706 * ttt2 + 0.006593 * ttt3 - 0.00003169 * ttt4) / 3600.0
            omega = 125.04455501 + (-6962890.5431 * ttt + 7.4722 * ttt2 + 0.007702 * ttt3 - 0.00005939 * ttt4) / 3600.0

            lonmer = 252.250905494 + 149472.6746358 * ttt
            lonven = 181.979800853 + 58517.8156748 * ttt
            lonear = 100.466448494 + 35999.3728521 * ttt
            lonmar = 355.433274605 + 19140.299314 * ttt
            lonjup = 34.351483900 + 3034.90567464 * ttt
            lonsat = 50.0774713998 + 1222.11379404 * ttt
            lonurn = 314.055005137 + 428.466998313 * ttt
            lonnep = 304.348665499 + 218.486200208 * ttt
            precrate = 1.39697137214 * ttt + 0.0003086 * ttt2

        # ---- iau 2000b theory
        if opt == '02':
            l = 134.96340251 + (1717915923.2178 * ttt) / 3600.0
            l1 = 357.52910918 + (129596581.0481 * ttt) / 3600.0
            f = 93.27209062 + (1739527262.8478 * ttt) / 3600.0
            d = 297.85019547 + (1602961601.2090 * ttt) / 3600.0
            omega = 125.04455501 + (-6962890.2665 * ttt + 7.4722 * ttt ** 2 + 0.007702 * ttt ** 3 - 0.00005939 * ttt ** 4) / 3600.0

            lonmer = 0.0
            lonven = 0.0
            lonear = 0.0
            lonmar = 0.0
            lonjup = 0.0
            lonsat = 0.0
            lonurn = 0.0
            lonnep = 0.0
            precrate = 0.0

        # iau 1996 theory
        if opt == '96':
            l = 134.96340251 + (1717915923.2178 * ttt + 31.8792 * ttt ** 2 + 0.051635 * ttt ** 3 - 0.00024470 * ttt ** 4) / 3600.0
            l1 = 357.52910918 + (129596581.0481 * ttt - 0.5532 * ttt ** 2 - 0.000136 * ttt ** 3 - 0.00001149 * ttt ** 4) / 3600.0
            f = 93.27209062 + (1739527262.8478 * ttt - 12.7512 * ttt ** 2 + 0.001037 * ttt ** 3 + 0.00000417 * ttt ** 4) / 3600.0
            d = 297.85019547 + (1602961601.2090 * ttt - 6.3706 * ttt ** 2 + 0.006593 * ttt ** 3 - 0.00003169 * ttt ** 4) / 3600.0
            omega = 125.04455501 + (-6962890.2665 * ttt + 7.4722 * ttt ** 2 + 0.007702 * ttt ** 3 - 0.00005939 * ttt ** 4) / 3600.0
            lonven = 181.979800853 + 58517.8156748 * ttt
            lonear = 100.466448494 + 35999.3728521 * ttt
            lonmar = 355.433274605 + 19140.299314 * ttt
            lonjup = 34.351483900 + 3034.90567464 * ttt
            lonsat = 50.0774713998 + 1222.11379404 * ttt
            precrate = 1.39697137214 * ttt + 0.0003086 * ttt ** 2

        # iau 1980 theory
        if opt == '80':
            l = ((((0.064) * ttt + 31.310) * ttt + 1717915922.6330) * ttt) / 3600.0 + 134.96298139
            l1 = ((((-0.012) * ttt - 0.577) * ttt + 129596581.2240) * ttt) / 3600.0 + 357.52772333
            f = ((((0.011) * ttt - 13.257) * ttt + 1739527263.1370) * ttt) / 3600.0 + 93.27191028
            d = ((((0.019) * ttt - 6.891) * ttt + 1602961601.3280) * ttt) / 3600.0 + 297.85036306
            omega = ((((0.008) * ttt + 7.455) * ttt - 6962890.5390) * ttt) / 3600.0 + 125.04452222
            lonmer = 252.3 + 149472.0 * ttt
            lonven = 179.9 + 58517.8 * ttt
            lonear = 98.4 + 35999.4 * ttt
            lonmar = 353.3 + 19140.3 * ttt
            lonjup = 32.3 + 3034.9 * ttt
            lonsat = 48.0 + 1222.1 * ttt
            lonurn  =   0.0
            lonnep  =   0.0
            precrate=   0.0

        l = np.remainder(l, 360.0) * deg2rad
        l1 = np.remainder(l1, 360.0) * deg2rad
        f = np.remainder(f, 360.0) * deg2rad
        d = np.remainder(d, 360.0) * deg2rad
        omega = np.remainder(omega, 360.0) * deg2rad
        lonmer = np.remainder(lonmer, 360.0) * deg2rad
        lonven = np.remainder(lonven, 360.0) * deg2rad
        lonear = np.remainder(lonear, 360.0) * deg2rad
        lonmar = np.remainder(lonmar, 360.0) * deg2rad
        lonjup = np.remainder(lonjup, 360.0) * deg2rad
        lonsat = np.remainder(lonsat, 360.0) * deg2rad
        lonurn = np.remainder(lonurn, 360.0) * deg2rad
        lonnep = np.remainder(lonnep, 360.0) * deg2rad
        precrate = np.remainder(precrate, 360.0) * deg2rad

        return l, l1, f, d, omega, lonmer, lonven, lonear, lonmar, lonjup, lonsat, lonurn, lonnep, precrate


    def sidereal(self, jdut1, deltapsi, meaneps, omega, lod, eqeterms):
        gmst = self.gstime(jdut1)

        if (jdut1 > 2450449.5) and (eqeterms > 0):
            ast = gmst + deltapsi * math.cos(meaneps) + 0.00264 * math.pi / (3600 * 180) * math.sin(omega) + 0.000063 * math.pi / (3600 * 180) * math.sin(2.0 * omega)
        else:
            ast = gmst + deltapsi * math.cos(meaneps)

        ast = ast % (2.0 * math.pi)
        thetasa = 7.29211514670698e-05 * (1.0 - lod / 86400.0)
        omegaearth = thetasa

        st = np.array([[math.cos(ast), -math.sin(ast), 0.0],
            [math.sin(ast), math.cos(ast), 0.0],
            [0.0, 0.0, 1.0]])

        stdot = np.array([[-omegaearth * math.sin(ast), -omegaearth * math.cos(ast), 0.0],
                [omegaearth * math.cos(ast), -omegaearth * math.sin(ast), 0.0],
                [0.0, 0.0, 0.0]])

        return st, stdot
    

    def gstime(self, jdut1):
        twopi = 2.0 * math.pi
        deg2rad = math.pi / 180.0

        tut1 = (jdut1 - 2451545.0) / 36525.0

        gst = -6.2e-6 * tut1 * tut1 * tut1 + 0.093104 * tut1 * tut1 + (876600.0 * 3600.0 + 8640184.812866) * tut1 + 67310.54841

        gst = (gst * deg2rad / 240.0) % (twopi)

        if gst < 0.0:
            gst += twopi

        return gst
    

    def polarm(self, xp, yp, ttt, opt):
        cosxp = math.cos(xp)
        sinxp = math.sin(xp)
        cosyp = math.cos(yp)
        sinyp = math.sin(yp)

        pm = np.zeros((3, 3))

        if opt == '80':
            pm[0, 0] = cosxp
            pm[0, 1] = 0.0
            pm[0, 2] = -sinxp
            pm[1, 0] = sinxp * sinyp
            pm[1, 1] = cosyp
            pm[1, 2] = cosxp * sinyp
            pm[2, 0] = sinxp * cosyp
            pm[2, 1] = -sinyp
            pm[2, 2] = cosxp * cosyp
        else:
            convrt = math.pi / (3600.0 * 180.0)
            sp = -47.0e-6 * ttt * convrt
            cossp = math.cos(sp)
            sinsp = math.sin(sp)

            pm[0, 0] = cosxp * cossp
            pm[0, 1] = -cosyp * sinsp + sinyp * sinxp * cossp
            pm[0, 2] = -sinyp * sinsp - cosyp * sinxp * cossp
            pm[1, 0] = cosxp * sinsp
            pm[1, 1] = cosyp * cossp + sinyp * sinxp * sinsp
            pm[1, 2] = sinyp * cossp - cosyp * sinxp * sinsp
            pm[2, 0] = sinxp
            pm[2, 1] = -sinyp * cosxp
            pm[2, 2] = cosyp * cosxp

        return pm



if __name__ == '__main__':

    # seguntos to arco segundos
    convrt = np.pi / (180.0 * 3600.0)

    r_eci = [5.1025, 6.1230, 6.3781]
    v_eci = [-4.7432, 0.7905, 5.5338]
    a_eci = [0.0026, 0.0001, 0.0030]

    ttt = 0.0426236319
    jdut1 = 2.4531e6
    lod = 0.001556300
    xp = -6.820455828585174e-07
    yp = 1.615927632369383e-06
    eqeterms = 2
    ddpsi = -2.530485008551223e-7
    ddeps = -1.878653014299452e-8
    conv = ECI2ECEF(r_eci, v_eci, a_eci, ttt, jdut1, lod, xp, yp, eqeterms, ddpsi, ddeps)

    recef, vecef, aecef = conv.eci2ecef()
    print(recef)  #coregir



