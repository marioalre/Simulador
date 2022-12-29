# Class with the main interplanetary methods 2D

import numpy as np

class Interplanetary2D:

    def __init__(self, body, r0, r1):
        self.mu = body.mu
        self.radius = body.radius
        self.r0 = np.array(r0)
        self.r1 = np.array(r1)

    def hohmann(self):
        '''Hohmann transfer between two circular orbits
        Returns
        -------
        dv1 : float
            Delta-V for the first burn in km/s
        dv2 : float
            Delta-V for the second burn in km/s
        '''
        r0 = np.linalg.norm(self.r0)
        r1 = np.linalg.norm(self.r1)

        # Compute the semi-major axis
        a = (r0 + self.radius) / 2

        # Compute the delta-V
        dv1 = np.sqrt(self.mu / r0) * (np.sqrt(2 * self.radius / (r0 + self.radius)) - 1)
        dv2 = np.sqrt(self.mu / self.radius) * (1 - np.sqrt(2 * a / (a + self.radius)))

        time = np.sqrt(a**3 / self.mu) * np.pi

        return dv1, dv2, a, time

    def bielliptic(self, rb):
        '''Bielliptic transfer between two circular orbits
        Parameters
        ----------
        rb : float
            Radius of the intermediate orbit in km
        Returns
        -------
        dv1 : float
            Delta-V for the first burn in km/s
        dv2 : float
            Delta-V for the second burn in km/s
        dv3 : float
            Delta-V for the third burn in km/s
        '''

        r0 = np.linalg.norm(self.r0)
        r1 = np.linalg.norm(self.r1)

        # Compute the semi-major axis
        a1 = (r0 + rb) / 2
        a2 = (rb + r1) / 2

        # Compute the delta-V
        ##########################################################################################
        dv1 = np.sqrt(self.mu* (2/r0-1/a1)) -np.sqrt(self.mu / r0)
        dv2 = np.sqrt(self.mu* (2/rb-1/a2)) -np.sqrt(self.mu* (2/rb-1/a1))
        dv3 = np.sqrt(self.mu / r1) - np.sqrt(self.mu* (2/r1-1/a2))

        dV = np.abs(dv1) + np.abs(dv2) + np.abs(dv3)

        time = np.sqrt(a1**3 / self.mu) * np.pi + np.sqrt(a2**3 / self.mu) * np.pi
        ##########################################################################################
        return dv1, dv2, dv3, time

    def One_target_burn(self, nu_trans_b, place):
        '''One-target-burn transfer between two circular orbits
        Parameters
        ----------
        v_trans_b : float
            Angle of the transfer in degrees
        place : str
            Place of the burn: 'periapsis' or 'apoapsis'
        Returns
        -------
        dv1 : float
            Delta-V for the first burn in km/s
        '''

        R_inv = self.r0 / self.r1

        if place == 'periapsis':
            e_trans = (R_inv - 1) / (np.cos(nu_trans_b) - R_inv)

            a_trans = self.r0 / (1 - e_trans)

        elif place == 'apoapsis':
            e_trans = (R_inv - 1) / (np.cos(nu_trans_b) + R_inv)

            a_trans = self.r0 / (1 + e_trans)
        else:
            Warning('Place of the burn not recognized')
            print('By default, the burn is at periapsis')
            e_trans = (R_inv - 1) / (np.cos(nu_trans_b) - R_inv)

            a_trans = self.r0 / (1 - e_trans)

        v_ini = np.sqrt(self.mu / self.r0)
        v_trans_a = np.sqrt(self.mu * (2 / self.r0 - 1/a_trans))
        v_fin = np.sqrt(self.mu / self.r1)
        v_trans_b = np.sqrt(self.mu * (2 / self.r1 - 1 /a_trans))
    
        dva =  v_trans_a - v_ini

        tanphi = (e_trans * np.sin(nu_trans_b)) / (1 + e_trans * np.cos(nu_trans_b))

        phi = np.arctan(tanphi)

        dvb = np.sqrt(v_trans_b**2 + v_fin**2 - 2*v_trans_b*v_fin*np.cos(phi))

        dv_otb = np.abs(dva) + np.abs(dvb)

        cosE = (e_trans * np.sin(nu_trans_b)) / (1 + e_trans * np.cos(nu_trans_b))
        E = np.arccos(cosE)

        time = np.sqrt(a_trans**3 / self.mu) * ((E - e_trans * np.sin(E)))  ## CHECK THIS


    def inclination(di, phi_fpa, v0):
        '''Inclination change
        Parameters
        ----------
        di : float
            Inclination change in degrees
        phi_fpa : float
            Final flight path angle in degrees
            (Elliptic orbit)
        v0 : float
            Initial velocity in km/s
        Returns
        -------
        dv : float
            Delta-V in km/s
        '''

        di = np.deg2rad(di)
        phi_fpa = np.deg2rad(phi_fpa)

        return 2 * v0 * np.sin(di / 2) * np.cos(phi_fpa)

    def change_ascension_node(draan, i0, v0):
        '''Change of the ascending node for a circular orbit
        Parameters
        ----------
        draan : float
            Change of the ascending node in degrees
        v0 : float
            Initial velocity in km/s
        Returns
        -------
        dv : float
            Delta-V in km/s
        '''

        draan = np.deg2rad(draan)
        i0 = np.deg2rad(i0)

        cosnu = np.cos(i0)**2 + np.sin(i0)**2 * np.cos(draan)**2
        nu = np.arccos(cosnu)

        return 2 * v0 * np.sin(nu / 2)

    def change_ascn_inc(i_0, i_f, draan, v0):
        '''Change of the ascending node and inclination for a circular orbit
        Parameters
        ----------
        i_0 : float
            Initial inclination in degrees
        i_f : float
            Final inclination in degrees
        draan : float
            Change of the ascending node in degrees
        v0 : float
            Initial velocity in km/s
        Returns
        -------
        dv : float
            Delta-V in km/s
        '''
        
        draan = np.deg2rad(draan)
        i_0 = np.deg2rad(i_0)
        i_f = np.deg2rad(i_f)

        cosnu = np.cos(i_0)*np.cos(i_f) + np.sin(i_0)*np.sin(i_f) * np.cos(draan)
        nu = np.arccos(cosnu)

        return 2 * v0 * np.sin(nu / 2)