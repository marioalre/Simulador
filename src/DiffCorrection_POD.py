from src.utilities import Utilities
from src.kepler_propagator import KeplerPropagator
from src.Coordinates import latlongH2ecef, rotation_matrix_2, rotation_matrix_3
import numpy as np

class DiffCorrection_POD:
    '''This class is used to perform the differential correction for Preliminary Orbit Determination (POD)'''

    def __init__(self) -> None:
        pass

    def Razel(r_ecef, v_ecef, phi_gd, lambda_, h):
        '''This function calculates the range, azimuth and elevation of a satellite
        Parameters
        ----------
        r_ecef : array_like
            Position vector in ECEF frame in km
        v_ecef : array_like
            Velocity vector in ECEF frame in km/s
        phi_gd : float
            Geodetic latitude in degrees
        lambda_ : float
            Longitude in degrees
        h : float
            Height in km
        Returns
        ----------
        rho : float
            Range in km
        az : float  
            Azimuth in degrees
        el : float
            Elevation in degrees
        diff_rho : float
            Derivative of range in km/s
        diff_az : float
            Derivative of azimuth in degrees/s
        diff_el : float
            Derivative of elevation in degrees/s
        '''

        # If we have r_eci and v_eci, we have to convert them to r_ecef and v_ecef by FK5

        # height km to m
        h *= 1000
        r_site_ecef = latlongH2ecef(phi_gd, lambda_, h)

        # convert to km
        r_site_ecef /= 1000

        rho_ecef = r_ecef - r_site_ecef
        drho_ecef = v_ecef

        # convert to radians
        phi_gd_rad = np.deg2rad(phi_gd)
        lambda_rad = np.deg2rad(lambda_)

        rho_sez = np.dot(rotation_matrix_2(np.pi/2-phi_gd_rad), np.dot(rotation_matrix_3(lambda_rad), rho_ecef))
        drho_sez = np.dot(rotation_matrix_2(np.pi/2-phi_gd_rad), np.dot(rotation_matrix_3(lambda_rad), drho_ecef))

        rho = np.linalg.norm(rho_sez)

        el = np.arcsin(rho_sez[2]/rho)

        if el != np.pi/2:
            beta = np.arcsin(rho_sez[1]/(np.sqrt(rho_sez[0]**2 + rho_sez[1]**2)))
        else:
            beta = np.arcsin(drho_sez[1]/(np.sqrt(drho_sez[0]**2 + drho_sez[1]**2)))

        diff_rho = np.dot(rho_sez, drho_sez)/rho
        diff_beta = (drho_sez[0]*rho_sez[1] - drho_sez[1]*rho_sez[0])/(rho_sez[0]**2 + rho_sez[1]**2)
        diff_el = (drho_sez[2] - np.linalg.norm(drho_sez)*np.sin(el))/(np.sqrt(rho_sez[0]**2 + rho_sez[1]**2))

        return rho, beta, el, diff_rho, diff_beta, diff_el

    

        