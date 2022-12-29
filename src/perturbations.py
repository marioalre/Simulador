# Module with the main orbital perturbations
# Path: src\perturbations.py

import numpy as np

####################################################################
# Comprobar
####################################################################

def atmopheric_drag(body, r, v, C_D, A, m):
    '''Compute the perturbation due to the atmospheric drag
    Parameters
    ----------
    body : Body
        Central body
    r : ndarray
        Position vector in km
    v : ndarray
        Velocity vector in km/s
    C_D : float
        Drag coefficient
    A : float
        Reference area in m^2
    m : float
        Mass in kg
    Returns
    -------
    dr : ndarray
        Perturbation in km
    '''
    # Compute the perturbation
    rho = body.atmospheric_density(r)
    acc = -0.5 * rho * C_D * A * np.linalg.norm(v)* v / m

    return acc

def atmospheric_density(body, r):
    '''Compute the atmospheric density for diferent altitudes
    Parameters
    ----------
    body : Body
        Central body
    r : ndarray
        Position vector in km
    Returns
    -------
    rho : float
        Atmospheric density in kg/m^3
    '''
    # Compute the atmospheric density for diferent altitudes ISA model
    h = np.linalg.norm(r) - body.radius
    
    if h < 11000:
        rho = 1.225 * (1 - 0.000006875 * h)**4.256
    elif h < 20000:
        rho = 0.36391 * np.exp(-0.0001577 * h)
    elif h < 32000:
        rho = 0.08803 * np.exp(-0.0001186 * h)
    elif h < 47000:
        rho = 0.01322 * np.exp(-0.0000457 * h)
    elif h < 51000:
        rho = 0.00143 * np.exp(-0.000011 * h)
    elif h < 71000:
        rho = 0.00086 * np.exp(-0.0000026 * h)
    else:
        rho = 0

    return rho

def third_body(body, r, m):
    '''Compute the perturbation due to a third body
    Parameters
    ----------
    body : Body
        Third body
    r : ndarray
        Position vector in km
    m : float
        Mass in kg
    Returns
    -------
    dr : ndarray
        Perturbation in km
    '''
    # Compute the perturbation
    acc = - body.mu * r / np.linalg.norm(r)**3

    return acc

def solar_pressure(body, r, C_S, A, m):
    '''Compute the perturbation due to the solar pressure
    Parameters
    ----------
    body : Body
        Central body
    r : ndarray
        Position vector in km
    C_S : float
        Solar pressure coefficient
    A : float
        Reference area in m^2
    m : float
        Mass in kg
    Returns
    -------
    dr : ndarray
        Perturbation in km
    '''
    # Compute the perturbation
    acc = - body.S * C_S * A * r / m

    return acc

def gravity(body, r):
    '''Compute the perturbation due to the gravity field
    Parameters
    ----------
    body : Body
        Central body
    r : ndarray
        Position vector in km
    Returns
    -------
    dr : ndarray
        Perturbation in km
    '''
    # Compute the perturbation
    acc = - body.mu * r / np.linalg.norm(r)**3

    return acc

def j2(body, r):
    '''Compute the perturbation due to the J2 effect
    Parameters
    ----------
    body : Body
        Central body
    r : ndarray
        Position vector in km
    Returns
    -------
    dr : ndarray
        Perturbation in km
    '''
    # Compute the perturbation
    r_norm = np.linalg.norm(r)
    r_xy = np.sqrt(r[0]**2 + r[1]**2)
    acc = - 3/2 * body.mu * body.j2 * body.radius**2 * r / r_norm**5 * (5 * r[2]**2 / r_norm**2 - 1)
    acc[2] = acc[2] - 3/2 * body.mu * body.j2 * body.radius**2 * r[2] / r_norm**5 * (5 * r[2]**2 / r_norm**2 - 3)

    return acc

def j3(body, r):
    '''Compute the perturbation due to the J3 effect
    Parameters
    ----------
    body : Body
        Central body
    r : ndarray
        Position vector in km
    Returns
    -------
    dr : ndarray
        Perturbation in km
    '''
    # Compute the perturbation
    r_norm = np.linalg.norm(r)
    r_xy = np.sqrt(r[0]**2 + r[1]**2)
    acc = 15/8 * body.mu * body.j3 * body.radius**3 * r / r_norm**7 * (7 * r[2]**2 / r_norm**2 - 3)
    acc[2] = acc[2] - 15/8 * body.mu * body.j3 * body.radius**3 * r[2] / r_norm**7 * (7 * r[2]**2 / r_norm**2 - 4)

    return acc

def j4(body, r):
    '''Compute the perturbation due to the J4 effect
    Parameters
    ----------
    body : Body
        Central body
    r : ndarray
        Position vector in km
    Returns
    -------
    dr : ndarray
        Perturbation in km
    '''
    # Compute the perturbation
    r_norm = np.linalg.norm(r)
    r_xy = np.sqrt(r[0]**2 + r[1]**2)
    acc = - 35/8 * body.mu * body.j4 * body.radius**4 * r / r_norm**9 * (3 - 42 * r[2]**2 / r_norm**2 + 63 * r[2]**4 / r_norm**4)
    acc[2] = acc[2] - 35/8 * body.mu * body.j4 * body.radius**4 * r[2] / r_norm**9 * (3 - 42 * r[2]**2 / r_norm**2 + 70 * r[2]**4 / r_norm**4)

    return acc

def j5(body, r):
    '''Compute the perturbation due to the J5 effect
    Parameters
    ----------
    body : Body
        Central body
    r : ndarray
        Position vector in km
    Returns
    -------
    dr : ndarray
        Perturbation in km
    '''
    # Compute the perturbation
    r_norm = np.linalg.norm(r)
    r_xy = np.sqrt(r[0]**2 + r[1]**2)
    acc = 63/8 * body.mu * body.j5 * body.radius**5 * r / r_norm**11 * (35 * r[2]**2 / r_norm**2 - 30 * r[2]**4 / r_norm**4 + 7 * r[2]**6 / r_norm**6)
    acc[2] = acc[2] - 63/8 * body.mu * body.j5 * body.radius**5 * r[2] / r_norm**11 * (35 * r[2]**2 / r_norm**2 - 60 * r[2]**4 / r_norm**4 + 35 * r[2]**6 / r_norm**6)

    return acc

def j6(body, r):
    '''Compute the perturbation due to the J6 effect
    Parameters
    ----------
    body : Body
        Central body
    r : ndarray
        Position vector in km
    Returns
    -------
    dr : ndarray
        Perturbation in km
    '''
    # Compute the perturbation
    r_norm = np.linalg.norm(r)
    r_xy = np.sqrt(r[0]**2 + r[1]**2)
    acc = - 231/16 * body.mu * body.j6 * body.radius**6 * r / r_norm**13 * (63 - 462 * r[2]**2 / r_norm**2 + 715 * r[2]**4 / r_norm**4 - 429 * r[2]**6 / r_norm**6 + 77 * r[2]**8 / r_norm**8)
    acc[2] = acc[2] - 231/16 * body.mu * body.j6 * body.radius**6 * r[2] / r_norm**13 * (63 - 462 * r[2]**2 / r_norm**2 + 1001 * r[2]**4 / r_norm**4 - 1001 * r[2]**6 / r_norm**6 + 429 * r[2]**8 / r_norm**8)

    return acc
