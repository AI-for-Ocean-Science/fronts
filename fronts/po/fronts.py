""" Fronts module """
import numpy as np
from skimage.transform import resize_local_mean


from fronts.utils import stats as front_stats

try:
    from gsw import density
except ImportError:
    print("gsw not imported;  cannot do density calculations")


def anly_cutout(item:tuple, fixed_km:float=None, field_size:int=None, 
                dx:float=None):
    """Simple function to measure front related stats
    for a cutout
    
    Enables multi-processing

    Args:
        item (tuple): Items for analysis
        field_size (int, optional): Field size. Defaults to None.
        dx (float, optional): Grid spacing in km

    Returns:
        tuple: int, dict if extract_kin is False
            Otherwise, int, dict, np.ndarray, np.ndarray (F_s, gradb)
    """
    # Unpack
    Theta_cutout, Salt_cutout, idx = item
    if Theta_cutout is None or Salt_cutout is None:
        return None, idx, None

    # Calculate
    gradb = calc_gradb(Theta_cutout, Salt_cutout, dx=dx)

    # Resize
    if fixed_km is not None:
        gradb = resize_local_mean(gradb, (field_size, field_size))

    # Meta
    meta_dict = front_stats.meta_stats(gradb)

    # Return
    return gradb, idx, meta_dict


def calc_gradb(Theta:np.ndarray, Salt:np.ndarray,
             ref_rho=1025., g=0.0098, dx=2.):
    """Calculate |grad b|^2

    Args:
        Theta (np.ndarray): SST field
        Salt (np.ndarray): Salt field
        ref_rho (float, optional): Reference density
        g (float, optional): Acceleration due to gravity
            in km/s^2
        dx (float, optional): Grid spacing in km

    Returns:
        np.ndarray: |grad b|^2 field
    """
    # Buoyancy
    rho = density.rho(Salt, Theta, np.zeros_like(Salt))
    b = g*rho/ref_rho

    # Gradient
    dbdx = -1*np.gradient(b, axis=1) / dx
    dbdy = -1*np.gradient(b, axis=0) / dx

    # Magnitude
    grad_b2 = dbdx**2 + dbdy**2

    return grad_b2