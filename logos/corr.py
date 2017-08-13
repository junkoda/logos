import logos._logos as c
import numpy as np

def compute_ucorr(ugrid, nr):
    """
    Args:
       ugrid (np.array): an array of velocity u
       nr (int): length of output correlation funtion

    Returns:
       corr (np.array): array of correlation function <u(x)u(y)>
    """

    ucorr = np.zeros(nr)
    
    c._corr_compute_ucorr(ugrid, ucorr)

    return ucorr

def compute_umoments(ugrid, nr):
    """
    Args:
       ugrid (np.array): an array of velocity u
       nr (int): length of output correlation funtion

    Returns:
       corr (np.array): array of correlation function <u(x)u(y)>
    """

    umoments = np.zeros((nr, 3))
    
    c._corr_compute_umoments(ugrid, umoments)

    return umoments
