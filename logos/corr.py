import logos._logos as c
import numpy as np

def compute_ucorr(ugrid):
    """
    Args:
       ugrid (np.array): a grid of velocity u

    Returns:
       corr (np.array): array of correlation function <u(x)u(y)>
    """

    n = ugrid.shape[0]
    ucorr = np.zeros(n)
    
    c._corr_compute_ucorr(ugrid, ucorr)

    return ucorr
