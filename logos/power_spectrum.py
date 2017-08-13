import math
import numpy as np
from logos.grid import Grid

class PowerSpectrum:
    """
    k (np.array):
    P (np.array):
    dk (float): bin width
    """
    def __init__(self):
        pass



def compute_power_spectrum(grid, *, subtract_shotnoise=True, correct_mas=True,
                           dk=0.01, k_max=None):
    """
    Compute power spectrum of delta array
    """

    boxsize = grid.boxsize
    if k_max is None:
        k_max = grid.nc / 2

    if grid.mode == 'real-space':
        grid.fft_forward()

    grid.interlace()
    
    # window function correction
    if correct_mas:
        nc = grid.nc
        nck = nc // 2 + 1
        i = np.arange(nck)
        w = np.sin(math.pi/nc*i)
        w[0] = 1.0
        w[1:nck] /= math.pi/nc*i[1:nck]
    else:
        w = np.ones(nck)

    if subtract_shotnoise:
        shotnoise = grid.shotnoise
    else:
        shotnoise = 0.0
    
    deltak = grid.fk
    delta2 = (deltak.real**2 + deltak.imag**2)/(boxsize*w**4)

    nbin = round(k_max / dk)
    dk_fundamental = 2.0*math.pi/boxsize

    k_grid = (1 + np.arange(nck))*dk_fundamental
    
    P = np.zeros(nbin)
    k = np.zeros(nbin)

    for ik in range(nbin):
        idx = np.logical_and(dk*ik <= k_grid, k_grid < dk*(ik + 1))
        k[ik] = dk*(ik + 0.5)
        P[ik]= np.mean(delta2[idx]) - shotnoise
        #P[ik]= np.mean(delta2[(1 + ik*binwidth):(1 + (ik + 1)*binwidth)]) - shotnoise
        #k[ik] = dk*(1 + (ik + 0.5)*binwidth)

    ps = PowerSpectrum()
    ps.k = k
    ps.P = P
    ps.dk = dk
    
    return ps

