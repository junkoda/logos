import math
import numpy as np
from logos.grid import Grid

class PowerSpectrum:
    def __init__(self, binwidth=100):
        self.binwidth = binwidth


def compute_power_spectrum(grid, *, subtract_shotnoise=True, correct_mas=True,
                           binwidth=100):
    """
    Compute power spectrum of delta array
    """

    boxsize = grid.boxsize

    if grid.mode == 'real-space':
        grid.fft_forward()
    
    # window function correction
    #i = np.arange(nc // 2 + 1)
    #w = np.sin(math.pi/nc*i)/(math.pi/nc*i)
    #w[0] = 1.0

    if subtract_shotnoise:
        shotnoise = grid.shot_noise
    else:
        shotnoise = 0.0
    
    deltak = grid.fk
    delta2 = (deltak.real**2 + deltak.imag**2)/boxsize

    nbin = grid.nc // (binwidth*2) - 1
    dk = 2.0*math.pi/boxsize

    P = np.zeros(nbin)
    k = np.zeros(nbin)

    for ik in range(nbin):
        P[ik]= np.mean(delta2[(1 + ik*binwidth):(1 + (ik + 1)*binwidth)]) - shotnoise
        k[ik] = dk*(1 + (ik + 0.5)*binwidth)

    ps = PowerSpectrum()
    ps.k = k
    ps.P = P
        
    return ps
