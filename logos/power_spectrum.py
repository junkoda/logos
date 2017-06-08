class PowerSpectrum:
    def __init__(self, binwidth=100):
        self.binwidth = binwidth


def compute_power_spectrum(grid, shot_noise, *, binwidth=100):
    """
    Compute power spectrum of delta array
    """

    boxsize = grid.boxize

    if grid.mode == 'real-space':
        grid.fft_forward()
    
    # window function correction
    i = np.arange(nc // 2 + 1)
    w = np.sin(math.pi/nc*i)/(math.pi/nc*i)
    w[0] = 1.0

    
    deltak = grid.fk
    delta2 = (deltak.real**2 + deltak.imag**2)/boxsize

    nk = nc // (binwidth*2)
    puu = np.zeros(nk)
    dk = 2.0*math.pi/boxsize

    k = np.zeros(nk)

    for ik in range(nk-1):
        ps[ik]= np.mean(delta2[(1 + ik*kwidth):(1 + (ik + 1)*kwidth)]) - shot_noise
        k[ik] = dk*(1 + (ik + 0.5)*kwidth)
                         
    
    return k, ps
