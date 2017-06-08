
def velocity_field(nc, boxsize, *, sigma_u=10.0, sigmap=1.0):
    """
    Generate a velocity field with power spectrum
    Pvel(k) = A*exp[ -(k*sigmap)**2 ]
    <u(x)^2> = sigma_u

    Returns u(x) field
    """
    
    nk = nc // 2 + 1
    uk = np.zeros(nk, dtype=complex)
    dk = 2.0*math.pi/boxsize

    k = dk*np.zeros(nc//2 + 1)

    # amplitude of u(k) = sqrt(V*Pu(k))
    ampfac = amp*sigma_u/(4*math.sqrt(math.pi))
    amp = math.sqrt(ampfac*boxsize)*np.exp(-0.5*(k*sigmap)**2)
    
    uk = amp*(np.random.normal(0.0, amp, nk)
              + 1j*np.random.normal(0.0, amp, nk))

    uk[0]= 0.0;
    uk[nc // 2] = 0.0

    ux = np.fft.irfft(uk)/boxsize
    print('sigma_u', np.std(ux))

    return ux
