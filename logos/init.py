import math
import numpy as np

def velocity_field(Pvel, nc, boxsize, *, sigma_u=10.0, sigmap=1.0):
    """Generate a velocity field with power spectrum

    Args:
        Pvel np.array -> np.array: velocity power spectrum (function)
        sigma_u (float): normalize Pvel(k) to this velocity rms

    Pvel(k) = A*exp[ -(k*sigmap)**2 ]
    <u(x)^2> = sigma_u

    Returns u(x) array, A
    """
    
    nk = nc // 2 + 1
    uk = np.zeros(nk, dtype=complex)
    dk = 2.0*math.pi/boxsize

    k = dk*np.arange(nk)

    pk = Pvel(k) #np.exp(-(k*sigmap)**2)
    
    if nc % 2 == 0:
        pk[nc // 2] = 0.0

    pksum = np.sum(pk)
    sigma_u_raw = math.sqrt(2.0/boxsize*pksum)

    mag = (sigma_u/sigma_u_raw)*np.sqrt(0.5*boxsize*pk)
    
    uk = mag*(np.random.normal(0.0, 1.0, nk)
                    + 1j*np.random.normal(0.0, 1.0, nk))    
    uk[0]= 0.0;

    ux = np.fft.irfft(uk)*(nc/boxsize)

    return ux, (sigma_u/sigma_u_raw)**2
