import math
import numpy as np
import logos.cosmology

def compute_projected_velocity_power(power, *, a=None):
    """Compute 1D power spectrum P_1D(k) from 3D power spectrum P(k)
    Args:
        power (2D array): k P

    Returns:
        P1 (array): 1D velocity power spectrum
        P1(k) = \int dky kz (kx/ky)
    """

    k = power[:, 0]
    P = power[:, 1]
    
    n = power.shape[0]
    assert(len(P) == n)
    P1 = np.zeros_like(P)

    k2 = k**2
    dk2 = k2[1:] - k2[:-1]
    Puu = P/k2**2
    f = 0.5*(Puu[1:] + Puu[:-1])
    
    for i in range(n):
        # trapezoidal integral
        P1[i] = math.pi*k2[i]*np.sum(f[i:]*dk2[i:])
    
    return P1

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
