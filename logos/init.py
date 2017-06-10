import math
import numpy as np
import logos.cosmology
from scipy.interpolate import interp1d

def compute_projected_velocity_power(power, *, a=None):
    """Compute 1D power spectrum P_1D(k) from 3D power spectrum P(k)
    Args:
        power (2D array): k P
        a: normalise the power spectrum by (f*D(a))^2

    Returns:
        P1 (array): 1D velocity power spectrum
        P1(k) = \int dky kz/(2pi)^2 (kx/k)^2 P(k)/k^2
        where k = sqrt(kx^2 + ky^2 + kz^2)
    """

    fac = 1.0/(4.0*math.pi)
    if a is not None:
        # Normalise by (D(a)*f)^2
        D = logos.cosmology.D_growth_factor(a)
        f = logos.cosmology.f_growth_rate(a)
        fac *= (f*D)**2

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
        P1[i] = k2[i]*np.sum(f[i:]*dk2[i:])
    
    return fac*P1

def velocity_field(Pvel, nc, boxsize, *, sigma_u=None, kmin=None):
    """Generate a velocity field with power spectrum

    Args:
        Pvel: velocity power spectrum (function k -> Pvel)
              or np.array of k Pvel
        sigma_u (float): normalize Pvel(k) to this velocity rms
                         <u(x)^2> = sigma_u

    Returns u(x) array, A
    """
    
    nk = nc // 2 + 1
    uk = np.zeros(nk, dtype=complex)
    dk = 2.0*math.pi/boxsize

    k = dk*np.arange(nk)
    pk = np.zeros(nk)

    if callable(Pvel):
        f = Pvel
    else:
        # If Pvel is a 2D array of k Pvel, interpolate
        kmin = Pvel[0, 0]
        f = interp1d(Pvel[:, 0], Pvel[:, 1])
    

    if kmin is not None:
        ikmin = math.ceil(kmin/dk)
    else:
        ikmin = 1

    pk[ikmin:] = f(k[ikmin:])
    
    if nc % 2 == 0:
        pk[nc // 2] = 0.0

    pksum = np.sum(pk)

    if sigma_u:
        sigma_u_raw = math.sqrt(2.0/boxsize*pksum)
        rescale = (sigma_u/sigma_u_raw)
    else:
        rescale=1.0

    mag = rescale*np.sqrt(0.5*boxsize*pk)
            
        
    uk = mag*(np.random.normal(0.0, 1.0, nk)
                    + 1j*np.random.normal(0.0, 1.0, nk))    
    uk[0]= 0.0;

    ux = np.fft.irfft(uk)*(nc/boxsize)

    return ux, rescale**2
