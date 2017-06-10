import numpy as np
import math
import scipy.integrate as integrate

def D_growth_factor(omega_m, a=None, *, z=None):
    if a is None:
        if z is None:
            raise TypeError('provide either a or z')
        
        a = 1.0/(1.0 + z)

    print('a=', a)
    print('omega_m=', omega_m)
    D0 = _growth_unnormalised(1.0, omega_m)

    return _growth_unnormalised(a, omega_m)/D0

def f_growth_rate(omega_m, a=None, *, z=None):
    """Linear growth rate f = dln D/dln a"""

    if a == 0.0:
        return 1.0
  
    d_un = _growth_unnormalised(a, omega_m);
    hf = _hubble_function(a, omega_m);

    return 1.0/(d_un*a*a*hf*hf) - 1.5*omega_m/(hf*hf*a*a*a)


def _hubble_function(a, omega_m):
    # H/H0= sqrt(Omega_m0*a^-3 + Omega_Lambda)
    return math.sqrt(omega_m/(a*a*a) + (1.0 - omega_m))

def _growth_integrand(a, omega_m):
    aHinv= 1.0/np.sqrt(omega_m/a + (1 - omega_m)*(a*a))
    return aHinv**3

def _growth_unnormalised(a, omega_m):
    """Compute unnormalised growth factor
       D(a) \propto \int_0^a (a H(a)/H0)^-3 da
    """

    result = integrate.fixed_quad(_growth_integrand, 0, a, args=(omega_m,), n=4)

    return _hubble_function(a, omega_m) * result[0]
