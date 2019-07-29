# -*- coding: utf-8 -*-
"""
Created in Spyder Editor
@description: This Module contains library functions taken from 
"Optimal Many-revolution Orbit Transfer,\" Alfano & Wiesel 1985.

Functions:
    eom_da, the change in semi-major axis (SMA) per delta-v increment.
    eom_di, the change in inclination per delta-v increment.

When executed rather than imported as a library,
this script plots its two functions as a function of u.  

@author: Colin Helms

@copyright Copyright Freelance Rocket Science, 2019

@change log:
    02/23/2018, baseline
    04/29/2019, updated to use AlfanoLib distribution
"""
import numpy as np
from alfano.utilities.AlfanoLib import cmp_ell_int_1st_kind
from alfano.utilities.AlfanoLib import cmp_ell_int_2nd_kind
from alfano.utilities.AlfanoLib import derivative_cmp_ell_int_1st
from alfano.utilities.AlfanoLib import derivative_cmp_ell_int_2nd
from alfano.utilities.AlfanoLib import alfano_P
from alfano.utilities.AlfanoLib import alfano_Pprime
from alfano.utilities.AlfanoLib import alfano_R
from alfano.utilities.AlfanoLib import alfano_Rprime
from alfano.utilities.AlfanoLib import alfano_phi

def eom_da(P: np.ndarray, sma = 6.6, mu = 1) -> np.ndarray:
    """
    This function computes the change in semi-major axis (SMA) per delta-v 
    increment.  See equation 10 in Alfano and Wiesel, 1985.
    
    Parameters:
    P: an array of a function of the elliptical integral of the first kind,
    appearing in the equations of motion for SMA change.
    sma: the orbit ratio in canonical units, defaults to standard GEO/LEO ratio
    mu: the canonical gravitational constant GM for the central body,
    by convenion defaults to 1 in canonical units.
    
    Returns:
    da/dtau: the differential relation between changes in semi-major axis
    with respect to the substitution variable, tau.  The tau is formed from
    the Rocket Equation and represents the change in delta-v per orbit.
    """
    a = np.sqrt(np.power(sma, 3))
    m = np.ones(P.shape)
    m = 4*a/(mu*np.pi)
    
    return m*P

def eom_di(R: np.ndarray, sma = 6.6, mu = 1) -> np.ndarray:
    """
    This function computes the change in inclination per delta-v 
    increment.  See equation 11 in Alfano and Wiesel, 1985.
    
    Parameters:
    R: an array of a function of the elliptical integral of the first and
    second kind, appearing in the equations of motion for inclination change.
    sma: the orbit ratio in canonical units, defaults to standard GEO/LEO ratio
    mu: the canonical gravitational constant GM for the central body,
    by convenion defaults to 1 in canonical units.
    
    Returns:
    di/dtau: the differential relation between changes in semi-major axis
    with respect to the substitution variable, tau.  The tau is formed from
    the Rocket Equation and represents the change in delta-v per orbit.
    """
    a = np.sqrt(sma)
    m = np.ones(R.shape)
    m = 2*a/(mu*np.pi)
    
    return m*R

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import matplotlib.transforms as mtransforms

    # dE(u)/du tends to infinity as u -> 1, and avoid divide by 0 in Rprime computaion
    u = np.arange(0.100, 0.85, 0.005)

    K = cmp_ell_int_1st_kind(u)
    E = cmp_ell_int_2nd_kind(u)
    
    dK = derivative_cmp_ell_int_1st(u, K, E)
    dE = derivative_cmp_ell_int_2nd(u, K, E)
    
    Pu = alfano_P(u, K)
    dPu = alfano_Pprime(u, K, dK)
    Ru = alfano_R(u, K, E)
    dRu = alfano_Rprime(u, Ru, E, dK, dE)
    phi = alfano_phi(Ru, Pu, dRu, dPu)
    
    fig = plt.figure(figsize=(10,10))
    ax = plt.subplot(2, 1, 1)
    
    trans_offset = mtransforms.offset_copy(ax.transData, fig=fig, x=-0.32, y=-0.05, units='inches')
    
    for a in range(2,57,5):
        ratio = a
        if a > 2:
            ratio = a -2
    
        da = eom_da(Pu, ratio)
        plt.plot(u, da, 'k-')
        plt.text(0.1, da[0], 'R=%d' % ratio, fontsize=8, transform=trans_offset)
    else:
        plt.title('Change of SMA with Orbit Ratio as Parameter')
        plt.xlabel('u')
        plt.ylabel('da/dTau (m/m/s)')
    
    ax = plt.subplot(2,1,2)
    trans_offset = mtransforms.offset_copy(ax.transData, fig=fig, x=-0.32, y=-0.05, units='inches')
    
    for a in range(2,57,5):
        ratio = a
        if a > 2:
            ratio = a - 2
    
        di = eom_di(Ru, ratio)
        plt.plot(u, di, 'k-')
        plt.text(0.1, di[0], 'R=%d' % ratio, fontsize=8, transform=trans_offset)
    else:
        plt.title('Change of Incl. with Orbit Ratio as Parameter')
        plt.xlabel('u')
        plt.ylabel('di/dTau (rad/m/s)')
        plt.show()
    
    plt.close()
