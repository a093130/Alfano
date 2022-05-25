# -*- coding: utf-8 -*-
"""
@description: This script plots Tau (delta-V) as a function of u, for various orbit ratios.
Reference "Optimal Many-revolution Orbit Transfer,\" Alfano & Wiesel 1985.

@version 0.1a1

@author: Colin Helms

@copyright Freelance Rocket Science, 2022

@license
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>

@change Log:
    02/23/2018: baseline
    04/29/2019: update to use alfano distribution
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
from AlfanoStateEqLib import eom_da
from AlfanoStateEqLib import eom_di

def alfano_tau(P: np.ndarray, sma = 6.6, mu = 1) -> np.ndarray:
    """
    This function computes the required delta-v to achieve a given Semi-major
    axis, expressed as an orbit ratio.
    The function is based on integration of the inverse of equation 10 
    in Alfano & Wiesel.
    
    Parameters:
    P: an array of a function of the elliptical integral of the first kind,
    appearing in the equations of motion for SMA change.
    sma: the orbit ratio in canonical units, defaults to standard GEO/LEO ratio
    mu: the canonical gravitational constant GM for the central body,
    by convenion defaults to 1 in canonical units.
    
    Returns:
        Tau, delta-v
    """
    forbit = 1- 1/np.sqrt(sma)
    Pinvert = 1/P
    
    return 0.25 * np.pi * np.sqrt(mu) * forbit * Pinvert

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import matplotlib.transforms as mtransforms
    
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
    
        tau = alfano_tau(Pu, ratio)
        plt.plot(u, tau, 'b-')
        plt.text(0.1, tau[0], 'R=%d' % ratio, fontsize=8, transform=trans_offset)
    else:
        plt.title('Change of Tau with Orbit Ratio as Parameter')
        plt.xlabel('u')
        plt.ylabel('Tau (m/s)')
    
    ax = plt.subplot(2,1,2)
    trans_offset = mtransforms.offset_copy(ax.transData, fig=fig, x=-0.32, y=-0.05, units='inches')
    
    for a in range(2,57,5):
        ratio = a
        if a > 2:
            ratio = a -2
    
        di = eom_di(Ru, ratio)
        plt.plot(u, di, 'k-')
        plt.text(0.1, di[0], 'R=%d' % ratio, fontsize=8, transform=trans_offset)
    else:
        plt.title('Change of Incl. with Orbit Ratio as Parameter')
        plt.xlabel('u')
        plt.ylabel('di/dTau (rad/m/s)')
        plt.show()
    
    plt.close()
