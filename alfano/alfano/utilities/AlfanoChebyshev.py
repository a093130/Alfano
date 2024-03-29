#! python
# # -*- coding: utf-8 -*-
"""
Created on Sat Apr 28 11:11:52 2018

@version 0.2a1

@description This module computes the control variable for an orbit transfer
with combined inclination and altitude change.
The algorithm is the Chebychev polynomial approximation from Appendix A of 
Wiesel and Alfano, "Optimal Many-Revolution Orbit Transfer".
Also see equation 6-55, Vallado "Fundamentals of Astrodynamcs", 4th ed.

@author: Colin Helms
@author_email: colinhelms@outlook.com

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

Change Log:
    05/11/2018, Correction to z_from_x(), ta = np.linspace()
    05/11/2019, Deprecated.
"""
import numpy as np
from alfano.utilities import AlfanoLib as alf
from numpy.polynomial.polynomial import polyval

def z_from_x(i_costate, a_current, mu = 1):
    """ Function computes the basis of the Chebyshev polynomial approximation
    from the current orbit ratio.  Parameter i_costate is an optimization constant.
    Reference equation 17, equation A1, Wiesel and Alfano.
    
    Use canonical variables, DU*, TU*, SU* where MU is defined as 1.
    
    Parameters:
        i_costate: the optimization constant
        a_current: the current orbit ratio
        
    Returns:
        Value of the argument to 
    """
    
    return a_current/mu * np.square(i_costate * 2/np.pi )
    
def alfano_cv (i_costate, a_current):
    """ Function uses the coefficients of the Chebyshev polynomial given in
    equation A-1, Wiesel and Alfano, 1985 to compute the approximate value 
    of the yaw law control variable .  
    
    Note that Vallado also includes this equation, however there is errata.
    The original Wiesel and Alfano paper is correct.
	
    Use canonical variables, DU*, TU*, SU*, (MU = 1), a_current as DU*.
    
    Parameters:
        i_costate: the optimization constant
        a_current: the current orbit ratio
    """
    
    alpha_coef = ([ \
    0.0,\
    2.467410607,\
    -1.907470562,\
    35.892442177,\
    -214.672979624,\
    947.773273608,\
    -2114.861134906,\
    2271.240058672,\
    -1127.457440108,\
    192.953875268,\
    8.577733773
    ])
    
    beta_coef = ([ \
    1.0,\
    0.4609698838,\
    13.7756315324,\
    -69.1245316678,\
    279.0671832500,\
    -397.6628952136,\
    -70.0139935047,\
    528.0334266841,\
    -324.9303836520,\
    20.5838245170,\
    18.8165370778
    ])
    
    z = z_from_x(i_costate, a_current)
    
    return polyval(z, alpha_coef)/polyval(z, beta_coef)

if __name__ == "__main__":
    """
    Test cases for AlfanoLib
    """    
    import matplotlib as mpl
    from mpl_toolkits.mplot3d import Axes3D 
    import matplotlib.pyplot as plt
    
    mpl.rcParams['legend.fontsize'] = 10
    
    a = np.linspace(1, 10, 100)
    costates = np.linspace(-0.1, -0.8, 100)
    ta = np.linspace(0, np.pi, 100)
    
    rad2angle = 180/np.pi
       
    fig1 = plt.figure()
    ax = Axes3D(fig1)

    """ Plot the z domain """
    z = z_from_x(costates, a)
    ax.plot3D(a, costates, z)
    ax.set_xlabel('orbit ratio')
    ax.set_ylabel('costate')
    ax.set_zlabel('z-domain')   
    plt.show()
    plt.close()
  
    fig = plt.figure(figsize=(10,10))
    ax = plt.subplot(2, 1, 1)
    
    """ Plot particular case of costate """
    cv = alfano_cv(-0.54, a)
    plt.plot(a, cv)
    plt.title('cv plot with costate -0.54')
    plt.xlabel('a')
    plt.ylabel('cv')
    plt.show()
    plt.close()

    """ Plot particular case of orbit ratio """
    cv = alfano_cv(costates, 10)
    plt.plot(costates, cv)
    plt.title('cv plot with orbit ratio 10')
    plt.xlabel('costates')
    plt.ylabel('cv')
    plt.show()
    plt.close()

    """ Plot particular case of orbit ratio """
    cv = alfano_cv(costates, 6.13)
    plt.plot(costates, cv)
    plt.title('cv plot with orbit ratio 6.13')
    plt.xlabel('costates')
    plt.ylabel('cv')
    plt.show()
    plt.close()

    """ Plot particular case of orbit ratio """
    cv = alfano_cv(costates, 1.1)
    plt.plot(costates, cv)
    plt.title('cv plot with orbit ratio 1.1')
    plt.xlabel('costates')
    plt.ylabel('cv')
    plt.show()
    plt.close()

    """Plot the control angle """
    theta = rad2angle * alf.yaw_angle(ta, alfano_cv(-0.54, 1.1))
    plt.plot(ta, theta)
    plt.title('Yaw angle for lambda -0.54 and orbit ratio 1.1')
    plt.xlabel('true anomaly')
    plt.ylabel('yaw angle')
    plt.show()
    plt.close()

    theta_h = rad2angle * alf.yaw_angle(ta, alfano_cv(-0.22, 1.1))
    plt.plot(ta, theta_h)
    plt.title('Yaw angle for lambda -0.22 and orbit ratio 1.1')
    plt.xlabel('true anomaly')
    plt.ylabel('yaw angle')
    plt.show()
    plt.close()

    theta = rad2angle * alf.yaw_angle(ta, alfano_cv(-0.54, 6.13))
    plt.plot(ta, theta_h)
    plt.title('Yaw angle for lambda -0.54 and orbit ratio 6.13')
    plt.xlabel('true anomaly')
    plt.ylabel('yaw angle')
    plt.show()
    plt.close()

    theta = rad2angle * alf.yaw_angle(ta, alfano_cv(-0.22, 6.13))
    plt.plot(ta, theta_h)
    plt.title('Yaw angle for lambda -0.22 and orbit ratio 6.13')
    plt.xlabel('true anomaly')
    plt.ylabel('yaw angle')
    plt.show()
    plt.close()
  