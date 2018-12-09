# -*- coding: utf-8 -*-
"""
Created in Spyder Editor
This Module contains library functions taken from 
"Optimal Many-revolution Orbit Transfer," Alfano & Wiesel 1985.

Functions
    cmp_ell_int_1st_kind is the complete elliptical integral of the first kind, K(u).
    cmp_ell_int_2nd_kind is the complete elliptical integral of the second kind, E(u).
    derivative_cmp_ell_int_1st is the derivative of K(u) wrt u.
    derivative_cmp_ell_int_2nd is the derivative of E(u) wrt u.
    alfano_P, is a substitution, P(u), containing terms of K(u).
    alfano_Pprime, the derivative of P(u) wrt u.
    alfano_R, is a substitution, R(u), containing terms of E(u) and K(u).
    alfano_Rprime, the derivative of R(u) wrt u.
    Note that Alfano & Wiesel is incorporate in Vallado, "Fundamentals of
        Astrodynamics," Section 6.7, where cv := u.
    The costate for Semi-Major Axis (SMA) is solved in terms of the costate for inclination resulting 
        in phi, a function of elliptic integrals and their derivatives.
    The costate for inclination is constant and its solution is optimal.
    The inverse of phi() is used to select the optimum value for the
        control variable, u. The optimal costate must be known.
    Given the current SMA and inclination, the u value solves the Edelbaum 
        control law for the many-revolution combined maneuver in inclination and  
        orbit raising.

If this module is executed as a script, the values of the above are computed
for a linear distribution of u, and figures containing plots are generated
for documentation.  

The costate uses the orbit ratio as a parameter, this functions is plotted
parametrically with R. 

Author: colinhelms@outlook.com

Version 1.12

Change Log:
(1.1) 23 Jan 2018 16:22 Error in sign corrected in function alfano_Pprime.
(1.12) 26 Apr 2018 13:14, costate correction.
(1.13) 05 May 2018 20:05, found error in derivative_cmp_ell_int_1st,
    corrected complex conjugate.
(1.14) 06 May 2018, added 3D plot of costate versus u and orbit ratio, R.  
"""
import numpy as np
from scipy import special

def cmp_ell_int_1st_kind(u: np.ndarray) -> np.ndarray:
    """
    K(u)
    Parameters:
    u: an array of floating point values between 0 and 1,

    Returns an array for the evaluation on the complete elliptical
    integral of the first kind around each of those points.
    """
    return special.ellipk(u)

def cmp_ell_int_2nd_kind(u: np.ndarray) -> np.ndarray:
    """
    E(u)
    Parameters:
    u: an array of floating point values between 0 and 1,

    Returns an array for the evaluation on the complete elliptical
    integral of the second kind around each of those points.
    """
    return special.ellipe(u)

def derivative_cmp_ell_int_1st(u: np.ndarray, K: np.ndarray, E: np.ndarray) -> np.ndarray:
    """
    K'(u), derivative of the complete elliptic integral of the first kind.
    
    Parameters:
    u: an array of floating point values between 0 and 1,
    K: an array of values of the complete elliptical integral of
    the first kind, evaluated at u
    E: an array of values of the complete elliptical integral of
    the second kind, evaluated at u

    Returns the first derivative of the complete elliptical
    integral of the first kind around each of those points, u.
    
    The integral is evaluated using the definition of the first
    derivative from dlmf.nist.gov, equ. 19.4.1 for the derivative of K:
        dK/du = [E - (u'**2)*K]/[u*(u'**2)]
        where u'**2 = sqrt(u**2) = u. 
    The complex conjugate of areal number is itself.
    """
    #if u.any() == 0 or u.any() == 1: return float(np.NaN)
    
    sq_compl_u = np.square(u)
    return (E - sq_compl_u * K)/(u*sq_compl_u)

def derivative_cmp_ell_int_2nd(u: np.ndarray, K: np.ndarray, E: np.ndarray) -> np.ndarray:
    """
    E'(u), derivative of the complete elliptic integral of the second kind.
    
    Parameters:
    u: an array of floating point values between 0 and 1,
    E: an array of values of the complete elliptical integral of
    the second kind, evaluated at u.
    
    Returns the first derivative of the complete elliptical
    integral of the second kind around each of the points, u.
    
    The integral is evaluated using the definition of the first
    derivative from dlmf.nist.gov, equ. 19.4.2 for the derivative of E:
        dE/du = [E - K]/u
    """
    #if u.any() == 0: return float(np.NaN)
    
    return (E - K)/u
    
def alfano_P(u: np.ndarray, K: np.ndarray) -> np.ndarray:
    """
    P(u) appears in the equations of motion and Hamiltonian, for the many
    revolution problem.  It is a factor in the Alfano Phi function.

    Parameters:
    u: an array of floating point values between 0 and 1,
    K: an array of values of the complete elliptical integral of
    the first kind, evaluated at u
    
    Returns P(u) evaluated around each of points u.
    P(u) is the substitution polynomial in K from Alfano equation 10:
        P = [(1-u)**1/2]*K
    """
    
    a = np.sqrt(1 - u)
    return  a * K

def alfano_Pprime(u: np.ndarray, K: np.ndarray, dK: np.ndarray) -> np.ndarray:
    """
    P'(u) appears as a factor in the Alfano Phi function.
    
    Parameters: 
    u: an array of floating point values between 0 and 1,
    K: an array of values of the complete elliptical integral of
    the first kind, evaluated at u,
    dK: an array of values of the derivative of the CE integral of the first kind,
    evaluated at u.

    Returns dP(u)/du around each of points, u.
    Derivation of Pprime is the original work of the author.
    dP = -[1/2*1/sqrt(1-u)]*K + sqrt(1-u)*dK
    """
    #if u.any() == 1: return float(np.NaN)
    
    a = np.sqrt(1 - u)
    return (0.5/a)*K + a*dK

def alfano_R(u: np.ndarray, K: np.ndarray, E: np.ndarray) -> np.ndarray:
    """
    R(u) appears as a factor in the Alfano Phi function.
    
    Parameters: 
    u: an array of floating point values between 0 and 1,
    K: an array of values of the complete elliptical integral of
    the first kind, evaluated at u,
    E: an array of values of the complete elliptical integral of
    the second kind, evaluated at u.

    Returns R(u) around each of points, u.
    R(u) is the substitution polynomial in K from Alfano equation 11:
        R = E/u + (sqrt(u) - 1/sqrt(u))*K
    """
    #if u.any() == 0 or u.any() == 1: return float(np.NaN)
    
    a = np.sqrt(u)
    return (1/u)*E + (a - 1/a)*K

def alfano_Rprime(u: np.ndarray, K: np.ndarray, E: np.ndarray, dK: np.ndarray, dE: np.ndarray) -> np.ndarray:
    """
    R'(u) appears as a factor in the Alfano Phi function.

    Parameters: 
    u: an array of floating point values between 0 and 1,
    K: an array of values of the complete elliptical integral of
    the first kind, evaluated at u,
    E: an array of values of the complete elliptical integral of
    the second kind, evaluated at u.
    dK: an array of values of the derivative of the CE integral of the first kind,
    evaluated at u.
    dE: an array of values of the derivative of the CE integral of the second kind,
    evaluated at u.
    
    Returns dR(u)/du around each of the points, u.
    Derivation of Rprime is the original work of the author.
    dR = 1/2*(1/sqrt(u) + 1/sqrt(u**3))*K - (1/u**2)*E + (1/u)*dE -1/2*(1/sqrt(u)*dK)
    """
   # if u.any() == 0: return float(np.NaN)
    
    a = np.sqrt(u)
    sq_u = np.square(u)
    u_to_3halves = np.power(a, 3)
    return 0.5*(1/a + 1/u_to_3halves)*K - (1/sq_u)*E + (1/u)*dE - 0.5*(1/a)*dK

def alfano_phi(R: np.ndarray, P: np.ndarray, dR: np.ndarray, dP: np.ndarray) -> np.ndarray:
    """
    This function returns the value of phi = P(u)*dR/dP - R(u),  Eqn. 18, 
    Alfano & Weisel, 1985
      
    Note there is a singularity where dP/dR = R(u)/P(u).
    
    Parameters:
    R: an array of a function of the elliptical integral of the first and
    second kind, appearing in the equations of motion for inclination change.
    P: an array of a function of the elliptical integral of the first kind,
    appearing in the equations of motion for SMA change.
    dR: the first derivative of R with respect to u
    dP: the first derivative of P with respect to u
    
    Returns:
    phi evaluated around the points, u.
    
    Note that phi evaluates to 0 as u -> 0.1.
    """
    
    f = P*(dR/dP) - R

    return (1/f)

def costate(phi: np.ndarray, sma = 6.6, mu = 1) -> np.ndarray:
    """
    This function returns the value of the Lagrangian multiplier as 
    costate for a given orbit ratio.  
    
    Parameters:
    inv_phi: substitution polynomial of elliptic integrals.
    sma: the orbit ratio in canonical units, defaults to standard GEO/LEO ratio
    mu: the canonical gravitational constant GM for the central body,
    by convenion defaults to 1 in canonical units.
    
    Returns:
    Array of lambda values as a function of phi
    """   
    a = np.sqrt(sma)
    m = np.sqrt(mu)
    
    return (np.pi/2)*(m/a)*phi

def yaw_scalefactor (u):
	""" Convenience function that returns the correct form of the denominator in the
		Edelbaum control law.
		
			sqrt[1/u - 1]
		
    where u is the return value from alfano_cv().
	Use canonical variables, DU*, TU*, SU*, (MU = 1), a_current as DU*.
    
    Parameters:
        u: Alfano Control Variable, AKA cv
	"""
	return np.sqrt((1/u) - 1)
	
def yaw_angle (TA, cv):
	""" Function implements the Edelbaum control law,
	
		yaw_angle = arctan[cos(TA)/sqrt(1/u - 1)]
        
    Return value: radians from -pi/2 to pi/2
		
	Parameters:
		TA: the true anomaly, or astronomical longitude for a circle orbit.
		cv: the Alfano control variable.
	"""
	sf = yaw_scalefactor (cv)
	
	return np.arctan(np.cos(TA)/sf)

if __name__ == "__main__":
    """
    Test case for AlfanoLib
    """    
    import matplotlib as mpl
    from mpl_toolkits.mplot3d import Axes3D 
    import matplotlib.pyplot as plt
    import matplotlib.transforms as mtransforms
        
    size=100
    u = np.linspace(0.100, 0.85, size)
    a = np.linspace(1, 10, size)

    K = cmp_ell_int_1st_kind(u)
    E = cmp_ell_int_2nd_kind(u)
    dK = derivative_cmp_ell_int_1st(u, K, E)
    dE = derivative_cmp_ell_int_2nd(u, K, E)
    P = alfano_P(u, K)
    dP = alfano_Pprime(u, K, dK)
    R = alfano_R(u, K, E)
    dR = alfano_Rprime(u, R, E, dK, dE)
    phi = alfano_phi(R, P, dR, dP)
    costates = costate(phi, a)

    mpl.rcParams['legend.fontsize'] = 10
    
    fig, axs = plt.subplots(2,1)

    Kplot, Eplot = axs[0].plot(u, K, 'r.', u, E, 'b.')
    dKplot, dEplot = axs[1].plot(u, dK, 'r-', u, dE, 'b-')

    fig.legend((Kplot, Eplot), ('K(u)', 'E(u)'), 'upper right')
    fig.legend((dKplot,dEplot), ('dK(u)/du', 'dE(u)/du'), 'right')
    plt.xlabel('u')
    plt.title('Cmpl Elliptic Integrals')
    
    plt.tight_layout()
    plt.show()
    plt.close()

    fig, axs = plt.subplots(2,1)

    Pplot, Rplot = axs[0].plot(u, P, 'm.', u, R, 'c.')
    dPplot, dRplot = axs[1].plot(u, dP, 'm-', u, dR, 'c-')

    fig.legend((Pplot, Rplot),('P(u)','R(u)'), 'upper right')
    fig.legend((dPplot, dRplot), ('dP(u)/du', 'dR(u)/du'), 'right')
    plt.xlabel('u')
    plt.title('P(u) and R(u)')
    
    plt.tight_layout()
    plt.show()
    plt.close()

    fig = plt.figure(figsize=(10,10))
    ax = plt.subplot(2, 1, 1)
    costates=costate(0.25, a)
    plt.plot(u, phi, 'bo')
    plt.ylabel('costate')
    plt.xlabel('Orbit Ratio')
    plt.title('Costate, u=0.25, R=1-10') 
    
    ax = plt.subplot(2, 1, 2)
    trans_offset = mtransforms.offset_copy(ax.transData, fig=fig, x=-0.5, y=-0.05, units='inches')
    
    costates=costate(phi)
    plt.plot(u, costates, 'b.')
    plt.ylabel('costate')
    plt.xlabel('control variable')
    plt.title('Costate, u=0-1, R=6.6')
    
    plt.tight_layout()
    plt.show()
    plt.close()
    
   
    """ Wireframe Plot of Costate """

    #Make arrays for axes  
    L = np.zeros((100,100))
    X = np.zeros((100,100))
    Y = np.zeros((100,100))
    Yt = np.zeros((100,100))

    """
    Neat trick to transpose the elements of a list from 
    https://docs.python.org/3/tutorial/datastructures.html?highlight=transpose
    zip(*X) makes an iterator that aggregates elements.  We don't use it here.
    """
  
    #Express the costates as functions of u and a, where i in a and j in u
    inx=-1
    iny=-1
    for sma in a:
        inx=inx+1
        iny=-1
        for cv in u:
            iny=iny+1
            
            X[inx,iny] = cv
            Y[inx,iny] = sma
                
            k = cmp_ell_int_1st_kind(cv)
            e = cmp_ell_int_2nd_kind(cv)
            dk = derivative_cmp_ell_int_1st(cv, k, e)
            de = derivative_cmp_ell_int_2nd(cv, k, e)
            p = alfano_P(cv, k)
            dp = alfano_Pprime(cv, k, dk)
            r = alfano_R(cv, k, e)
            dr = alfano_Rprime(cv, r, e, dk, de)
            f = alfano_phi(r, p, dr, dp)
            
            L[inx,iny]=costate(f, sma)
    
    fig3d2 = plt.figure()
    ax = Axes3D(fig3d2)
    ax.plot_wireframe(X, Y, L)
    
    ax.set_xlabel('ctl variable')
    ax.set_ylabel('orbit ratio')
    ax.set_zlabel('costates') 
    plt.show()
    plt.close()
