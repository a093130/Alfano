# -*- coding: utf-8 -*-
"""
@name: AlfanoLib
@version 0.4a

Created in Spyder Editor
This Module contains library functions derived from 
"Optimal Many-revolution Orbit Transfer," Alfano & Wiesel 1985.

@description: The following functions are incorporated in this library.
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

@author: Colin Helms

@change log:
23 Jan 2018 - Error in sign corrected in function alfano_Pprime.
26 Apr 2018 - costate correction.
05 May 2018 - found error in derivative_cmp_ell_int_1st,
    corrected complex conjugate.
06 May 2018 - added 3D plot of costate versus u and orbit ratio, R.
29 Mar 2019 - corrected alfano_phi() to return phi() rather than its reciprocal.
    updated costate to account for this change.
28 Apr 2019 - version 0.2a deployed to C:/Users/chelm/Anaconda3/Lib/site-packages
    for integration of new, complete costate table with GMAT.
21 May 2019 - version 0.3a, factored out lin_interp() to AlfanoLib.  Used in YawAngles and GenerateControlTable.
21 Jun 2019 - version 0.4a, fixed computational errors.
"""
import logging
import traceback
import platform
import getpass
import numpy as np
#import base64
import json as js
import matplotlib as mpl
import matplotlib.pyplot as plt
#    import matplotlib.transforms as mtransforms
from scipy import special
from mpl_toolkits.mplot3d import Axes3D

PREC_ORBITR = 0.01
""" Precision in orbit ratio must be the same as GenerateControlTable.py """
PREC_COSTATE = 0.001
""" Precision in lambda must be the same as GenerateControlTable.py """
MU = 1
""" Use canonical values for orbit elements.  Solve for time and velocity using mu = 1 """
nrows = int(round(1 + (10 - 1)/PREC_ORBITR, 0))
""" Orbit Ratio varies from 1 - 10 """
ncols = int(round((np.pi/2 - 0.1)/PREC_COSTATE, 0))
""" Lambda, l, varies from 0.1 to pi/2 """
halfpi = np.pi/2

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
    
    sq_compl_u = 1 - np.square(u)
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
    
    su = np.sqrt(u)
    #Fix 18Jun2019, was 1/u*E, is 1/su*E
    return (1/su) * E + (su - 1/su) * K
    #return (1/u) * E + (u - 1/u) * K
    #after Edelbaum

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
    
    sr_u = np.sqrt(u)
    #sq_u = np.square(u)
    #Fix power function, was np.power(a, 3) is np.power(a, 3/2)
    u_to_3halves = np.power(u, 3/2)
    #return 0.5*(1/a + 1/u_to_3halves)*K - (1/sq_u)*E + (1/u)*dE - 0.5*(1/a)*dK
    return -(1/u_to_3halves)*E + (1/sr_u)*dE + 0.5*(1/sr_u + 1/u_to_3halves)*K - (sr_u - 1/u)*dK
    #After Edelbaum

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
    
    return (P*(dR/dP) - R)

def costate(phi: np.ndarray, sma = 6.6, mu = 1) -> np.ndarray:
    """
    This function returns an array of values of the Lagrangian multiplier as 
    costate for a given orbit ratio.  
    
    Parameters:
    inv_phi: substitution polynomial of elliptic integrals.
    sma: the orbit ratio in canonical units, defaults to standard GEO/LEO ratio
    mu: the canonical gravitational constant GM for the central body,
    by convenion defaults to 1 in canonical units.
    
    Returns:
    Array of lambda values as a function of phi
    """     
    return ((np.pi/2) * np.sqrt(mu/sma) * 1/phi)

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

def dumps(*args, **kwargs):
    """ Overload dumps to use NumpyEncoder. """
    kwargs.setdefault('cls', NumpyEncoder)
    return js.dumps(*args, **kwargs)

def dump(*args, **kwargs):
    """ Overload dump use NumpyEncoder. """
    kwargs.setdefault('cls', NumpyEncoder)
    return js.dump(*args, **kwargs)

def loads(*args, **kwargs):
    """ Overload loads to use the decoder callback. """
    kwargs.setdefault('object_hook', cb_json_to_ndarray)    
    return js.loads(*args, **kwargs)

def load(*args, **kwargs):
    """ Overload load to use the decoder callback. """
    kwargs.setdefault('object_hook', cb_json_to_ndarray)
    
    return js.load(*args, **kwargs)

def cb_json_to_ndarray(dct):
    """ Callback to decode a JSON encoded numpy ndarray.
    The shape and dtype is stored in the dictionary returned from NumpyDecoder
    Returns ndarray.
    
    base64.b64decode credit to Adam Hughes, and hpaulj on Stack Overflow,
    https://stackoverflow.com/questions/27909658/
    json-encoder-and-decoder-for-complex-numpy-arrays/27948073#27948073
    Status: not working, as a workaround, Controls.json stores a python list 
    rather than ndarray.

    Parameters:
        dct: (dict) json encoded ndarray
    """
#    if isinstance(dct, dict) and '__ndarray__' in dct:
#        data = base64.b64decode(dct['__ndarray__'])        
#        return np.frombuffer(data, dct['dtype']).reshape(dct['shape'])
#        return np.frombuffer(data, dtype=dct['dtype']).reshape(dct['shape'])
    
    if (dct, dict):
        return(dct)

def lin_interp(l_hi, l_lo, lamb, u_lo, u_hi):
    """ linear interpolation returning value between u_hi and u_lo proportional to 
    lamb between l_hi and l_lo.
    
    Calling procedure may need to round.
    """
    return u_lo + (u_hi - u_lo) * (l_hi - lamb)/(l_hi - l_lo)

try:
    """ Compute Global arrays - this constitutes an interface agreement with AlfanoLib"""
    u = np.round(0.1 * np.linspace(1, 10, ncols, endpoint=False), 4)
    """ This is the 1470 element domain of cv from 0 - 1, excluding singularities. """
    a = np.round(np.linspace(1, 10, nrows), 2)
    """ This is the 901 element domain of orbit ratio. """
    
    vec_k = cmp_ell_int_1st_kind(u)
    vec_e = cmp_ell_int_2nd_kind(u)
    vec_dk = derivative_cmp_ell_int_1st(u, vec_k, vec_e)
    vec_de = derivative_cmp_ell_int_2nd(u, vec_k, vec_e)
    vec_p = alfano_P(u, vec_k)
    vec_dp = alfano_Pprime(u, vec_k, vec_dk)
    vec_r = alfano_R(u, vec_k, vec_e)
    vec_dr = alfano_Rprime(u, vec_r, vec_e, vec_dk, vec_de)  
    vec_phi = alfano_phi(vec_r, vec_p, vec_dr, vec_dp)
    
    Lambda = np.round(halfpi * (MU/np.sqrt(1)) * 1/vec_phi, 4)
    """ Canonical Lambda - all other values are multiples of this row vector. """
    
    UbyRbyL = {l: np.zeros(nrows) for l in Lambda}
    """ This dictionary is the main interface between YawAngles.py nad GenerateControlTable.py 
    The structure of this dictionary is {lambda:array(cv)},
    where cv is in order of orbit ratio, R.
    
    It is important that undefined elements of u are zero because  
    in the steering law this means tan(u) = 0, no yaw.
    """
    
except Exception as e:
    lines = traceback.format_exc().splitlines()
    
    logging.error("Exception: %s, %s, %s\n%s\n%s", e.__class__, e.__doc__, e.__cause__, lines[0], lines[-1])
    print("Error in AlfanoLib:{0}".format(e.__doc__))
    

class NumpyEncoder(js.JSONEncoder):
    """ Encode a Numpy ndarray. """
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            """ Object is ndarray. Convert into a dict holding 
            dtype, shape and the data, base64 encoded.
            """
            #data_b64 = base64.b64encode(np.ascontiguousarray(obj).data)
            #data_b64 = io.StringIO(base64.b64encode(np.ascontiguousarray(obj).data))
            """ Proximous fix for a bug in Adam Hughes' code """
            
            #return dict(__ndarray__=data_b64, dtype=str(obj.dtype), shape=obj.shape)
            return obj.tolist()
        
        logging.debug("NumpyEncoder fell through for type {0}.".format(type(obj)))
        super(NumpyEncoder, self).default(obj)

    """ Thanks to https://stackoverflow.com/users/3768982/tlausch on Stack Overflow,
    https://stackoverflow.com/questions/3488934/simplejson-and-numpy-array/24375113#24375113
    """

if __name__ == "__main__":
    """
    Test case for AlfanoLib
    """    
    __spec__ = None
    """ Necessry tweak to get Spyder IPython to execute this code. 
    See:
    https://stackoverflow.com/questions/45720153/
    python-multiprocessing-error-attributeerror-module-main-has-no-attribute
    """
    
    logging.basicConfig(
            filename='./AlfanoLib.log',
            level=logging.INFO,
            format='%(asctime)s %(filename)s %(levelname)s:\n%(message)s', datefmt='%d%B%Y_%H:%M:%S')

    logging.info("!!!!!!!!!! AlfanoLib Plots Started !!!!!!!!!!")
    
    host_attr = platform.uname()
    logging.info('User Id: %s\nNetwork Node: %s\nSystem: %s, %s, \nProcessor: %s', \
                 getpass.getuser(), \
                 host_attr.node, \
                 host_attr.system, \
                 host_attr.version, \
                 host_attr.processor)

    logging.info('Plotspace consists of {0} rows and {1} columns.'.format(nrows, ncols))
      
    K = vec_k 
    E = vec_e 
    dK = vec_dk 
    dE = vec_de 
    P = vec_p 
    dP = vec_dp 
    R = vec_r 
    dR = vec_dr   
    phi = vec_phi 
    
    logging.info('Plot instance of Phi(u):\n{0}'.format(phi))
   
    costates = costate(phi)
        
    logging.info('Plot instance of costates:\n{0}'.format(costates))
    
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
    
    fig = plt.figure(figsize=(10,10))
    ax = plt.subplot(2, 1, 1)
    
    plt.plot(u, P, 'm.')
    
    plt.ylabel('P')
    plt.xlabel('u')
    plt.title('Function P(u), u=0.1 to 1.0')
    
    ax = plt.subplot(2, 1, 2)
    
    plt.plot(u, R, 'c.')
    
    plt.ylabel('R')
    plt.xlabel('u')
    plt.title('Function R(u), u=0.1 to 1.0')
    
    plt.tight_layout()
    plt.show()
    plt.close()

    fig = plt.figure(figsize=(10,10))
    ax = plt.subplot(2, 1, 1)

    plt.plot(u, dP, 'm-')
    
    plt.ylabel('dP')
    plt.xlabel('u')
    plt.title('Function dP(u), u=0.1 to 1.0')

    ax = plt.subplot(2, 1, 2)

    plt.plot(u, dR, 'c-')
    
    plt.ylabel('dR')
    plt.xlabel('u')
    plt.title('Function dR(u), u=0.1 to 1.0')

#    fig, axs = plt.subplots(2,1)

#    Pplot, Rplot = axs[0].plot(u, P, 'm.', u, R, 'c.')
#    dPplot, dRplot = axs[1].plot(u, dP, 'm-', u, dR, 'c-')

#    fig.legend((Pplot, Rplot),('P(u)','R(u)'), 'upper right')
#    fig.legend((dPplot, dRplot), ('dP(u)/du', 'dR(u)/du'), 'right')
#    plt.xlabel('u')
#    plt.title('P(u) and R(u)')
    
    plt.tight_layout()
    plt.show()
    plt.close()

    fig = plt.figure(figsize=(10, 8))
    
    plt.plot(u, phi, 'bo')
    
    plt.ylabel('Phi')
    plt.xlabel('u')
    plt.title('Function Phi(u), u=0.1 to 1.0')

    plt.tight_layout()
    plt.show()
    plt.close()

#    fig = plt.figure(figsize=(10,10))
#    ax = plt.subplot(2, 1, 1)
    
#    costates=costate(0.25, a)
    
#    plt.plot(u, phi, 'bo')
#    plt.ylabel('costate')
#    plt.xlabel('Orbit Ratio')
#    plt.title('Costate, u=0.25, R=1-10') 
    
#    ax = plt.subplot(2, 1, 2)
#    trans_offset = mtransforms.offset_copy(ax.transData, fig=fig, x=-0.5, y=-0.05, units='inches')
    
#    costates=costate(phi)
    
#    plt.plot(u, costates, 'b.')
#    plt.ylabel('costate')
#    plt.xlabel('control variable')
#    plt.title('Costate, u=0-1, R=6.6')
    
#    plt.tight_layout()
#    plt.show()
#    plt.close()
    
   
    """ Wireframe Plot of Costate """

    #Make arrays for axes  
    L = np.zeros((nrows, ncols))
    X = np.zeros((nrows, ncols))
    Y = np.zeros((nrows, ncols))
    Yt = np.zeros((nrows, ncols))

    """
    TODO: Neat trick to transpose the elements of a list from 
    https://docs.python.org/3/tutorial/datastructures.html?highlight=transpose
    zip(*X) makes an iterator that aggregates elements.  We don't use it here.
    """
  
    #Express the costates as functions of u and a, where i in a and j in u
    inx=-1
    for sma in a:
        iny = -1
        inx = inx + 1
        for cv in u:
            iny = iny + 1
            X[inx, iny] = cv
            Y[inx, iny] = sma
                
            k = cmp_ell_int_1st_kind(cv)
            e = cmp_ell_int_2nd_kind(cv)
            dk = derivative_cmp_ell_int_1st(cv, k, e)
            de = derivative_cmp_ell_int_2nd(cv, k, e)
            p = alfano_P(cv, k)
            dp = alfano_Pprime(cv, k, dk)
            r = alfano_R(cv, k, e)
            dr = alfano_Rprime(cv, r, e, dk, de)
            f = alfano_phi(r, p, dr, dp)
            
            L[inx, iny] = costate(f, sma)
                       
    fig3d2 = plt.figure(figsize=(12, 8))
    ax = Axes3D(fig3d2)
    ax.plot_wireframe(X, Y, L)
    
    ax.set_xlabel('control variable u')
    ax.set_ylabel('orbit ratio')
    ax.set_zlabel('lambda_i') 
    plt.show()
    plt.close()
