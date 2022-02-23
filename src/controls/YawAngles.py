#! python
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 28 11:11:52 2018
@Version 0.3a

@Description:
    This module performs a lookup on a precomputed table of control law scale factors 
    by orbit ratio.  It is designed to be called from within a GMAT mission sequence.

    The result of the lookup is multiplied by the cosine of the Argument of Latitude and the 
    yaw angle is returned in degrees.
    
    Yaw is defined here as as an angle normal to the orbital plane.
    This can be confusing as some older authors define this angle as pitch, following the
    convention of aircraft wherein pitch is normal to the horizon.

    The computation of the yaw control angle is based on Edelbaum's control law, see
    Wiesel and Alfano, "Optimal Many-Revolution Orbit Transfer" and is implemented 
    in package Alfano.AlfanoLib.py.

    The argment to the control law, the control variable, is read from a 901x1471 dictionary
    formatted as a JSON file:
    
    {
    Lambda_0:[cv1[0], cv1[1], ..., cv1[900]], 
    Lambda_1:[cv2[0], cv2[1], ..., cv2[900]], 
    ..., 
    Lambda_1470:[cv1470[0], cv1470[1], ..., cv1470[900]]
    }
    
    where the cv are in order of orbit ratio from 1.00 to 10.00.

    Internal Dependencies:
        get_control_onrev() ->
            get_control() ->
                get_yaw_angle() -> 
                    get_yaw_cv()  

@author: Colin Helms

@copyright Copyright Freelance Rocket Science, 2018, 2019

@change log:
    04/28/2018: Initial baseline.
    04/30/2018: AlfanoChebyshev.py used to generate u values.
    05/11/2018, Rewritten to use a pre-generated linear array of u values.
    05/18/2018, Added test cases and adjusted signatures for GMAT call.
    05/20/2018, Changed name of input files, integrated with GMAT.
    06/01/2018, Added logic for Kluever eclipse weighting, and test case.
    03/06/2019, Added costate parameter, enables variable inclinations under external control.
    04/08/2019, Controls.json elaborated for all costates. 
    04/10/2019, Added JSON codec for dictionary of ndarray.
    04/28/2019, Deployed version 0.2a to C:/Users/chelm/Anaconda3/Lib/site-packages
    for integration of new Controls.json with GMAT.
    05/21/2019, Renamed get_yaw_sf() to get_yaw_cv(). Added interpolation of
    input costate value.

"""
import numpy as np
#import scipy.special as sp
import json as js
from alfano import AlfanoLib as alf

_fp = None
""" file pointer in Global scope """


def get_ecc_control (TA, only = True, more = -1):
    """ Per Edelbaum equ. 22, change eccentricity without change to SMA or Inclination.
    
    Parameters:
        The argument TA must be True Anomaly, measured from the perigee.
        The default value eonly = True provides and eccentricity change only,
            a Value of False will provide a combined maneuver with delta-V approximately
            Vini * 0.649 * delta-ecc.
        The default value more = -1 will decrease eccentricity.
    """
    
    phase = np.radians(TA)
        
#    B = more * np.sin(phase)/abs(np.sin(phase))
#    alpha = alf.halfpi * np.tan(phase)/abs(np.tan(phase))
    
    if only:
        """ Eccentricity change only """
#        alpha = alf.halfpi * np.sign(np.tan(0.5*phase))
        alpha = alf.halfpi * np.sign(np.sin(phase))        
    else:
        """ Combined eccentricity and altitude change """
#        alpha = np.arctan(0.5 * np.tan(phase))
        alpha = np.arctan(np.tan(0.5*phase))
#        alpha = np.arctan(0.5 * sp.cotdg(TA))
   
    return [more * np.cos(alpha), 0.0, more * np.sin(alpha)]

def get_control_onrev (costate, AOL, SMA, SMA_init = 6838.1366, more = -1):
    """ Function provides a wrapper to perform conversions from SMA in km
    to the orbit ratio; the given SMA is divided by SMA_init.  SMA_init is
    set to the canonical Earth radius to avoid divide by zero, but the 
    caller should pass in the SMA of the starting orbit to get correct results.
    
    Function checks for completed revolution and calls get_yaw_angle() with 
    given AOL and computed orbit_r, sets AOL_init equal to the current AOL, else
    returns original yaw angle.
    
    Future: eclipse effects will accounted for in function set_apsides().
    
    Returns:
        [Thrust vector components in Velocity-Normal-Bi-normal coordinates]
        
    Parameters:
        AOL: argument of longitude in degrees
        SMA: the current SMA in kilometers
        SMA_init: the SMA of the initial orbit in kilometers
        more: defines the direction of the yaw angle, -1 decreases inclination.
        costate: the lambda inclination solution.
    """
    if (costate >= -1.570) and (costate <= -0.100):
        """ costate is essentially a tangent, causes computation problems in the neighborhood of 0 """
        if SMA_init >= 6338.1366:
            orbit_r = np.round(SMA/SMA_init, 2)
        else:
            raise ValueError("SMA_init {0} is invalid.".format(SMA_init))
    else:
        raise BadCostate("Value {0} is out of range.".format(costate))
        
    if (orbit_r >= 1) and (orbit_r <= 10):
        B = 0.0
        """ For combined orbit-raising and inclination, the binormal thrust angle (pitch) is 0. """
        
        V, N = get_control(costate, AOL, orbit_r, more)
         
        return [V, N, B]
    
    else:
        raise ValueError ("Orbit ratio {0} is invalid value.".format(orbit_r))

        
def get_control (costate, AOL, orbit_r, nmore):
    """ Returns the Alfano control for each orbit increment.
    
            yaw_angle = arctan[cos(AOL)/sqrt(1/u - 1)].
        
        where u is indexed by orbit_r and costate in the controls.json file
        
        The algorithm is a reverse-lookup of precomputed values
        of the Edelbaum control variable, u (aka cv) which are computed in 
        another script using the Alfano method for solution of the 
        multi-revolution orbit transfer (see AlfanoLib.py)
        
        The precomputed values are first written into an Excel Workbook and
        a table of u is created by look-up from the optimization costate.
        This table is called the Trajectory table, because values of u for
        any given costate define trajectories in the inclination and 
        semi-major axis plane.  The Trajectory table is provided both in
        the generating Excel workbook, and in a JSON table.  
        
        The Excel workbook access may provide better performance since the
        Excel file remains open as an active object following the first call, 
        values are returned via the active object socket protocol.
        
        The JSON file method is more portable, however values for only one
        value of the optimization costate are available per file.
        
        Set useJSON = False in order to use the Excel active object.
        
        Returns:
            Components of thrust in VNB coordinates.
        
        Parameters:
            AOL: the angular position from the line of nodes, in degrees 
            orbit_r: the current orbit ratio
            nmore: +/1, defines the direction of the yaw angle, 
                defined negative for decreasing inclination.
    """
    beta = nmore * get_yaw_angle (AOL, orbit_r, costate)
    
    return [np.cos(beta), np.sin(beta)]


def get_yaw_angle (AOL, orbit_r, costate):
    """ Function implements the Edelbaum control law.  Returns the yaw angle
    in degrees.  This function is good for plots.
        
        Returns:
            Thrust angle in radians.
        
        Parameters:
            AOL: the angular position from the line of nodes, in degrees 
            orbit_r: the current orbit ratio  
    """
    
    AOL = np.radians(AOL)
    """ GMAT provides degrees, np.cos expects radians """
    
    cv = get_cv(orbit_r, costate)
    
    if cv != 1:
        sf = np.sqrt(cv/(1-cv))
        
        beta = np.arctan(np.cos(AOL) * sf)
    else:
        beta = alf.halfpi * np.sign(np.cos(AOL))
        """
        The trajectory for any given value of λ_i can be visualized as a plane parallel to the u, 
        R axes cutting through the Φ surface at a vertical offset equal to λ_i.  
        Starting at a λ_ivalue of -0.496 the trajectory is cutoff at the right edge,
        where the u value approaches 1. This agrees with Alfano and Wiesel’s original figure 
        as reprinted in Vallado Figure 6-24 [6], which shows that trajectories for λ_i< -0.5 
        are predominately inclination change.
        
        """
              
    return beta
        
def get_cv(orbit_r, costate) :
    """ Function looks up control variable in JSON Control file based on costate. 
    returns the denominator of the control function.
    Changed 08Apr2019: complete table of costates is searched.
    """
    global _fp
    
    if _fp == None:
        read_controlfile()
    
    row = int(round(1 + orbit_r/alf.PREC_ORBITR, 0)) - 100
    
    """ The given costate may not be an exact match for the key to UbyRbyL.  For
    example, a given value of -0.36 will fall between the values -0.3599 and
    -0.3604.  The array AlfanoLib.Lambda contains the exact keys.  Since the intervals
    between the values of Lambda are small in the costate dictionary (UbyRbyL),
    linear interpolation of the return value is feasible.
    """

    if costate > -0.1186 or costate < -1.5686:
        """ The costate should be a negative number and is the argument to a tangent. 
        A tangent has singularities at zero and pi/2 (1.57). 
        """
        raise KeyError('Costate value {0} is not a negative argument to a tangent.'.format(costate))
    
    isfound = np.where(alf.Lambda <= costate)
    """ This algorithm is similar to that used at line 290 in GenerateControlTable.py.
    
    Note: np.where returns an ndarray in element [0] with a sequence of indices to values
    which meet the condition.
    Costate values decrease (become more negative) from left to right.  Therefore the 
    first element which meets the condition is the index of the first number more negative than
    the given costate.  Remaining elements all point to more negative values than the given costate.
    The syntax for extracting the first element: "isfound[0][0]".
    """
    
    if np.size(isfound[0], 0) > 0:
        """ Possibly the given costate does not exist in Lambda. Obvious errors are 
        a number between zero and -0.1186, or a number less than -1.5686.
        """                
        found_index = isfound[0][0]
        if found_index == 0:
            raise KeyError('The Key {0} is out of bounds for the Control Table.'.format(costate))
            
        l_found = alf.Lambda[found_index]
        l_before = alf.Lambda[found_index - 1]

        if costate == l_found:
            """ By chance or design the given costate exactly equals a canonical
            costate value.
            """
            return alf.UbyRbyL[l_found][row]
        else:
            """ Interpolate the value of cv """
            l_lo = l_found
            l_hi = alf.Lambda[found_index - 1]
            u_hi = alf.UbyRbyL[l_found][row]
            u_lo = alf.UbyRbyL[l_before][row]
            
            return alf.lin_interp(l_hi, l_lo, costate, u_lo, u_hi)
    else:
        """ Some other problem with the passed in costate value. """
        raise KeyError('The costate value {0} is not valid.'.format(costate))


def read_controlfile(ctlfil = r'..\userfunctions\python\Controls.json'):
    """ Reads the Controls.json file.  
    One trick is that the default value of the path to the Controls.json file 
    is always used for GMAT.  However, for testing, the file is usually in the local
    directory.  The global file pointer is used to indicate when the file is other than 
    the default path.
    """    
    global _fp
    
    try:
        with open(ctlfil, 'r') as _fp:
            dct = alf.load(_fp)
                            
    except OSError as e:
        raise OSError("Invalid JSON filepath: {0} {1}".format(e.filename, e.strerror))

    except Exception as e:
        raise  FileNotFoundError("Exception reading JSON file: {0}".format(e.__doc__))
    
    for l in alf.UbyRbyL:
        """ Parameter dct is a dictionary of lists, key is a string.
        convert to ndarray with key as float64. Loop is elaborated for debug."""
                    
        try:
            U = np.array(dct[str(l)])
            #print ('np.array(dct[str(l)]) = {0}'.format(U))
            alf.UbyRbyL[l] = U
            
        except Exception as e:
            raise RuntimeError('Exception loading UbyRbyL dictionary: {0} for costate {1}.'.format(e.__doc__), l)        

   
def shadow_arc(beta, P, RMAG, MU = 1, RINIT = 6378.136):
    """ This function computes the shadow arc computation from Vallado,
    p.305.  Vallado uses a cylindrical shadow model.
    
    Parameters: 
        beta, the beta angle from Earth to the sun for any season.
        P, the orbital period.
        RMAG, the magnitude of the orbit radius vector
    Returns:
        The shadow arc in radians
    """
    deg2rad = np.pi/180
    
    beta_rads = beta*deg2rad
    
    shadow_size = np.sqrt(1 - (6378.136/RMAG)**2)
    
    earth_angle = shadow_size/np.cos(beta_rads)
    
    if earth_angle <= 1:
        #sharc = np.arccos(earth_angle)*P/np.pi
        sharc = np.arccos(earth_angle)
    else:
        sharc = 0
    #return np.arccos(earth_angle/np.cos(beta_rads))*P/np.pi
    return sharc


def eclipse_weight(beta, P, RMAG):
    """ This function is after the treatment in Journal of Guidance and Control,
    Vol 34, No 1, p.300, Kluever.
    Function computes the time in shadow from sun Beta angle, RMAG, Period and the 
    radius of the Earth.
    Parameters:
        beta, P and RMAG passed to shadow_arc function
    Returns:
        Proportion of orbit that is in sunlight.
    """
    sharc = shadow_arc(beta, P, RMAG)
    return 1 - sharc/(2*np.pi)


class BadCostate(Exception):
    def __init__(self, message):
        self.message = message
        self.__doc__ = "Bad Costate"
  
    
if __name__ == "__main__":
    """    Test cases     """
    import logging
    import platform
    import getpass
        
    logging.basicConfig(
            filename='./GenControls.log',
#            level=logging.INFO,
            level=logging.DEBUG,
            format='%(asctime)s %(filename)s %(levelname)s:\n%(message)s', datefmt='%d%B%Y_%H:%M:%S')

    logging.info("!!!!!!!!!! Control Table Generation Started !!!!!!!!!!")
    
    host_attr = platform.uname()
    logging.info('User Id: %s\nNetwork Node: %s\nSystem: %s, %s, \nProcessor: %s', \
                 getpass.getuser(), \
                 host_attr.node, \
                 host_attr.system, \
                 host_attr.version, \
                 host_attr.processor)
    
    """ Test Case 0: initialized dictionary structure from file. Fundamental to other test cases."""
    read_controlfile(r'.\Controls.json')
       
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    
    AOL = np.linspace(0, 360, 60)
    
    mpl.rcParams['legend.fontsize'] = 10
    
#    fig, axs = plt.subplots(2,1)

    """ Test Case 1: Plot the values of Yaw angles at extremes."""
#    angles21 = get_yaw_angle (AOL, 1.1, -0.1186)
#    axs[0].set_title('Test Case 1: Yaw Angle at R=1.1, Costate -0.1186')
#    axs[0].set_xlabel('Arg of Latitude')
#    axs[0].set_ylabel('Yaw(radians)')
#    plot21 = axs[0].plot(AOL, angles21)
   
#    angles61 = get_yaw_angle (AOL, 6.13, -0.1186)
#    axs[1].set_title('Test Case 1: Yaw Angle at R=6.13, Costate -0.1186')
#    axs[1].set_xlabel('Arg of Latitude')
#    axs[1].set_ylabel('Yaw(radians)')
#    plot61 = axs[1].plot(AOL, angles61)
 
#    plt.tight_layout()
#    plt.show()
#    plt.close()

    fig, axs = plt.subplots(2,1)

    angles21 = get_yaw_angle (AOL, 1.1, -0.4284)
    axs[0].set_title('Test Case 1: Yaw Angle at R=1.1, Costate -0.4284')
    axs[0].set_xlabel('Arg of Latitude')
    axs[0].set_ylabel('Yaw(radians)')
    plot21 = axs[0].plot(AOL, angles21)
   
    angles61 = get_yaw_angle (AOL, 6.6, -0.4284)
    axs[1].set_title('Test Case 1: Yaw Angle at R=6.6, Costate -0.4284')
    axs[1].set_xlabel('Arg of Latitude')
    axs[1].set_ylabel('Yaw(radians)')
    plot61 = axs[1].plot(AOL, angles61)
    
    plt.tight_layout()
    plt.show()
    plt.close()

    fig, axs = plt.subplots(2,1)

    angles21 = get_yaw_angle (AOL, 1.1, -0.55)
    axs[0].set_title('Test Case 1: Yaw Angle at R=1.1, Costate -0.55')
    axs[0].set_xlabel('Arg of Longitude')
    axs[0].set_ylabel('Yaw(radians)')
    plot21 = axs[0].plot(AOL, angles21)
   
    angles61 = get_yaw_angle (AOL, 6.6, -0.55)
    axs[1].set_title('Test Case 1: Yaw Angle at R=6.6, Costate -0.55')
    axs[1].set_xlabel('Arg of Longitude')
    axs[1].set_ylabel('Yaw(radians)')
    plot61 = axs[1].plot(AOL, angles61)
    
    plt.tight_layout()
    plt.show()
    plt.close()

    fig, axs = plt.subplots(2,1)

    angles21 = get_yaw_angle (AOL, 1.1, -0.6214)
    axs[0].set_title('Test Case 1: Yaw Angle at R=1.1, Costate -0.6214')
    axs[0].set_xlabel('Arg of Latitude')
    axs[0].set_ylabel('Yaw(radians)')
    plot21 = axs[0].plot(AOL, angles21)
   
    angles61 = get_yaw_angle (AOL, 6.6, -0.6214)
    axs[1].set_title('Test Case 1: Yaw Angle at R=6.6, Costate -0.6214')
    axs[1].set_xlabel('Arg of Latitudee')
    axs[1].set_ylabel('Yaw(radians)')
    plot61 = axs[1].plot(AOL, angles61)
    
    plt.tight_layout()
    plt.show()
    plt.close()
    
    fig, axs = plt.subplots(2,1)

    angles21 = get_yaw_angle (AOL, 1.1, -0.9695)
    axs[0].set_title('Test Case 1: Yaw Angle at R=1.1, Costate -0.9692')
    axs[0].set_xlabel('Arg of Latitude')
    axs[0].set_ylabel('Yaw(radians)')
    plot21 = axs[0].plot(AOL, angles21)
   
    angles61 = get_yaw_angle (AOL, 6.6, -0.9695)
    axs[1].set_title('Test Case 1: Yaw Angle at R=6.6, Costate -0.9692')
    axs[1].set_xlabel('Arg of Latitude')
    axs[1].set_ylabel('Yaw(radians)')
    plot61 = axs[1].plot(AOL, angles61)
    
    plt.tight_layout()
    plt.show()
    plt.close()


    """ Interpolate the costate and plot the values of Yaw angles."""
    fig, axs = plt.subplots(2,1)

    angles21 = get_yaw_angle (AOL, 1.1, -0.38)
    axs[0].set_title('Test Case 1: Yaw Angle at R=1.1, Interpolate Costate -0.38')
    axs[0].set_xlabel('Arg of Latitude')
    axs[0].set_ylabel('Yaw(radians)')
    plot21 = axs[0].plot(AOL, angles21)
   
    angles61 = get_yaw_angle (AOL, 6.6, -0.38)
    axs[1].set_title('Test Case 1: Yaw Angle at R=6.6, Interpolate Costate -0.38')
    axs[1].set_xlabel('Arg of Latitude')
    axs[1].set_ylabel('Yaw(radians)')
    plot61 = axs[1].plot(AOL, angles61)
    
    plt.tight_layout()
    plt.show()
    plt.close()

 
    """ All Inclination Change """
    fig, axs = plt.subplots(2,1)
    
    inc_thrust_plus = get_yaw_angle(AOL, 6.09, -0.6345)
    axs[0].set_title('Inclination Control - 6.09 with costate -0.6345')
    axs[0].set_xlabel('Argument of Latitude')
    axs[0].set_ylabel('Yaw(radians)')
    plotplus = axs[0].plot(AOL, inc_thrust_plus)
    
    inc_thrust_comb = get_yaw_angle(AOL, 6.11, -0.6345)
    axs[1].set_title('Inclination Control - 6.11 with costate -0.6345')
    axs[1].set_xlabel('Argument of Latitude')
    axs[1].set_ylabel('Yaw(radians)')
    plotplus = axs[1].plot(AOL, inc_thrust_comb)
    
    plt.tight_layout()
    plt.show()
    plt.close()
    
    """ Test Case 2: Eccentricity Change """
    TA = AOL
    fig, axs = plt.subplots(2,1)
    
    ecc_thrust_plus = get_ecc_control(TA)
    axs[0].set_title('Test Case 2: ECC Control - ECC only')
    axs[0].set_xlabel('True Anomaly')
    axs[0].set_ylabel('Pitch(radians)')
    plotplus = axs[0].plot(AOL, ecc_thrust_plus[2])
    
    ecc_thrust_minus = get_ecc_control(TA, False)
    axs[1].set_title('Test Case 2: ECC Control - Combined')
    axs[1].set_xlabel('True Anomaly')
    axs[1].set_ylabel('Pitch(radians)')
    plotplus = axs[1].plot(AOL, ecc_thrust_minus[2])
    
    plt.tight_layout()
    plt.show()
    plt.close()
    

    """ Test Case 3: shadow_arc() and eclipse_weight() """
    betas = np.linspace(-60, 60, 60)
    sharc = np.zeros(60)
    weights = np.ones(60)
    
    P = 18000
    RMAG = 15000
    
    for n in range(0,60):
        beta = betas[n]
        sharc[n] =   shadow_arc(beta, P, RMAG,)
        weights[n] = eclipse_weight(beta, P, RMAG,)
        
    fig, axs = plt.subplots(2,1)

    """ Plot the values of Kleuver Weights. """
    axs[0].set_title('Shadow Angles at 15000km')
    axs[0].set_xlabel('Beta')
    axs[0].set_ylabel('Shadow Arc')
    plot21 = axs[0].plot(betas, sharc)
   
    axs[1].set_title('Kleuver Weights at 15000km')
    axs[1].set_xlabel('Beta')
    axs[1].set_ylabel('Weights')
    plot61 = axs[1].plot(betas, weights)
    
    plt.tight_layout()
    plt.show()
    plt.close()
      

    """ Test Case 4: Thrust Components per Revolution.  
    This test case uses a logic similar to that incorporated in the 
    GMAT AlfanoXfer script.  
    Reference "Simulating Alfano Trajectory with GMAT", author's report.
    
    State:
        AOL = np.linspace(0, 360, 60)
        AOL_init = 99.88774933204886
        Thrust_init = [1,0,0]
        SMA, variable
        SMA_init = 6838.1366
    """
    
    SMA_init = 6838.1366

    """ canonical gravitational constant """
    MU = 1
    """canonical dsistance unit - Earth radius"""
    DU = 6378.1366
    """ canonical time unit - Solar second"""
    TU = 806.81112382429
    """ non-dimensional DU 1.072121"""
    DUstar = SMA_init/DU
    """canonical TU"""
    TUstar = TU * np.sqrt(np.power(DUstar,3)/MU)
    """Orbit Ratio """
    
    R_init = 1
    R_final = 10
    AOL_init = 0
    
    Thrust=[0, 1.0, 0]
        
    with open('thrustlog.log', 'w+') as log:
        log.write('Test Case 4, Thrust.\n')
        log.write('Using costate = {0}.\n'.format(-0.4284))
        log.write('Columns are: \nAOL, Tangental, Yaw, Pitch\n')
        
        for R in range(R_init, R_final, 1):
            SMA = R * SMA_init
            log.write('\nOrbit SMA: {0}\n'.format(SMA))
                
            
            for theta in AOL:                
                Thrust = get_control_onrev(-0.4284, theta, SMA)
                
                log.write('AOL = {0}, Thrust angles: '.format(theta))
                js.dump(Thrust, log)
                log.write('\n')
    
    print('See file "thrustlog.log" for results of Test Case 4.')
    
  
        

 
