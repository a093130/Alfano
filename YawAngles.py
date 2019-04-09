# -*- coding: utf-8 -*-
"""
Created on Sat Apr 28 11:11:52 2018

@Description:
    This module performs a lookup on a precomputed table of control law scale factors 
    by orbit ratio.  It is designed to be called from within a GMAT mission sequence.

    The result of the lookup is multiplied by the cosine of the true anomaly and the 
    yaw angle is returned in degrees.
    
    Yaw is defined here as in Vallado, p.381 as an angle within the plane normal
    to the radial vector.  This can be confusing as Alfano uses an LVLH coordinate
    system (RSW) in which the control angle is defined as pitch.

    The computation of the control law is based on Edelbaum's control law, see
    Wiesel and Alfano, "Optimal Many-Revolution Orbit Transfer" and is mainly
    implemented in module AlfanoLib.py.

    The control law scale factor is read from a JSON formatted file as a 2x901 array:
        [[Ri, SFi], 
         [Ri+1, SFi+1], 
         ..., 
         [Rn, SFn]]
    where Rn is an orbit ratio from 1.00 to 10.00, and the SFn are the control
    law a scale factors.

    Internal Dependencies:
        get_control_onrev() ->
            get_control() ->
                get_yaw_angle() -> 
                    scale_fm_json() ->
                        get_yaw_sf()
                    or
                    scale_fm_xlsx() -> 
                        get_yaw_sf()
    
    Windows dependency:
        scale_fm_xlsx() -> xlwings -> Excel 2007 or later
    This library may be used without Windows and Excel.

@author: Colin Helms, chelms@socal.rr.com

@Change Log:
    05/11/2018, Completely rewritten to use a pre-generated table of yaw values.
    Previous version with Chebyshev polynomial saved as AlfanoChebyshev.py.
    05/18/2018, added test cases and adjusted signatures for GMAT call.
    05/20/2018, changed name of input files, integration with GMAT.
    06/01/2018, added logic for Kluever eclipse weighting, and test case.
    03/06/2019, added costate parameter to automate variable inclinations.
    04/08/2019, Controls.json greatly expanded for all costates.
    04/09/2019, AOL correction in get_yaw_angle() to align maximum angle pi/2 from nodes.
"""

#import sys
import platform
#import base64
import getpass
import logging
import traceback
import numpy as np
import json as js

import AlfanoLib as alf

PREC_ORBITR = delta_a = 0.01
""" Precision in orbit ratio must be the same as GenerateControlTable.py """
PREC_COSTATE = delta_l = 0.001
""" Precision in lambda must be the same as GenerateControlTable.py """
mu = 1
""" Use canonical values for orbit elements.  Solve for time and velocity using mu = 1 """
nrows = int(round(1 + (10 - 1)/delta_a, 0))
""" Orbit Ratio varies from 1 - 10 """
ncols = int(round((np.pi/2 - 0.1)/delta_l, 0))
""" Lambda, l, varies from 0.1 to pi/2 """
halfpi = np.pi/2

def get_control_onrev (costate, AOL, SMA, SMA_init = 6838.1366, more = True):
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
        more: defines the direction of the yaw angle, 
            True increases inclination, False decreases inclination.
        costate: the lambda inclination solution.
    """
    if (costate >= -1.570) and (costate <= -0.100):
        if SMA_init >= 6338.1366:
            orbit_r = SMA/SMA_init
        else:
            raise ValueError("SMA_init {0} is invalid.".format(SMA_init))
    else:
        raise BadCostate("Value {0} is out of range.".format(costate))
        
    if (orbit_r >= 1) and (orbit_r <= 10):
        B = 0.0
        """ For combined orbit-raising and inclination change the binormal thrust angle (pitch) is 0. """
        
        if more:
            V, N = get_control(costate, AOL, orbit_r, 1)
            """ increasing inclination """
        else:
            V, N = get_control(costate, AOL, orbit_r, -1)
            """ decreasing inclination """
        
        return [V, N, B]
    
    else:
        raise ValueError ("Orbit ratio {0} is invalid value.".format(orbit_r))
        


def get_control (costate, AOL, orbit_r, nmore):
    """ Function implements the Edelbaum control law.
    
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
    theta = nmore * get_yaw_angle (AOL, orbit_r, costate)
    
    return [np.cos(theta), np.sin(theta)]

def get_yaw_angle (AOL, orbit_r, costate):
    """ Function implements the Edelbaum control law.  Returns the yaw angle
    in degrees.  This function is good for plots.
        
        Returns:
            Thrust angle in radians.
        
        Parameters:
            AOL: the angular position from the line of nodes, in degrees 
            orbit_r: the current orbit ratio  
    """
    
    AOL = AOL*(np.pi/180)
    """ GMAT provides degrees, np.cos expects radians """
    
    cv = CControls().get_yaw_sf(orbit_r, costate)
    
    if cv != 1 :
        sf = np.sqrt(1/cv - 1)

        #theta = np.arctan(np.cos(AOL)/sf)
        theta = np.arctan(np.sin(AOL)/sf)
        """ Make a pi/2 correction to align maximum yaw pi/2 from nodes. """
    else:
        theta = 0
    
    return theta

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

#def loads(*args, **kwargs):
#    """ Overload loads to use the decoder callback. """
#    kwargs.setdefault('object_hook', cb_json_to_ndarray)    
#    return js.loads(*args, **kwargs)

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

class BadStructure(Exception):
    def __init__(self, message):
        self.message = message
        self.__doc__ = "Shape mismatch between structures."

class BadCostate(Exception):
    def __init__(self, message):
        self.message = message
        self.__doc__ = "Bad Costate"

class CControls():
    """ This class is designed as a Singleton instance.  Manage the JSON file of controls. 
    UbyRbyL, delta_a, delta_l, mu, nrows, ncols, halfpi are shared with GenerateControlTable.py.
    
    TODO: factor out all the shared datastructures and variables into a common module.
    """
    _instances = {}
    
    _fp = None
    """ File pointer """
    _u = np.round(0.1 * np.linspace(1, 10, ncols, endpoint=False), 4)
    """ This is the 1470 element domain of cv from 0 - 1, excluding singularities. """
    _a = np.round(np.linspace(1, 10, nrows), 2)
    """ This is the 901 element domain of orbit ratio. """    

    _vec_k = alf.cmp_ell_int_1st_kind(_u)
    _vec_e = alf.cmp_ell_int_2nd_kind(_u)
    _vec_dk = alf.derivative_cmp_ell_int_1st(_u, _vec_k, _vec_e)
    _vec_de = alf.derivative_cmp_ell_int_2nd(_u, _vec_k, _vec_e)
    _vec_p = alf.alfano_P(_u, _vec_k)
    _vec_dp = alf.alfano_Pprime(_u, _vec_k, _vec_dk)
    _vec_r = alf.alfano_R(_u, _vec_k, _vec_e)
    _vec_dr = alf.alfano_Rprime(_u, _vec_r, _vec_e, _vec_dk, _vec_de)  
    _vec_phi = alf.alfano_phi(_vec_r, _vec_p, _vec_dr, _vec_dp)

    _lcanon = np.round(halfpi * (mu/np.sqrt(1)) * 1/_vec_phi, 4)
    
    _UbyRbyL = {l: np.zeros(nrows) for l in _lcanon}
  
    def __new__(CControls, *args, **kwargs):
        """ Make this class a Singleton """
        if CControls not in CControls._instances:
            CControls._instances[CControls] = super().__new__(CControls)
            
        return CControls._instances[CControls]
             
    @classmethod
    def read_controlfile(cls, ctlfil = r'..\userfunctions\python\Controls.json'):
        
        try:
            with open(ctlfil, 'r') as cls._fp:
                dct = load(cls._fp)
                                
        except OSError as e:
            print("OSError reading JSON file: {0} {1}".format(e.filename, e.strerror))
            #exit
    
        except Exception as e:
            lines = traceback.format_exc().splitlines()
            print("Exception reading JSON file: {0}\n{1}\n{2}".format(e.__doc__, lines[0], lines[-1]))
            #exit
        
        for l in cls._UbyRbyL:
            """ Parameter dct is a dictionary of lists, key is a string.
            convert to ndarray with key as float64. Loop is elaborated for debug."""
                        
            try:
                U = np.array(dct[str(l)])
                #print ('np.array(dct[str(l)]) = {0}'.format(U))
                cls._UbyRbyL[l] = U
                
            except Exception as e:
                lines = traceback.format_exc().splitlines()
                print("Exception loading UbyRbyL dictionary: {0}\n{1}\n{2}".format(e.__doc__, lines[0], lines[-1]))
                break
            
    @classmethod        
    def get_yaw_sf(cls, orbit_r, costate) :
        """ Function looks up control variable in JSON file based on costate. 
        returns the denominator of the control function.
        Changed 08Apr2019: complete table of costates is searched.
        """
        if cls._fp == None:
            cls.read_controlfile()
            
        costate = np.round(costate, 4)
        r = np.round(orbit_r, 2)
        
        row = int(round(1 + r/delta_a, 0)) - 100

        return cls._UbyRbyL[costate][row]             

read_controlfile = CControls().read_controlfile
""" This is equivalent to a static method in C++ """        
get_yaw_sf = CControls().get_yaw_sf
""" This is equivalent to a static method in C++ """
    
if __name__ == "__main__":
    """    Test cases     """
        
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
    
    fig, axs = plt.subplots(2,1)

    """ Test Case 1: Plot the values of Yaw angles useJSON = True """
    angles21 = get_yaw_angle (AOL, 1.1, -0.1186)
    axs[0].set_title('Yaw Angle at R=1.1, Costate -0.1186')
    axs[0].set_xlabel('Arg of Latitude')
    axs[0].set_ylabel('Yaw(radians)')
    plot21 = axs[0].plot(AOL, angles21)
   
    angles61 = get_yaw_angle (AOL, 6.13, -0.1186)
    axs[1].set_title('Yaw Angle at R=6.13, Costate -0.1186')
    axs[1].set_xlabel('Arg of Latitude')
    axs[1].set_ylabel('Yaw(radians)')
    plot61 = axs[1].plot(AOL, angles61)
 
    plt.tight_layout()
    plt.show()
    plt.close()

    fig, axs = plt.subplots(2,1)

    angles21 = get_yaw_angle (AOL, 1.1, -0.4284)
    axs[0].set_title('Yaw Angle at R=1.1, Costate -0.4284')
    axs[0].set_xlabel('Arg of Longitude')
    axs[0].set_ylabel('Yaw(radians)')
    plot21 = axs[0].plot(AOL, angles21)
   
    angles61 = get_yaw_angle (AOL, 6.13, -0.4284)
    axs[1].set_title('Yaw Angle at R=6.13, Costate -0.4284')
    axs[1].set_xlabel('Arg of Longitude')
    axs[1].set_ylabel('Yaw(radians)')
    plot61 = axs[1].plot(AOL, angles61)
    
    plt.tight_layout()
    plt.show()
    plt.close()

    """ Test Case 3: Thrust Components per Revolution.  
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
    SMA_final= 42159.48

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
        log.write('Test Case 3\n')
        log.write('Data structure: [True Anomaly, cosine(yaw), sine(yaw), 0]\n')
        
        for R in range(R_init, R_final, 1):
            SMA = R * SMA_init
            subheading = '\nOrbit SMA: '+ str(SMA)
            log.write(subheading)
            log.write('\n')
            
            for theta in AOL:
                V, N, B = get_control_onrev (-0.3289, theta, SMA, more=False)
                """costate, AOL, SMA, SMA_init = 6838.1366, more = True"""
                js.dump([theta, V, N], log)
                log.write('\n')
    
                V, N, B = get_control_onrev (-0.4284, theta, SMA, more=False)
                """costate, AOL, SMA, SMA_init = 6838.1366, more = True"""
                js.dump([theta, V, N], log)
                log.write('\n')
    
    print('See file "thrustlog.log" for results of Test Case 3.')
    
    """ Test Case 4: shadow_arc() and eclipse_weight() """
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

    """ Test Case 1: Plot the values of Yaw angles useJSON = True """
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
    
        

 
