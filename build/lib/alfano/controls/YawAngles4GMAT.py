#! python
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 28 11:11:52 2018
@Version 0.2a

@Description:
    This module performs a lookup on a precomputed table of control law scale factors 
    by orbit ratio.  It is designed to be called from within a GMAT mission sequence.

    The result of the lookup is multiplied by the cosine of the true anomaly and the 
    yaw angle is returned in degrees.
    
    Yaw is defined here as as an angle normal to the orbital plane.
    This can be confusing as Alfano defines this angle as pitch, probably following the
    convention of aircraft wherein pitch is normal to the horizon.

    The computation of the yaw control angle is based on Edelbaum's control law, see
    Wiesel and Alfano, "Optimal Many-Revolution Orbit Transfer" and is implemented 
    in package Alfano.AlfanoLib.py.

    The control law scale factor is read from a JSON formatted file as a 901x1470 dictionary:
    {Lambda_1:[SF_i, SF_i+1, ...], Lambda_2:[SF_j, SF_j+1, ...], ..., Lambda_n:[SF_m, SF_m+1, ...]}
    where the SF are in order of orbit ratio from 1.00 to 10.00.

    Internal Dependencies:
        get_control_onrev() ->
            get_control() ->
                get_yaw_angle() -> 
                    scale_fm_json() ->
                        get_yaw_sf()  

@author: Colin Helms

@copyright Copyright Freelance Rocket Science, 2019

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
    05/06/2019, GMAT can't import AlfanoLib, so just copy all the necessary defs here.  This is
    hopefully temporary just to get work done.
"""
import sys
import numpy as np
import json as js

#from alfano import AlfanoLib as alf

PREC_ORBITR = 0.01
""" Precision in orbit ratio must be the same as GenerateControlTable.py """
PREC_COSTATE = 0.001
""" Precision in lambda must be the same as GenerateControlTable.py """
nrows = int(round(1 + (10 - 1)/PREC_ORBITR, 0))
""" Orbit Ratio varies from 1 - 10 """
ncols = int(round((np.pi/2 - 0.1)/PREC_COSTATE, 0))
""" Lambda, l, varies from 0.1 to pi/2 """
UbyRbyL = {}
""" This dictionary is the main interface to the data stored by GenerateControlTable.py 
The structure of this dictionary is {lambda:array(cv)}, where cv is in order of orbit ratio, R.
"""

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
            orbit_r = SMA/SMA_init
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
    
    cv = get_yaw_sf(orbit_r, costate)
    
    if cv != 1 :
        sf = np.sqrt(1/cv - 1)

        #Not sure how this makes sense, but it is consistent with Alfano and Edelbaum.
        theta = np.arctan(np.cos(AOL)/sf)
        """ Make a pi/2 correction to align maximum yaw pi/2 from nodes. """
    else:
        theta = 0
    
    return theta

_fp = None
""" file pointer in Global scope """
         
def get_yaw_sf(orbit_r, costate) :
    """ Function looks up control variable in JSON file based on costate. 
    returns the denominator of the control function.
    Changed 08Apr2019: complete table of costates is searched.
    """
    global _fp
    
    if _fp == None:
        read_controlfile()
        
    costate = np.round(costate, 4)
    r = np.round(orbit_r, 2)
    
    row = int(round(1 + r/PREC_ORBITR, 0)) - 100

    return UbyRbyL[costate][row]             

def read_controlfile(ctlfil = r'..\userfunctions\python\Controls.json'):
    """ Reads the Controls.json file.  
    One trick is that the default value of the path to the Controls.json file 
    is always used for GMAT.  However, for testing, the file is usually in the local
    directory.  The global file pointer is used to indicate when the file is other than 
    the default path.  Because YawAngles.py is called as a library from GMAT, the file
    cannot be located using a GUI.
    """    
    global _fp
    
    try:
        with open(ctlfil, 'r') as _fp:
            dct = load(_fp)
                            
    except OSError as e:
        print("OSError reading JSON file: {0} {1}".format(e.filename, e.strerror))
        sys.exit(-1)

    except Exception as e:
        print("Exception reading JSON file: {0}".format(e.__doc__))
        sys.exit(-1)
    
    for l in UbyRbyL:
        """ Parameter dct is a dictionary of lists, key is a string.
        convert to ndarray with key as float64. Loop is elaborated for debug."""
                    
        try:
            U = np.array(dct[str(l)])
            #print ('np.array(dct[str(l)]) = {0}'.format(U))
            UbyRbyL[l] = U
            
        except Exception as e:
            print("Exception loading UbyRbyL dictionary: {0}".format(e.__doc__))
            break        
   
class BadCostate(Exception):
    def __init__(self, message):
        self.message = message
        self.__doc__ = "Bad Costate"
