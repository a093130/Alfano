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
"""
import numpy as np
import json as js
#import xlwings as wings

""" The following path is specific to use within GMAT """
#jfile = r'..\userfunctions\python\Controls.json'

""" Use this alternatepath for debug within the development environment """
jfile = r'.\Controls.json'

""" GMAT is having trouble locating the win32api module
Windows features must be disabled for use in GMAT """
#xfile = r'.\Controls.xlsx'


#def get_control_onrev (TA, TA_init, SMA, SMA_init=6338.1366, Control=[6878.136,1,0,0]):
#    """ Function provides a wrapper to perform conversions from SMA in km
#    to the orbit ratio; the given SMA is divided by SMA_init.  SMA_init is
#    set to the canonical Earth radius to avoid divide by zero, but the 
#    caller should pass in the SMA of the starting orbit to get correct results.
#    
#    Function checks for completed revolution and calls get_yaw_angle() with 
#    given TA and computed orbit_r, sets TA_init equal to the current TA, else
#    returns original yaw angle.
#    
#    Future: eclipse effects will accounted for in function set_apsides().
#   
#    Returns:
#        [TA_init, Thrust]
#    See below for definition of Thrust and TA.
#        
#    Parameters:
#        TA: TA: the current true anomaly, in degrees
#        TA_init: the true anomaly at start of rev
#        Thrust_init: 3 Thrust unit components in VNB coordinates [T1, T2, T3]
#        SMA: the current SMA in kilometers
#        SMA_init: the SMA of the initial orbit in kilometers
#    """   
#    if SMA_init == 0:
#        return Control
#    
#    if TA > TA_init:
#        orbit_r = SMA/SMA_init
#        TA_init = TA - 0.00000000000001
#        Thrust = get_control(TA, orbit_r, True)
#        Control = [TA_init, Thrust[0], Thrust[1], Thrust[2]]
#        print(Control)
#    """GMAT can only accept one dimensional arrays """
#    return Control

def get_control_onrev (TA, SMA, SMA_init = 6838.1366):
    """ Function provides a wrapper to perform conversions from SMA in km
    to the orbit ratio; the given SMA is divided by SMA_init.  SMA_init is
    set to the canonical Earth radius to avoid divide by zero, but the 
    caller should pass in the SMA of the starting orbit to get correct results.
    
    Function checks for completed revolution and calls get_yaw_angle() with 
    given TA and computed orbit_r, sets TA_init equal to the current TA, else
    returns original yaw angle.
    
    Future: eclipse effects will accounted for in function set_apsides().
    
    Returns:
        [TA_init, Thrust]
    See below for definition of Thrust and TA.
        
    Parameters:
        TA: the current true anomaly, in degrees
        TA_init: the true anomaly at start of rev
        Thrust: Tangential and Normal Thrust unit components in VNB coordinates.
        SMA: the current SMA in kilometers
        SMA_init: the SMA of the initial orbit in kilometers
    """
    
    orbit_r = SMA/SMA_init
    
    V, N = get_control(TA, orbit_r, True)
    B = 0.0
        
    return [V, N, B]

def get_control (TA, orbit_r, useJSON):
    """ Function implements the Edelbaum control law.
    
            yaw_angle = arctan[cos(TA)/sqrt(1/u - 1)].
        
        This function is good for call from GMAT.
        
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
            TA: the true anomaly, in degrees 
            orbit_r: the current orbit ratio
            useJSON: two methods of finding sf are used.
    """
    
    """ The Yaw Angle is defined negative for reducing inclination. """
    theta = -1 * get_yaw_angle (TA, orbit_r, useJSON)
    
    return [np.cos(theta), np.sin(theta)]

def get_yaw_angle (TA, orbit_r, useJSON):
    """ Function implements the Edelbaum control law.  Returns the yaw angle
    in degrees.  This function is good for plots.
        
        Returns:
            Thrust angle in radians.
        
        Parameters:
            TA: the true anomaly, in degrees 
            orbit_r: the current orbit ratio  
    """
    
    """ GMAT provide degrees, np.cos expects radians """
    """ GMAT provides degrees, np.cos expects radians """
    TA = TA*(np.pi/180)
    
    """ Use method selected """
    if useJSON:
        sf = scale_fm_json(orbit_r)
    else:
#        sf = scale_fm_xlsx(orbit_r)
        print('Not implemented because of GMAT issue with win32api.')
    
    theta = np.arctan(np.cos(TA)/sf)
    
    return theta

#def scale_fm_xlsx(orbit_r):
    """ Function will find or open an instance of an Excel server and 
    retrieve the yaw angle for a given orbit ratio.  The orbit ratio
    is interpolated to 0.01 precision.
    
    Requirements: Windows platform and an instance of Excel 2007 or later.
    Note that in the unlikely occasion that two Excel servers have the
    Controls workbook open, this code will select the first instance.
    
    Note that this code may have better performance than
    """    
#    wb = wings.Book(xfile)
#    sht = wb.sheets('trajectory')
#    data = sht.range('A2:B902').value
    
#    return get_yaw_sf(orbit_r, data)

def scale_fm_json(orbit_r):
    """ Function looks up scale factor in JSON file. """
    
    with open(jfile, 'r+') as fp:
        data = js.load(fp)
        
    return get_yaw_sf(orbit_r, data)
  
def get_yaw_sf(orbit_r, data):
    """ Function is used to parse a 2x901 array
    of computed Alfano yaw control scale factors by rows of orbit ratio.
    
    The structure of the data is, 
    [[Ri, Yi], 
    [R1+1, Yi+1],
    ..., 
    [Rn, Yn]]
    where Ri is the current orbit ratio and Yi is the corresponding value of yaw 
    scale factor at cos(TA) = 1.
    
    This function assumes an external python call to read the JSON file 
    has been made in the outer scope.
    """
    #If the shape of the data array changes, then revise this computation for a.
    a = int(np.round(orbit_r - 1, 2) * 100)
    
    return data[a][1]

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

if __name__ == "__main__":
    """    Test cases     """
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    
    TA = np.linspace(0, 360, 60)
    
    mpl.rcParams['legend.fontsize'] = 10
    
    fig, axs = plt.subplots(2,1)

    """ Test Case 1: Plot the values of Yaw angles useJSON = True """
    angles21 = get_yaw_angle (TA, 1.1, True)
    axs[0].set_title('Yaw Angle at R=1.1, Costate -0.64')
    axs[0].set_xlabel('Arg of Latitude')
    axs[0].set_ylabel('Yaw(radians)')
    plot21 = axs[0].plot(TA, angles21)
   
    angles61 = get_yaw_angle (TA, 6.13, True)
    axs[1].set_title('Yaw Angle at R=6.13, Costate -0.64')
    axs[1].set_xlabel('Arg of Latitude')
    axs[1].set_ylabel('Yaw(radians)')
    plot61 = axs[1].plot(TA, angles61)
    
    plt.tight_layout()
    plt.show()
    plt.close()

    """ Test Case 2: Plot the values of Yaw angles using Excel, useJSON = False """
#    fig, axs = plt.subplots(2,1)
    
#    angles21 = get_yaw_angle (TA, 1.1, False)
#    axs[0].set_title('Excel Test Case: Yaw Angle at R=1.1')
#    axs[0].set_xlabel('TA')
#    axs[0].set_ylabel('Yaw(radians)')
#    plot21 = axs[0].plot(TA, angles21)
   
#    angles61 = get_yaw_angle (TA, 6.1, False)
#    axs[1].set_title('Excel Test Case: Yaw Angle at R=6.1')
#    axs[1].set_xlabel('TA')
#    axs[1].set_ylabel('Yaw(radians)')
#    plot61 = axs[1].plot(TA, angles61)
    
#    plt.tight_layout()
#    plt.show()
#    plt.close()

    """ Test Case 3: Thrust Components per Revolution.  
    This test case uses a logic similar to that incorporated in the 
    GMAT AlfanoXfer script.  
    Reference "Simulating Alfano Trajectory with GMAT", author's report.
    
    State:
        TA = np.linspace(0, 360, 60)
        TA_init = 99.88774933204886
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
    TA_init = 99.8877493
    
    Thrust=[0,1.0,0]
    
    with open('thrustlog.log', 'w+') as log:
        log.write('Test Case 3\n')
        log.write('Data structure: [True Anomaly, cosine(yaw), sine(yaw), 0]\n')
        
        for R in range(R_init, R_final, 1):
            SMA = R * SMA_init
            subheading = '\nOrbit SMA: '+ str(SMA)
            log.write(subheading)
            log.write('\n')
            
            for theta in TA:
                V, N, B = get_control_onrev (theta, SMA, SMA_init)
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
    
        

 
