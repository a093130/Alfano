#! python
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 28 11:11:52 2018
@Version 0.2a

@Description:
    This module is a simple wrapper to call the library function get_control_onrev. 

@author: Colin Helms

@copyright Copyright Freelance Rocket Science, 2019

"""
from alfano import YawAngles as ctl

def get_control_onrev (costate, AOL, SMA, SMA_init = 6838.1366, more=-1):
    """ Function provides a wrapper to library model YawAngles.
    
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
    return ctl.get_control_onrev(costate, AOL, SMA, SMA_init, more)
        

if __name__ == "__main__":
    """    Test cases     """
    import logging
    import platform
    import getpass
    import numpy as np
    import json as js
        
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
    ctl.read_controlfile(r'.\Controls.json')
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
    AOL = np.linspace(0, 360, 60)
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
                V, N, B = get_control_onrev (-0.3289, theta, SMA, -1)
                """costate, AOL, SMA, SMA_init = 6838.1366, more = True"""
                js.dump([theta, V, N], log)
                log.write('\n')
    
                V, N, B = get_control_onrev (-0.4284, theta, SMA, -1)
                """costate, AOL, SMA, SMA_init = 6838.1366, more = True"""
                js.dump([theta, V, N], log)
                log.write('\n')
    
    print('See file "thrustlog.log" for results of Test Case 3.')
    

