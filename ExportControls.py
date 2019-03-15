# -*- coding: utf-8 -*-
"""
Created on Tue May 22 19:25:09 2018

@description:
    Generate JSON control tables from the YawAngleGenerator.xlsx workbook.
    
@author: colinhelms@outlook.come
"""
import os
import platform
import getpass
import logging
import traceback
import json as js
import xlwings as wings
from PyQt5.QtWidgets import(QApplication, QFileDialog)

def export_yaw_sf(xfile ='.\\YawAngleGenerator.xlsx', jfile = '.\\Controls-costate.json'):
    """ [Deprecated] Function is used to write out a JSON file containing a 2x901 array
    of computed Alfano yaw control scale factors by orbit ratio.
    
    The user must enter the desired value of costate into the Inv Phi tab of the 
    YawAngleGenerator.xlsx workbook.
    
    The structure of the json file is,     
    [
    [Ri, Yi], [Ri+1, Yi+1],...,[Rn, Yn]
    ]
    where Ri is the current orbit ratio and Yi is the corresponding value of the control
    variable.
    """

    wb = wings.Book(xfile)
    sht = wb.sheets('traj')
    """ Unfortunatesly json does not handle np.array """
    data = sht.range('A2:B902').value
    
    """ Dump data to JSON formatted text file """
    with open(jfile, 'w+') as fp:
        js.dump(data,fp)
    
    return data
 
def export_controls(xfile = '.\\YawAngleGenerator.xlsx', jfile = '.\\Controls.json'):
    """ Function is used to write out a JSON file containing a 296x901 array
    of computed Alfano yaw control scale factors by orbit ratio for values of lambda
    costate from -0.1 to 1.57 in increments of 0.005.
    
    The structure of the json file is,     
    [
    Lj[[Ri, Yi], [Ri+1, Yi+1],...,[Rn, Yn]], 
    Lj+1[[Ri, Yi], [Ri+1, Yi+1],...,[Rn, Yn]],
    ...,
    Lj+m[[Ri, Yi], [Ri+1, Yi+1],...,[Rn, Yn]]
    ]
    where Lj is the costate, Ri is the current orbit ratio, and 
    Yi is the corresponding value of the control variable.
    """

    wb = wings.Book(xfile)
    sht = wb.sheets('Transfers')
    #data = sht.range('A2:B902').options(np.array, ndim=2).value
    #Unfortunatesly json does not handle np.array
    data = sht.range('Controls').options(ndim=2).value
    
    """ Dump data to JSON formatted text file """
    with open(jfile, 'w+') as fp:
        js.dump(data,fp)
    
    return data

if __name__ == "__main__":
    logging.basicConfig(
            filename='./XportCtls.log',
            level=logging.INFO,
            format='%(asctime)s %(filename)s %(levelname)s:\n%(message)s', datefmt='%d%B%Y_%H:%M:%S')

    logging.info("!!!!!!!!!! JSON File Generation Started !!!!!!!!!!")
    
    host_attr = platform.uname()
    logging.info('User Id: %s\nNetwork Node: %s\nSystem: %s, %s, \nProcessor: %s', \
                 getpass.getuser(), \
                 host_attr.node, \
                 host_attr.system, \
                 host_attr.version, \
                 host_attr.processor)
    
    app = QApplication([])
    
    fname = QFileDialog().getOpenFileName(None, 'Open YawAngleGenerator.xlsx.', 
                       os.getenv('USERPROFILE'))
        
    logging.info('Path to YawAngleGenerator.xlsx is %s', fname[0])
    
  
    try:
        export_controls(fname[0])
    
        print('Created Controls.json in the current directory.')
        logging.info('Created Controls.json in the current directory.')

    except Exception as e:
        lines = traceback.format_exc().splitlines()
        logging.error("Exception: %s, %s, %s\n%s\n%s", e.__class__, e.__doc__, e.__cause__, lines[0], lines[-1])
    
      
