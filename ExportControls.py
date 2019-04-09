# -*- coding: utf-8 -*-
"""
Created on Tue May 22 19:25:09 2018

@description:
    Generate JSON control tables from the YawAngleGenerator.xlsx workbook.

@changes
    5/11/2018 - Baseline
    4/3/2019 - modified for Version 2 Excel file and Qt file dialog.
    
@author: colinhelms@outlook.com
"""
import os
import platform
import getpass
import logging
import traceback
import json as js
import xlwings as wings
from PyQt5.QtWidgets import(QApplication, QFileDialog)

def export_controls(xfile = '.\\YawAngleGeneratorV2.xlsx', jfile = '.\\Controls.json'):
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
    
    sht = wb.sheets('cv')

    data = sht.range('cv').options(ndim=2).value
    
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
    
      
