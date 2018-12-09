# -*- coding: utf-8 -*-
"""
Created on Tue May 22 19:25:09 2018

Description:
    Simple call to GenerateControlTable.export_yaw_sf().
    
@author: chelms@socal.rr.com
"""
import GenerateControlTable as gct

xfile = r'.\Controls.xlsx'
jfile = r'.\Controls.json'

if __name__ == "__main__":
    
    gct.export_yaw_sf()
    print('Controls.json file created in current directory.')