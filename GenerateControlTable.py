# -*- coding: utf-8 -*-
"""
Created on Sun May  6 18:36:48 2018
Reference Wiesel and Alfano, 1983, "Optimal Many-Revolution Orbit Transfer".
Reference Sec. 6.7, Vallado, "Fundamentals of Astrodynamics and Applications," 
4th edition.

Description:
This Module contains code to output tables that can be joined to determine
the Alfano control variable, u (aka cv in Vallado).  It's purpose is to define 
the thrust profile for an optimum combined inclination change and orbit raising
maneuver.  Its file products are intended to be used with an astrodynamcs
simulation tool such as GMAT.  The output files should be copied to the Python 
user functions directory of GMAT.

An Excel workbook named YawAngleGenerator.xlsx is created in the current directory 
when the module is executed. 

BEWARE - the program will overwrite YawAngleGenerator.xlsx in the current directory.
  
Functions from AlfanoLib.py are used to compute the analytical values of lambda
given linear arrays of cv and orbit ratio.  The goal is to find the value of cv 
corresponding to the written values of lambda and orbit ratio.  This is the 
inverse function of Phi() in Alfano. The workbook is designed to perform a "Join" 
of the current orbit ratio, cv, and the computed lambda costate using the Excel 
Match functions.  

The optimum value of the lambda costate must be input by the mission planner
in order to obtain the feasible values for his mission. 

A default value of -0.54 is written to the workbook, corresponding to a 
transfer from 500 km to 35786 km. Graphical means of determining the optimum 
lambda value are given in Vallado, Figure 6-24.

Given the lambda costate, and Join, the workbook generates values of the cv 
scale factor by orbit ratio in sheet "trajectory'.  As of Version 1.00 these values 
are generated for orbit ratios from 1-to-10 in increments of 0.01, yielding
901 rows, written in array B2:B902.  The scale factor is used to compute the
vehicle yaw thrust angle, which is periodic with True Anomaly (TA).

Function export_yaw_sf() is provided to read the first two columns of 
sheet "trajectory" and write it out as a JSON file.  It is intended that 
the astrodynamics simulation call a JSON function to read the appropriate 
element from this file for each revolution.

Error analysis of this algorithm shows that for the 500x500 table size,
the max error in lambda is 0.004 and in orbit ratio 0.01.  The error
is computed on inv_phi as the difference in given value and table value.  

Note that the Match function introduces a systematic error since it locates the 
value that is less than or equal to the given value.

Author: Colin Helms, chelms@socal.rr.com

Version 1.03

Change Log:
    06 May 2018, baseline version, generates 3D plot using AlfanoLib.
    13 May 2018, added mainline code to generate YawAngleGenerator.xlsx.
    15 May 2018, added functions to read workbook, write ControlScaleFactors.json.
    20 May 2018, changed names of output files to be more general and terse.
"""
import numpy as np
import json as js
import xlsxwriter as xl
import xlwings as wings
import AlfanoLib as alf

xfile = '.\\Controls.xlsx'
jfile = '.\\Controls.json'

def export_yaw_sf():
    """ Function is used to write out a JSON file containing a 2x901 array
    of computed Alfano yaw control scale factors by orbit ratio.
    
    The structure of the json file is,     
    [[Ri, Yi], 
    [R1+1, Yi+1], 
    ..., 
    [Rn, Yn]]
    where Ri is the current orbit ratio and Yi is the corresponding value of yaw 
    angle at cos(TA) = 1.
    """

    wb = wings.Book(xfile)
    sht = wb.sheets('trajectory')
    #data = sht.range('A2:B902').options(np.array, ndim=2).value
    #Unfortunatesly json does not handle np.array
    data = sht.range('A2:B902').value
    
    """ Dump data to JSON formatted text file """
    with open(jfile, 'w+') as fp:
        js.dump(data,fp)
    
    return data
    
if __name__ == "__main__":
    """
    Test case for AlfanoLib
    """    
    import matplotlib as mpl
    from mpl_toolkits.mplot3d import Axes3D 
    import matplotlib.pyplot as plt
        
    size=500
    u = np.linspace(0, 1, size)
    a = np.linspace(1, 10, size)
    
    wb = xl.Workbook(xfile, {'nan_inf_to_errors': True})
    cell_format1 = wb.add_format() 
    cell_format1.set_text_wrap()
    cell_format1.set_align('top')
    cell_format1.set_align('left')
    
    cell_bold = wb.add_format() 
    cell_bold.bold
    cell_bold.set_text_wrap()
    cell_bold.set_align('top')
    cell_bold.set_align('center')

    """ Documentation of Controls.xlsx workbook, to be written to the 
    summary sheet.  
    """
    summary_textB1 = ('Calculations per Wiesel and Alfano,'
                      '1983, "Optimal Many-Revolution Orbit Transfer"')
    summary_textA2 = ('Description:')
    summary_textB3 = ('This workbook is designed to permit interpolation '
                      'of the Alfano yaw control variable from numerical '
                      'computations for the inclination lambda value (costate) '
                      'from the orbit ratio (Orbit R) and the control variable (cv).')
    summary_textB4 = ('Due to the form of the differential equations of motion, ' 
                      'which contain elliptic integrals and their derivatives, ' 
                      'it is not possible to directly solve for the control variable. ' 
                      'However, it is possible to back-compute the lambda values. '
                      'from a knowledge of cv and the orbit ratio.')
    summary_textB5 = ('The orbit ratio has been written to rows in sheet '
                      '"orbit R" as a linear series from 1-to-10 in intervals '
                      'of 0.01.')
    summary_textB6 = ('Values of the control variable have been written to '
                      'columns in sheet "cv" as a linear series equal in size '
                      'to orbit ratio.  The cv is computed with seven digit precision.')
    summary_textB7 = ('The values of lambda computed from the series of orbit R and cv '
                      'are written to the rows and columns of sheet "costate" in '
                      'cells corresponding to the domain values of orbit ratio and cv.')
    summary_textB8 = ('Knowing lambda, and the orbit ratio, the value of cv '
                      'that corresponds can be found with a Join of the three, '
                      'tables.  This is done on sheet "trajectory" using Excel '
                      'Match macros to find the appropriate row and column from '
                      'sheets "orbit R" and "costate".')
    summary_textB9 = ('Obviously there is error in the calculations, '
                      ' especially since Match only works by finding the column '
                      'of the value of costate less than or equal to the given value. '
                      'This error can be decreased by increasing the resolution '
                      'of the tables.')
    summary_textB10 = ('The current resolution of 500 x 500 has been chosen because '
                       ' there are roughly 500 revolutions in the astrodynamics '
                       'simulation to achieve the geosynchronous orbit ratio (6.13).')
    summary_textB11 = ('The first two columns of sheet "Trajectory" are computed '
                       'in order to be written out to the astrodynamics computation. '
                       'The value of max yaw is used to adjust the values of the '
                       'yaw thrust angle for each revolution, where the orbit ratio '
                       'increases for each revolution.')
    summary_textB12 = ('It is expected that a 500x500 table of costate has adequate '
                       'precision for the astrodynamics computation since the thruster '
                       'direction is only accurate to a few tens of arc seconds and '
                       'the precision of the yaw angle in the table is 36 arcseconds.')
        
    summarySht = wb.add_worksheet('summary')
    orbitSht = wb.add_worksheet('orbit_R')
    costateSht = wb.add_worksheet('costate')
    controlSht = wb.add_worksheet('cv')
    phiSht = wb.add_worksheet('inv_phi')
    trajSht = wb.add_worksheet('trajectory')
    
    """ Write out the documentation """    
    summarySht.write_string('B1', summary_textB1, cell_format1)
    summarySht.write_string('A2', summary_textA2, cell_format1)
    summarySht.write_string('B3', summary_textB3, cell_format1)
    summarySht.write_string('B4', summary_textB4, cell_format1)
    summarySht.write_string('B5', summary_textB5, cell_format1)
    summarySht.write_string('B6', summary_textB6, cell_format1)
    summarySht.write_string('B7', summary_textB7, cell_format1)
    summarySht.write_string('B8', summary_textB8, cell_format1)
    summarySht.write_string('B9', summary_textB9, cell_format1)
    summarySht.write_string('B10', summary_textB10, cell_format1)
    summarySht.write_string('B11', summary_textB11, cell_format1)
    summarySht.write_string('B12', summary_textB12, cell_format1)
  
    summarySht.set_column('A:A',12)
    summarySht.set_column('B:B',76)
  
    """ write out the costates as functions of cv and a, where i in a and j in cv """
    L = np.zeros((size,size))
    X = np.zeros((size,size))
    Y = np.zeros((size,size))
    
    inx=0
    for sma in a:
        iny=0
        for cv in u:
            
            X[inx,iny] = cv
            Y[inx,iny] = sma
            
            orbitSht.write_number(inx, iny, sma)
            controlSht.write_number(inx, iny, cv)
               
            k = alf.cmp_ell_int_1st_kind(cv)
            e = alf.cmp_ell_int_2nd_kind(cv)
            dk = alf.derivative_cmp_ell_int_1st(cv, k, e)
            de = alf.derivative_cmp_ell_int_2nd(cv, k, e)
            p = alf.alfano_P(cv, k)
            dp = alf.alfano_Pprime(cv, k, dk)
            r = alf.alfano_R(cv, k, e)
            dr = alf.alfano_Rprime(cv, r, e, dk, de)
            f = alf.alfano_phi(r, p, dr, dp)
            
            L[inx,iny] = alf.costate(f, sma)
            
            costateSht.write_number(inx, iny, L[inx,iny])
            
            iny += 1
        inx += 1
    
    
    """ Write the phiSht  """
    phiSht.write_number('B1', -0.54)
    phiSht.write_string('A1','Lambda')
    wb.define_name('Lambda','=inv_phi!$B$1:$B$1')

    phiSht.write_number('B2', 6.13)
    phiSht.write_string('A2','Orbit R')
    wb.define_name('Orbit_R', '=inv_phi!$B$2:$B$2')
    
    # Maintenance note: change $A$500 and $SF$500 for greater precision. 
    # Orbit ratio and costate ranges must match.
    wb.define_name('cv_values','cv!$A$1:$SF$500')
    wb.define_name('orbit_ratios','orbit_R!$A$1:$A$500')
    wb.define_name('scale_factors','trajectory!$A$2:$B$902')

    """ Write the trajSht """ 
    """ Orbit Ratios in Column A, Col 0"""
    trajSht.write_string('A1', 'Orbit R')
    R_values = np.linspace(1,10,901)
    trajSht.write_column('A2',R_values)
    
    """ 
    XLXWriter does not support editing capability, fill-down,
    so we have to generate the formulas with correct row-column
    addresses manually.
    """
    trajSht.write_string('B1', 'Scale')
    trajSht.write_string('C1', 'CV')
    trajSht.write_string('D1', 'Row')
    trajSht.write_string('E1', 'Row Ref')
    trajSht.write_string('F1', 'Column')

    for row in range(1, 902, 1):
        """Formula in Column B, Col 1: SQRT(1/C2 - 1)"""
        formulaB = '=SQRT(1/C'
        formulaB += str(row+1)
        formulaB += ' - 1)'
        trajSht.write_formula(row,1,formulaB)
    
        """Formula in Column C, Col 2: INDEX(cv_values,D2,D2)"""
        formulaC = '=INDEX(cv_values, D'
        formulaC += str(row+1)
        formulaC += ',F'
        formulaC += str(row+1)
        formulaC += ')'
        trajSht.write_formula(row,2,formulaC)
    
        """Formula in Column D, Col 3: MATCH(A2,orbit_ratios,1)"""
        formulaD = '=MATCH(A'
        formulaD += str(row+1)
        formulaD += ',orbit_ratios,1)'
        trajSht.write_formula(row,3,formulaD)

        """ Formula in Column E, Col 4: CONCAT("costate!",D2,":",D2) """
        cell_ref = 'D'
        cell_ref += str(row+1)
        formulaE = '=CONCATENATE("costate!",'
        formulaE += cell_ref
        formulaE += ',":",'
        formulaE += cell_ref
        formulaE += ')'
        trajSht.write_formula(row,4,formulaE)
        
        """Formula in Column F, Col 5: MATCH(Lambda,INDIRECT(),-1)"""
        formulaF = '=MATCH(Lambda,INDIRECT(E'
        formulaF += str(row+1)
        formulaF += '),-1)'
        trajSht.write_formula(row,5,formulaF)
    
    trajSht.set_column('E:E',12)
    trajSht.activate()    
    """ Do wireframe Plot of costate """
    
    mpl.rcParams['legend.fontsize'] = 10
    
    fig3d2 = plt.figure()
    ax = Axes3D(fig3d2)
    ax.plot_wireframe(X, Y, L)
    
    ax.set_xlabel('ctl variable')
    ax.set_ylabel('orbit ratio')
    ax.set_zlabel('costates') 
    plt.show()
    plt.close()
    
    wb.close()
        