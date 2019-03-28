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
    09 Mar 2019, added capabilty to generate full controls in Transfer sheet using macros
    15 Mar 2019, New high precision algorithm, per "Simulating Alfano Trajectory ...r0.4.docx"
"""
import os
import platform
import getpass
import logging
import traceback

import numpy as np
import xlsxwriter as xl
import AlfanoLib as alf

xfile = '.\\YawAngleGeneratorV2.xlsx'
jfile = '.\\Controls.json'
mu = 1

def phi(cv):
    k = alf.cmp_ell_int_1st_kind(cv)
    e = alf.cmp_ell_int_2nd_kind(cv)
    dk = alf.derivative_cmp_ell_int_1st(cv, k, e)
    de = alf.derivative_cmp_ell_int_2nd(cv, k, e)
    p = alf.alfano_P(cv, k)
    dp = alf.alfano_Pprime(cv, k, dk)
    r = alf.alfano_R(cv, k, e)
    dr = alf.alfano_Rprime(cv, r, e, dk, de)
    
    return alf.alfano_phi(r, p, dr, dp)

if __name__ == "__main__":
    """ Build and export the control table for simulation and flight. """
    
    import matplotlib as mpl
    from mpl_toolkits.mplot3d import Axes3D 
    import matplotlib.pyplot as plt
            
    logging.basicConfig(
            filename='./GenControlTable.log',
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

    summarySht = wb.add_worksheet('summary')
#    orbitSht = wb.add_worksheet('orbit_R')
    costateSht = wb.add_worksheet('costate')
#    controlSht = wb.add_worksheet('cv')
#    phiSht = wb.add_worksheet('inv_phi')
    trajSht = wb.add_worksheet('traj')
#    xferSht = wb.add_worksheet('transfers')

    """ Documentation of Controls.xlsx workbook, to be written to the summary sheet.   
    """
    summary_textB1 = ('Reference Wiesel and Alfano,'
                      '"Optimal Many-Revolution Orbit Transfer"')
    summary_textA2 = ('Description: ')
    summary_textB3 = ('This workbook is designed to store a table of values of the '
                      'Alfano yaw control variable (cv) ordered by rows of orbit ratio'
                      'and columns of the Alfano lambda costate. '
                      'The table is generated by GenerateControlTable.py.')
    summary_textB4 = ('ExportControls.py is used to read the table and exports it, ' 
                      'as a JSON file for use by mission simulations and flight software. ') 
    summary_textB5 = ('The orbit ratio is elaborated as a linear series from 1-to-10 ' 
                      'in intervals of 0.01. These values are written to the first column '
                      'in sheet "costate", rows 2 to 902. ')
    summary_textB6 = ('The control variable cv is elaborated as a linear series '
                      'in 0.001 steps from 0.1000 to 0.9995. These values are '
                      'written  to the first row, columns 1 - 1471 of sheet "costate".'
                      'The Alphano Phi function is computed from cv. These values are '
                      'written  to the second row, columns 1 - 1471 of sheet "costate".'
                      'Lambda costates are generated as a function of Phi and orbit ratio. ')
    summary_textB7 = ('The values of lambda are written to the rows and columns of sheet "costate". '
                      'The rows of lambda derive from the first row values by a factor, the '
                      'inverse square root of orbit ratio. '
                      'Thus only the first row values are independent.')
    summary_textB8 = ('Knowing the orbit ratio, and a value of lambda, the costate table '
                      'is searched for a value within 0.0006 of the given lambda '
                      'along the row representing the next increment of orbit ratio '
                      'The tolerance 0.0006 results from the expansion of cv and orbit ratio in
                      'linear arrays.'
                      'The column where the given value of lambda is located identifies the, '
                      'the value of cv in row 1. Interpolation is required in most cases.')
    summary_textB9 = ('The sheet "traj" contains columns of cv ordered by rows of orbit ratio. '
                      'Each column is a set of cv associated with a particular value of lambda. '
                      'The columns are formed by collecting the set of cv during the above search. '
                      'This table represents the Alfano Inverse Phi function. ')
    summary_textB10 = ('The current resolution of 901 x 1470 has been chosen based upon the'
                       'precision required in simulation to achieve a geosynchronous ' 
                       'orbit ratio 6.13 and 28.5 degree inclination change. As previously
                       'stated, the precision works out to be +/-0.0006 in Lambda.')

    """ Write out the documentation """    
    summarySht.set_column('A:A',12)
    summarySht.set_column('B:B',76)
    
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
    

    """ write out the costates as functions of cv and a, where i in a and j in cv """
    nrows=901
    ncols=1470
    """ The resultant costate table: 
    901 rows of orbit ratio from 1.00 to 10.00 by steps of 0.01 and  
    1471 columns of lambda from -0.100 to -1.570 by steps of 0.0001.
    """
    wb.define_name('cv_values','cv!$A$1:$KI$1')
    wb.define_name('orbit_ratios','orbit_R!$A$1:$A${0}'.format(nrows+1))
    wb.define_name('costates','costate!$A$1:$KI${0}'.format(nrows+1))

    #u = np.linspace(0, 1, ncols)
    u = np.round(0.1 * np.linspace(1, 10, ncols, endpoint=False), 4)
    """ This is the 1470 element domain of cv from 0 - 1, excluding singularities. """
    a = np.round(np.linspace(1, 10, nrows), 2)
    """ This is the 901 element domain of orbit ratio. """    
    costate = np.zeros((nrows, ncols))
    """ This stores the 901 x 1471 element range of lambda. """
    
    UbyRbyL = {}
    """ The structure of this dictionary is like {costate:{R:cv}}. """
    
    CvPlot = list()
    RPlot = list()
    LPlot = list()
    """ These lists are for plotting the costate matrix """
    
    """ Precompute several useful arrays and factors. """
    halfpi = np.pi/2
    vec_k = alf.cmp_ell_int_1st_kind(u)
    vec_e = alf.cmp_ell_int_2nd_kind(u)
    vec_dk = alf.derivative_cmp_ell_int_1st(u, vec_k, vec_e)
    vec_de = alf.derivative_cmp_ell_int_2nd(u, vec_k, vec_e)
    vec_p = alf.alfano_P(u, vec_k)
    vec_dp = alf.alfano_Pprime(u, vec_k, vec_dk)
    vec_r = alf.alfano_R(u, vec_k, vec_e)
    vec_dr = alf.alfano_Rprime(u, vec_r, vec_e, vec_dk, vec_de)
    vec_phi = alf.alfano_phi(vec_r, vec_p, vec_dr, vec_dp)
    
    row=1
    for R in a:
        col=1
        
        costateSht.write_number(row, 0, R)
        """ a column vector of orbit ratio """
       
        for cv in u:
            """ Calculate the costate for each combination of cv and orbit ratio """                
            costateSht.write_number(0, col, cv)
            """ a row vector of cv """
                           
            L = np.round(alf.costate(phi(cv), R), 4)
            """ Costate matrix is the range of a mapping from domain cv and R. 
            Must limit the precision of this number because it's dictionary hash must
            match the column number. 
            """
            
            costateSht.write_number(row, col, L)
            """ We save the costate matrix in the workbook for verification. """
           
            if L in UbyRbyL:
                """ Continue building u vector for this costate """
                if R in UbyRbyL[L]:
                    UbyRbyL[L][R] = cv
                    
                else:
                    UbyRbyL[L] = {R:cv}
                    
            else:
                """ New costate value. """
                UbyRbyL[L] = {R:cv}
                
                logging.debug("New Costate Value %s for cv %s at orbit ratio %s.", L, cv, R)
            
            col += 1
        row += 1
    
    logging.info("Completed Calculation of costates. Rows: %d, Columns: %d. Costate worksheet written. ", row, col)
    
    row = 1
    col = 1
    
    for L in UbyRbyL:
        """ The 'traj' worksheet will be written with cv as the range 
        of function with domain consisting of a row vector of lambda values
        and column vector of orbit ratio.  
        This matrix is formed by writing out the UbyRbyL dictionary as rows and columns.
    
        The deprecated excel macro is preserved here:
        '=INDEX(cv!$A$1:$SF$500,MATCH($D2,'orbit R'!$A:$A,1),MATCH(E$1,INDIRECT(CONCAT("costate!",$A2,":",$A2)),-1))'
        
        ExportControls.py will read in then write this sheet out as a JSON file.
        """
        trajSht.write(col, 0, L)
        
        for R in UbyRbyL:
            trajSht.write(0, row, R)
            
            cv = UbyRbyL[L][R]
            trajSht.write(col, row, cv)
            
            row += 1
        col += 1

    logging.info("Completed Inversion of costate function with cv. Rows: %d, Columns: %d. Trajectory worksheet written.", row, col)
    
    """
    TODO: from the U by R by L data, calculate the incl change for each value of cv 
    by columns of costate and rows of R.  Write to workbook.
    TODO: compute the delta-v required for each value of cv by columns of costate and 
    rows of R. Write to workbook.
    TODO: Recreate the Alfano nomograph, figure 6-4 in Vallado. Plot the costate
    value and delta-v as function of orbit ratio and inclination change.  Requires
    rotating the cv axis into the inclination axis, likewise into the delta-V axis.
    """
                   
    """ Do wireframe plot of costate matrix """

#    mpl.rcParams['legend.fontsize'] = 10
    
#    fig3d2 = plt.figure()
#    ax = Axes3D(fig3d2)

#    ax.plot_wireframe(CvPlot, RPlot, LPlot)
    
#    ax.set_xlabel('ctl variable')
#    ax.set_ylabel('orbit ratio')
#    ax.set_zlabel('costates') 
#    plt.show()
#    plt.close()
    
    wb.close()
        