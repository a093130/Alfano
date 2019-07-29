# -*- coding: utf-8 -*-
"""
Created on Sun May  6 18:36:48 2018
@version 0.3a



@description:
Reference Wiesel and Alfano, 1983, "Optimal Many-Revolution Orbit Transfer".
Reference Sec. 6.7, Vallado, "Fundamentals of Astrodynamics and Applications," 
4th edition.

This Module contains code to output tables that can be joined to determine
the Alfano control variable, u (aka cv in Vallado).  It's purpose is to define 
the thrust profile for an optimum combined inclination change and orbit raising
maneuver.  Its file products are intended to be used with an astrodynamcs
simulation tool such as GMAT.  The output files should be copied to the Python 
user functions directory of GMAT.

An Excel workbook named YawAngleGenerator.xlsx is created in the current directory 
when the module is executed. The workbook is not used further, but is provided for 
inspection and numerical analysis.

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

@author: Colin Helms, chelms@socal.rr.com


@change Log:
    06 May 2018, baseline version, generates 3D plot using AlfanoLib.
    13 May 2018, added mainline code to generate YawAngleGenerator.xlsx.
    15 May 2018, added functions to read workbook, write ControlScaleFactors.json.
    20 May 2018, changed names of output files to be more general and terse.
    09 Mar 2019, added capabilty to generate full controls in Transfer sheet using macros
    15 Mar 2019, New high precision algorithm, per "Simulating Alfano Trajectory ...r0.4.docx"
    21 May 2019, Factored out lin_interp() to AlfanoLib.  Used in YawAngles.py.
"""
import os
import platform
import getpass
import logging
import traceback
import numpy as np
import xlsxwriter as xl
from xlsxwriter import utility as xlut
from PyQt5.QtWidgets import(QApplication, QFileDialog)

from alfano import AlfanoLib as alf

jfilename = 'Controls.json'
""" This file is named following an interface convention with "YawAngles.py". """

def altitude(a, cv):
    """ Compute phi all at once.
    Note that if computing the costate, the costate function should take the
    reciprocal of alfano.phi().
    """
    k = alf.cmp_ell_int_1st_kind(cv)
    e = alf.cmp_ell_int_2nd_kind(cv)
    
    p = alf.alfano_P(cv, k)
    r = alf.alfano_R(cv, k, e)
    
    r_over_p = np.true_divide(r, p)
    
    with np.errstate(divide='ignore', invalid='ignore'):
        """ Avoid divide by zero warning. Courtesy of Stack Overflow.       
        https://stackoverflow.com/questions/26248654/numpy-return-0-with-divide-by-zero 
        """
        alt = alf.PREC_ORBITR * np.true_divide(r_over_p, 2*a)
        alt[alt == np.inf] = 0
        """ alt[ alt = np.inf( alt )] = 0 means find the positions where alt is infinite  
        and set the not-finite values to 0.
        """
        
        alt = np.nan_to_num(alt)
        """ Invalid sets to NaN """
        
    return alt

def lin_interp(l_hi, l_lo, lamb, u_lo, u_hi):
    """ linear interpolation returning value between u_hi and u_lo proportional to 
    lamb between l_hi and l_lo.
    
    Calling procedure may need to round.
    """
    return u_lo + (u_hi - u_lo) * (l_hi - lamb)/(l_hi - l_lo)

def export_controls(data: dict, jfile = ".\\{0}".format(jfilename)):
    """ Function is used to write out a JSON file containing computed Alfano 
    control values by rows of orbit ratio and columns of lambda costate.
    """
    try:       
        """ Dump data to JSON formatted text file """
        with open(jfile, 'w+') as fp:
            #js.dump(data, fp)
            alf.dump(data, fp)
            """ Use overloaded dump function.  Fix to serialize ndarray. """
            
    except OSError as e:
        logging.error("Unable to write file: %s. %s", e.filename, e.strerror)

    except Exception as e:
        lines = traceback.format_exc().splitlines()
        logging.error("Exception writing JSON file: %s\n%s\n%s", e.__doc__, lines[0], lines[-1])
        
    logging.info("{0} trajectory data file written.".format(jfile))
        
if __name__ == "__main__":
    """ Build and export the control table for simulation and flight. """
    __spec__ = None
    """ Necessry tweak to get Spyder IPython to execute this code. 
    See:
    https://stackoverflow.com/questions/45720153/
    python-multiprocessing-error-attributeerror-module-main-has-no-attribute
    """
                
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
   

    app = QApplication([])
    
    xfile = QFileDialog().getSaveFileName(None, 
                       'Select Excel File to Record Generated Costates.', 
                       os.getenv('USERPROFILE'),
                       filter='Excel files(*.xlsx *.xlsm)')

    #xfile = r'./YawAngleGeneratorV2.xlsx'    
    logging.info('Selected Excel Costate file is %s', xfile[0])
    
    wb = xl.Workbook(xfile[0], {'nan_inf_to_errors': True})
    
    cell_format1 = wb.add_format() 
    cell_format1.set_text_wrap()
    cell_format1.set_align('top')
    cell_format1.set_align('left')
    
    cell_bold = wb.add_format() 
    cell_bold.set_bold()
    cell_bold.set_text_wrap()
    cell_bold.set_align('top')
    cell_bold.set_align('center')

    summarySht = wb.add_worksheet('summary')
    costateSht = wb.add_worksheet('costate')
    trajSht = wb.add_worksheet('cv')
    iaSht = wb.add_worksheet('dida')
    xferSht = wb.add_worksheet('transfers')

    """ Documentation of Controls.xlsx workbook, to be written to the summary sheet.   
    """
    summary_textB1 = ('Reference Wiesel and Alfano,'
                      '"Optimal Many-Revolution Orbit Transfer"')
    summary_textA2 = ('Description: ')
    summary_textB3 = ('This workbook is generated by "GenerateControlTable.py". '
                      'The table include in worksheet "cv" is used to generate trajectories '
                      'of combined inclination and orbit raising for spacecraft flight software.')
    summary_textB4 = ('"ExportControls.py" is reads the "cv" table and exports it ' 
                      'as a JSON file.  "YawAngles.py" reads the JSON file and generates thrust '
                      'yaw commands.') 
    summary_textB5 = ('The orbit ratio is elaborated as a linear series from 1-to-10 ' 
                      'in intervals of 0.01. These values are written to the first column, '
                      'rows 2 to 902 in sheets "costate", and "cv".')
    summary_textB6 = ('The control variable cv is elaborated as a linear series '
                      'in 0.001 steps from 0.1000 to 0.9995. These values are '
                      'written  to the first row, columns 1 - 1471 of sheet "costate".')
    summary_textB7 = ('Sheet "costate" contains values of lambda ordered by rows of orbit ratio '
                      'and columns of cv. Lambda is computed using Phi(cv). '
                      'Successive values of lambda are proportional to the first row values '
                      'by the inverse square root of orbit ratio. '
                      'Only the first row values are independent.')
    summary_textB8 = ('Knowing the orbit ratio, and given a value of lambda, the costate table '
                      'is searched in successive rows for value less than or equal to the given lambda.')
    summary_textB9 = ('Sheet "cv" contains columns of cv by rows of orbit ratio and columns of lambda. '
                      'The columns are formed by searching along the rows of sheet "costate" for '
                      'each value of lambda, and identifying each cv value in the same column as found. '
                      'Sheet "cv" represents the Alfano Inverse Phi function. ')
    summary_textB10 = ('The current resolution of 901 x 1470 has been chosen based upon the '
                       'precision required in lambda to achieve a geosynchronous ' 
                       'orbit ratio 6.13 and 28.5 degree inclination change with less than 0.01 '
                       'error in eccentricity.')
    summary_textB11 = ('Sheet "transfers" contains the values of inclination change by rows of '
                       'orbit ratio and columns of lambda. This sheet maps direction to the columns '
                       'of the "cv" sheet.')
    summary_textB12 = ('Text file of cv by rows of Orbit Ratio by columns of Lambda is written out '
                       'in JSON format at conclusion.')


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
    summarySht.write_string('B11', summary_textB11, cell_format1)
    summarySht.write_string('B12', summary_textB12, cell_format1)
    

    wb.define_name('cv','cv!$A$1:$BDO${0}'.format(alf.nrows+1))
    wb.define_name('costates','costate!$A$1:$BDO${0}'.format(alf.nrows+1))

    costates = np.zeros((alf.nrows, alf.ncols))
    """ This stores the 901 x 1471 element range of lambda. """
                
    costateSht.write_string('A1', 'Lambda = f(R,U)', cell_bold)
    costateSht.write_row('B1', alf.u, cell_bold)
    costateSht.write_column('A2', alf.a, cell_bold)
    
    logging.info('Key costate values are:\n{0}'.format(alf.Lambda))
    
    row=0
    mu = 1
    for R in alf.a:
        """ R is the current orbit ratio """
        col=0
        
        L = np.round((mu/np.sqrt(R)) * alf.Lambda, 4)
        """ Compute each row of lambda as a linear combination of the first row. """

        costateSht.write_row(row + 1, col + 1, L)
        
        costates[row] = L
        
        msg = "Developing costates for orbit ratio {0}".format(R)
        print(msg)
        logging.debug(msg)
        
        for lamb in alf.Lambda:
            """ Collect the set of u that is mapped to this costate, each element
            is an instance of u for a particular orbit ratio.  
            There are two tricks to this algorithm:
            (1) Each row of costates is the Lambda multiplied by a factor np.sqrt(mu/R).
            Those values appearing in the first row are therefore canonical. It is only
            necessary to search for the costate values that appear in the first row.
            (2) Costate values in rows decrease (become more negative) from left
            to right.
            
            Lambda and costate are synonyms, Lambda refers to the solution to the 
            Lagrangian multiplier for the constraint on the boundary value problem.
            Control specialists refer to these as costates.
            """
            isfound = np.where(costates[row] <= lamb)
            
            if np.size(isfound[0], 0) > 0:
                """ Values of lambda do not exist in all rows. """
                
                found_index = isfound[0][0]
                """ np.where returns an array with sequence in elemnt [0]
                in which sequence the index of the location is in element [0]
                """
                l_found = costates[row][found_index]
                
                if lamb == l_found:
                    alf.UbyRbyL[lamb][row] = alf.u[found_index]
                else:
                    l_lo = costates[row][found_index]
                    l_hi = costates[row][found_index - 1]
                    u_hi = alf.u[found_index]
                    u_lo = alf.u[found_index - 1]
                    
                    alf.UbyRbyL[lamb][row] = np.round(alf.lin_interp(l_hi, l_lo, lamb, u_lo, u_hi), 4)
            else:
                logging.warning("Costate Table row %d stops at lambda value = %d.", row, lamb)
                break
                
            col += 1
        row += 1
        
    
    logging.info("Completed Calculation of costates. Rows: %d, Columns: %d. Costate worksheet written. ", row, col)
    
    trajSht.write_string('A1', 'cv = f(R,Lambda)', cell_bold)
    trajSht.write_row('B1', alf.Lambda, cell_bold)
    trajSht.write_column('A2', alf.a, cell_bold)
    
    iaSht.write_string('A1', 'di/da = f(R, U)', cell_bold)
    iaSht.write_row('B1', alf.Lambda, cell_bold)
    iaSht.write_column('A2', alf.a, cell_bold)

    xferSht.write_row('B1', alf.Lambda, cell_bold)
    xferSht.write_column('A2', alf.a, cell_bold)

    col = 0  
    for lamb in alf.UbyRbyL:
        """ Write out the 'cv' worksheet.  This is essentially the Alfano Inverse Phi().
        The costate value and orbit ratio are the domain, u as the range. 
        This function is formed by writing out the UbyRbyL dictionary as rows and columns.
    
        The deprecated excel macro is preserved here:
        '=INDEX(cv!$A$1:$SF$500,MATCH($D2,'orbit R'!$A:$A,1),MATCH(E$1,INDIRECT(CONCAT("costate!",$A2,":",$A2)),-1))'
        
        ExportControls.py will read in then write this sheet out as a JSON file.
        """     
        
        UforL = alf.UbyRbyL[lamb]
        """ UforL an array of u. """
        
        try:
            delta_i = altitude(alf.a, UforL)
            
        except Exception as e:
            lines = traceback.format_exc().splitlines()
            
            logging.error("Exception: %s, %s, %s\n%s\n%s", e.__class__, e.__doc__, e.__cause__, lines[0], lines[-1])
            print("Error calculating altitude change: {0}. Continuing.".format(e.__doc__))
            

        """ From the U by R by L data, calculate the incl change for each value of cv 
        by columns of costate and rows of R.  Write to workbook.
        Equation 23 of "Simulating and Alfano Transfer with GMAT is used. 
        """
        
        trajSht.write_column(1, col + 1, UforL)
        iaSht.write_column(1, col + 1, delta_i)
        
        msg ="Writing U column vector for costate {0}".format(lamb)

        print(msg)
        logging.debug(msg)
           
        col += 1

    col = 0
    for cv in alf.u:
        """ Fill in altitude summation formula """
        row = 0
        for R in alf.a:
            #xferSht.write_formula('=SUM(dida!${0}2:{0}:{1}'.format(row + 1, col + 1))
            start_cell = xlut.xl_rowcol_to_cell(2, col + 1)
            sum_cell = xlut.xl_rowcol_to_cell(row + 1, col + 1)
            xferSht.write_formula(sum_cell,'=SUM(dida!${0}:{1})'.format(start_cell, sum_cell))
            
            row +=1
        col+=1
    
    logging.info("Completed Inverse Phi(). Columns: %d. Trajectory worksheet written.", col)
                                
    jfile = QFileDialog().getSaveFileName(None, 
                       'Select JSON Control File for Output.', 
                       os.getenv('USERPROFILE'),
                       filter='Text files(*.json)')
        
    logging.info('Selected JSON Control file is %s', jfile[0])
    
    export_controls(alf.UbyRbyL, jfile[0])

#    mpl.rcParams['legend.fontsize'] = 10
    
#    fig3d2 = plt.figure()
#    ax = Axes3D(fig3d2)

#    ax.set_xlabel('orbit ratio')
#    ax.set_ylabel('control variable')
#    ax.set_zlabel('costates') 

#    row = 0
#    for R in a:
#        col = 0
#        for cv in u:
#            ax.plot_wireframe(R, cv, costates)
#            col +=1
#        row += 1
    
#    plt.show()
#    plt.close()
    """
    TODO: Recreate the Alfano nomograph, figure 6-4 in Vallado. Plot the costate
    value and delta-v as function of orbit ratio and inclination change.  Requires
    rotating the cv axis into the inclination axis, likewise into the delta-V axis.
    """
    
    print("Cleaning up and writing files.")
    wb.close()
        