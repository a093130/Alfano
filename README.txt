The Alfano package is used with the Goddard Mission Analysis ToolKit (GMAT) to implement continuous low-thrust orbit transfer trajectory simulations.  The primary module for use with GMAT is the YawAngles.py procedure found in the Alfano/controls subpackage.  YawAngles.py depends on the AlfanoLib.py from the utilities package. AlfanoLib.py implements the Alfano optimal control equations.  Rather than compute the trajectory functions in real time, GenerateControlTable.py is provided in the controls subpackage and uses AlfanoLib.py to generate a JSON file containing thrust angles by orbit ratio.  YawAngles.py reads the created JSON file to iteratively provide the optimal low-thrust yaw angles to GMAT for each revolution of a circle-to-circle orbit transfer.

