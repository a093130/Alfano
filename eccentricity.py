# name: eccentricity
# description: various functions for conic sections dependent on eccentricity
import math

def scale(ecc: int = 1): # create the scaling factor
    """ Scaling Factor

    Function, scale() relates the semi-major axis to the semi-minor axis.
    Arguments:
    ecc -- eccentricity, always a positive integer
    """
    return math.sqrt(1 - ecc**2)
