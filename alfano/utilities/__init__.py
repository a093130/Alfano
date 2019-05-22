# -*- coding: utf-8 -*-
"""
@Name: AlfanoLib 
@Version: 2.0

Created on Sun Apr 28 15:28:56 2019
Deployed to <PYTHONPATH>/alfano

@author: Colin Helms
@author_email: colinhelms@outlook.com

@Description: This package contains library functions derived from 
"Optimal Many-revolution Orbit Transfer," Alfano & Wiesel 1985.

"""
from __future__ import absolute_import
import sys

if sys.version_info[:2] < (3, 4):
    m = "Python 3.4 or later is required for Alfano (%d.%d detected)."
    raise ImportError(m % sys.version_info[:2])
del sys

#from .AlfanoLib import *