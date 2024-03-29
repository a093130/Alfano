# -*- coding: utf-8 -*-
"""
@name: alfano 
@version: 0.5b1

Created on Sun Apr 28 15:28:56 2019

@author: Colin Helms
@author_email: colinhelms@outlook.com

"""
from __future__ import absolute_import

import sys

if sys.version_info[:2] < (3, 4):
    m = "Python 3.4 or later is required for Alfano (%d.%d detected)."
    raise ImportError(m % sys.version_info[:2])
del sys

from alfano import controls
from alfano import utilities
