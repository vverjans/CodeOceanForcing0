#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functions 
@author: vincent
"""

import os
import sys
import copy
import numpy as np
import netCDF4 as nc
import csv
from scipy import interpolate

def freezingPoint(sal,depth):
    '''Freezing point temperature of seawater from Cowton et al. (2015)'''
    depth = abs(depth) #make sure to use depth values increasing downwards
    lb1 = -5.73e-2 #[°C psu-1] Cowton et al. (2015) Table 1
    lb2 = 8.32e-2 #[°C] Cowton et al. (2015) Table 1
    lb3 = -7.61e-4 #[°C m-1] Cowton et al. (2015) Table 1 (negative -> more negative as depth increases downwards)
    tempfr = lb1*sal+lb2+lb3*depth #Cowton et al. (2015) Eq.7
    return(tempfr)



    
