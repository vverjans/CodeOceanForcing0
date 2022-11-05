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

def calcDpAveraged(values,depth,dmin,dmax):
    '''
    Calculates depth-averaged values
    '''
    if np.all(depth<=0): #function must use depth positive downwards
        depth = -1*depth
    if depth[0]>0:
        depth  = np.append(0,depth)
        values = np.append(values[0],values) #approxiumation for value at surface

    # Limit domain to [dmin:dmax] range # 
    dp0  = depth[depth<abs(dmax)]
    vl0  = values[depth<abs(dmax)]
    dp1  = dp0[dp0>abs(dmin)]
    vl1  = vl0[dp0>abs(dmin)]
    # Interpolate values at dmin and dmax #
    if np.any(depth<=abs(dmin)):
        dp1  = np.append(dmin,dp1)
        vl1  = np.append(linearInterpolation(depth,values,[dmin]),vl1)
    if np.any(depth>=abs(dmax)):
        dp1 = np.append(dp1,dmax)
        vl1 = np.append(vl1,linearInterpolation(depth,values,[dmax]))

    deltaz      = dp1[1:]-dp1[0:-1]
    valuesav    = (vl1[1:]+vl1[0:-1])/2 #mean value over interval
    averagedval = np.sum(deltaz*valuesav)/np.sum(deltaz) #values integrated over [dmin:dmax]
    return(averagedval)

    
