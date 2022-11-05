#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Depth-averaging TF in given DepthRange at all grid points of given Model
Choose scenario of interest

@author: vincent
"""

import os
import sys
import copy
import csv
import numpy as np
import netCDF4 as nc

from codeFunctions import freezingPoint
from codeFunctions import calcDpAveraged

savingTF         = True
cwd              = os.getcwd()+'/'

SelModel         = 'MIROCES2L'
DepthRange       = [0,500] #depth range of interest, [200:500] is from Slater et al. (2019, 2020)
ShallowThreshold = 100 #bathymetry threshold: if bathymetry is shallower, gridpoint is discarded

DirModNC   = f'{cwd}InputOutput/rawCMIP6files/'
DirSave    = f'{cwd}InputOutput/'

To2015hist                 = False
To2100histssp585           = False
To2100histssp126           = True

if(To2015hist):
    partname = 'hist'
elif(To2100histssp585):
    partname = 'hist2100ssp585'
elif(To2100histssp126):
    partname = 'hist2100ssp126'
    
if(SelModel=='MIROCES2L'):
    dim2d              = True
    if(To2015hist):
        ls_members     = [f'r{id}' for id in range(1,30+1)]
    elif(To2100histssp585 or To2100histssp126):
        ls_members     = [f'r{id}' for id in range(1,10+1)]
    elif(To2300histssp585ssp534over):
        ls_members     = ['r1']
elif(SelModel=='IPSLCM6A'):
    dim2d              = True
    if(To2015hist):
        ls_members     = [f'r{id}' for id in range(1,32+1)]
        ls_members.remove('r2') #no r2 member for IPSLCM6A
    elif(To2100histssp585 or To2100histssp126):
        ls_members     = ['r1','r3','r4','r6','r14']
    elif(To2300histssp585ssp534over):
        ls_members     = ['r1']

### Load coordinates from first member ###
file0         = f'thetao_tf_{SelModel}{partname}_{ls_members[0]}.nc' 
ds            = nc.Dataset(DirModNC+file0)
latsfull      = np.array(ds.variables['lat'])
lonsfull      = np.array(ds.variables['lon'])
timefull      = np.array(ds.variables['time'])
depthfull     = np.array(ds.variables['depth'])
ds.close()
adjdepthfull  = np.sort(np.append(DepthRange,depthfull))

if(DepthRange[0]<depthfull[0]):
    iz0 = 0
else:
    iz0 = np.where(depthfull<=DepthRange[0])[0][0]
iz1  = np.where(depthfull>=DepthRange[1])[0][0]
izth = np.where(depthfull>=ShallowThreshold)[0][0]

if(len(np.shape(latsfull))==2):
    dim2 = True
else:
    dim2 = False
if(dim2):
    nny = np.shape(latsfull)[0]
    nnx = np.shape(latsfull)[1]
else:
    nny = len(latsfull)
    nnx = len(lonsfull)


### Depth indices ###
dpmin = DepthRange[0]
### Compute 1 member at a time ###
for mm,member in enumerate(ls_members):
    print(f'Member: {member}')
    tfdpavg    = np.zeros((len(timefull),nny,nnx))
    memberfile = f'thetao_tf_{SelModel}{partname}_{member}.nc'
    ds         = nc.Dataset(DirModNC+memberfile)
    tfprofilefull = np.array(ds.variables['thermalforcing'])
    ds.close()
    for indy in range(nny):
        print(f'indy: {indy}')
        for indx in range(nnx):
            # Extract entire tf profile #
            tf0   = tfprofilefull[:,:,indy,indx]
            # Check if bathy is at least deeper than ShallowThreshold #
            if(tf0[0,izth]<1e10):
                # Find depth index of bathymetry #
                izmax = np.where(tf0[0,:]<1e10)[0][-1]
                # Constrain max depth by DepthRange or bathymetry #
                dpmax = min(DepthRange[1],depthfull[izmax])
                # Calculate average TF over depth range #
                for tt in range(len(timefull)):
                    tfdpavg[tt,indy,indx] = calcDpAveraged(tf0[tt,:],depthfull,dmin=dpmin,dmax=dpmax)
            else:
                # Bathymetry does not go deep enough #
                tfdpavg[:,indy,indx] = 1.1e20
            
    if(savingTF):
        nameout = f'ensemble{SelModel}_{partname}_M{member}_TFdpavg_Dp{DepthRange[0]}to{DepthRange[1]}_bathymin{ShallowThreshold}.nc'
        ### Open netcdf ###
        outnc        = nc.Dataset(DirSave+nameout,'w',format='NETCDF4')
        timedim      = outnc.createDimension('timeDim',size=len(timefull)) 
        zdim         = outnc.createDimension('depthDim',size=len(depthfull)) 
        latdim       = outnc.createDimension('latDim',nny) 
        londim       = outnc.createDimension('lonDim',nnx) 
        
        time_nc      = outnc.createVariable('time','f4',('timeDim',))
        depth_nc     = outnc.createVariable('depth','f4',('depthDim',))
        if(dim2d==True):
            lat_nc   = outnc.createVariable('lat','f4',('latDim','lonDim',))
            lon_nc   = outnc.createVariable('lon','f4',('latDim','lonDim',))
        elif(dim2d==False):
            lat_nc   = outnc.createVariable('lat','f4',('latDim',))
            lon_nc   = outnc.createVariable('lon','f4',('lonDim',))
        tfdpavg_nc   = outnc.createVariable(f'tfdpavg{DepthRange[0]}to{DepthRange[1]}_bathymin{ShallowThreshold}','f4',('timeDim','latDim','lonDim',))
            
        time_nc[:]          = timefull
        depth_nc[:]         = depthfull
        if(dim2d==True):
            lat_nc[:,:]     = latsfull
            lon_nc[:,:]     = lonsfull
        elif(dim2d==False):
            lat_nc[:]       = latsfull
            lon_nc[:]       = lonsfull
        tfdpavg_nc[:,:,:]   = tfdpavg
        
        depth_nc.units     = 'meter'
        time_nc.units      = 'yr'
        tfdpavg_nc.units   = 'degC'
        outnc.close()

            
     
print('End of python job')
os._exit(0) 
    




