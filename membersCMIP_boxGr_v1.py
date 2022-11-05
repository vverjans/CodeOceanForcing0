#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Boxing in Greenland domain
Computing the 1850-2300 TF and thetao time series of CMIP6 models

@author: vincent
"""

import os
import sys
import copy
import csv
import numpy as np
import netCDF4 as nc

from codeFunctions import freezingPoint

saveBoxGreenlandNC = True
cwd                = os.getcwd()+'/'

SelModel = 'MIROCES2L'

DirThetaoNC = f'{cwd}InputOutput/rawCMIP6files/'
DirSoNC     = f'{cwd}InputOutput/rawCMIP6files/'
DirSaveNC   = f'{cwd}InputOutput/'

### Select experiment ###
To2015hist                 = False
To2100histssp585           = False
To2100histssp126           = True

### Limits of Greenland domain ###
limN           = 86.0
limS           = 57.0
limE           = 4.0
limW           = 274.0
limDp          = 1200.0
depthSubSample = 1

if(To2015hist):
    Experiments = ['hist']
    DatesCut    = [2015]
elif(To2100histssp585): 
    Experiments = ['hist','ssp585']
    DatesCut    = [2015,2100]
elif(To2100histssp126): 
    Experiments = ['hist','ssp126']
    DatesCut    = [2015,2100]
nExp          = len(Experiments)
depthUnitConv = 1.0 #initialize depth unit converter

if(SelModel=='MIROCES2L'):
    dim2d              = True
    if(To2015hist):
        ls_members     = [f'r{id}' for id in range(1,30+1)]
    elif(To2100histssp585 or To2100histssp126):
        ls_members     = [f'r{id}' for id in range(1,10+1)]
    namelat            = 'latitude'
    namelon            = 'longitude'
    namez              = 'lev'
    datesendhist       = np.array(['201412'])
    if(To2100histssp585):
        datesendssp585     = np.array(['210012'])
    if(To2100histssp126):
        datesendssp126     = np.array(['210012'])
        
if(SelModel=='IPSLCM6A'):
    dim2d              = True
    if(To2015hist):
        ls_members     = [f'r{id}' for id in range(1,32+1)]
        ls_members.remove('r2') #no r2 member for IPSLCM6A
    elif(To2100histssp585 or To2100histssp126):
        ls_members     = ['r1','r3','r4','r6','r14']
    namelat            = 'nav_lat'
    namelon            = 'nav_lon'
    namez              = 'olevel'
    datesRef           = [1850.0,2015.0,2040.0] 
    datesendhist       = np.array(['194912','201412'])
    if(To2100histssp585):
        datesendssp585     = np.array(['210012'])
    if(To2100histssp126):
        datesendssp126     = np.array(['210012'])

else:
    print(f'Error script not implemented yet for {SelModel}')

fullfilesthetao = os.listdir(DirThetaoNC)
fullfilesso     = os.listdir(DirSoNC)
nMemb           = len(ls_members)

### Read all files of each experiment for each member ###
if(To2015hist):
    lsFilesThetaExps = [[[]] for mm in range(nMemb)]
    lsFilesSoExps    = [[[]] for mm in range(nMemb)]
    for mm,memb in enumerate(ls_members):
        for dd,dend in enumerate(datesendhist):
            file = [fl for fl in fullfilesthetao if ((fl[0:6]=='thetao') and (Experiments[0] in fl) and (f'_{memb}i' in fl) and (dend in fl))][0]
            lsFilesThetaExps[mm][0].append(file)
            file = [fl for fl in fullfilesso if ((fl[0:2]=='so') and (Experiments[0] in fl) and (f'_{memb}i' in fl) and (dend in fl))][0]
            lsFilesSoExps[mm][0].append(file)
elif(To2100histssp585):
    lsFilesThetaExps = [[[],[]] for mm in range(nMemb)]
    lsFilesSoExps    = [[[],[]] for mm in range(nMemb)]
    for mm,memb in enumerate(ls_members):
        for dd,dend in enumerate(datesendhist):
            file = [fl for fl in fullfilesthetao if ((fl[0:6]=='thetao') and (Experiments[0] in fl) and (f'_{memb}i' in fl) and (dend in fl))][0]
            lsFilesThetaExps[mm][0].append(file)
            file = [fl for fl in fullfilesso if ((fl[0:2]=='so') and (Experiments[0] in fl) and (f'_{memb}i' in fl) and (dend in fl))][0]
            lsFilesSoExps[mm][0].append(file)
        for dd,dend in enumerate(datesendssp585):
            file = [fl for fl in fullfilesthetao if ((fl[0:6]=='thetao') and (Experiments[1] in fl) and (f'_{memb}i' in fl) and (dend in fl))][0]
            lsFilesThetaExps[mm][1].append(file)
            file = [fl for fl in fullfilesso if ((fl[0:2]=='so') and (Experiments[1] in fl) and (f'_{memb}i' in fl) and (dend in fl))][0]
            lsFilesSoExps[mm][1].append(file)
elif(To2100histssp126):
    lsFilesThetaExps = [[[],[]] for mm in range(nMemb)]
    lsFilesSoExps    = [[[],[]] for mm in range(nMemb)]
    for mm,memb in enumerate(ls_members):
        for dd,dend in enumerate(datesendhist):
            file = [fl for fl in fullfilesthetao if ((fl[0:6]=='thetao') and (Experiments[0] in fl) and (f'_{memb}i' in fl) and (dend in fl))][0]
            lsFilesThetaExps[mm][0].append(file)
            file = [fl for fl in fullfilesso if ((fl[0:2]=='so') and (Experiments[0] in fl) and (f'_{memb}i' in fl) and (dend in fl))][0]
            lsFilesSoExps[mm][0].append(file)
        for dd,dend in enumerate(datesendssp126):
            file = [fl for fl in fullfilesthetao if ((fl[0:6]=='thetao') and (Experiments[1] in fl) and (f'_{memb}i' in fl) and (dend in fl))][0]
            lsFilesThetaExps[mm][1].append(file)
            file = [fl for fl in fullfilesso if ((fl[0:2]=='so') and (Experiments[1] in fl) and (f'_{memb}i' in fl) and (dend in fl))][0]
            lsFilesSoExps[mm][1].append(file)

def datedecMIROCES2L(daysar):
    '''
    MIROCES2L time comes in units of days since 01Jan1850
    '''
    output  = 1850+daysar/365.25
    return(output)

def datedecIPSLCM6A(daysar,dateref):
    '''
    IPSLCM6A time comes in units of days since a reference date
    Provide the reference date as dateref
    '''
    output = dateref+daysar/365.25
    return(output)

### Load 1 file for geographical processing ###
path0      = DirThetaoNC+lsFilesThetaExps[0][0][0]
ds = nc.Dataset(path0)
lat1   = np.array(ds.variables[namelat])
lon1   = np.array(ds.variables[namelon])
depth1 = np.array(ds.variables[namez])
ds.close()
# Convert depth units #
depth1 = depthUnitConv*depth1 #units of [m]

### Select regions and prepare spatial outputs ###
if(dim2d):
    indsSel = np.logical_and(np.logical_or(\
               lon1>limW,lon1<limE),
               np.logical_and(\
               lat1>limS,lat1<limN))
    jj1,jj2 = np.where(indsSel)
    latout   = lat1[min(jj1):max(jj1)+1,min(jj2):max(jj2)+1].astype(np.float32)
    lonout   = lon1[min(jj1):max(jj1)+1,min(jj2):max(jj2)+1].astype(np.float32)
    zzMax    = np.where(abs(depth1)<limDp)[0][-1]+1
    inoloadz = np.ones(len(depth1)).astype(bool)
    inoloadz[0:zzMax:depthSubSample] = False
    inoloadz[zzMax]                  = False
    iloadz     = inoloadz==False
    depthout   = depth1[iloadz]
    ndpout     = len(depthout)
    nlatout    = np.shape(latout)[0]
    nlonout    = np.shape(latout)[1]
else:
    latGr = []
    lonGr = []
    for ii,latval in enumerate(lat1):
        if(latval>=limS and latval<=limN):
            for jj,lonval in enumerate(lon1):
                if(lonval<=limE or lonval>=limW):
                    latGr.append(ii)
                    lonGr.append(jj)
    iilatmin  = min(latGr)     
    iilatmax  = max(latGr) 
    zzMax    = np.where(abs(depth1)<limDp)[0][-1]+1
    inoloadz = np.ones(len(depth1)).astype(bool)
    inoloadz[0:zzMax:depthSubSample] = False
    inoloadz[zzMax]                  = False
    iloadz     = inoloadz==False
    depthout   = depth1[iloadz]
    sellon     = np.unique(lonGr)
    lonout     = lon1[sellon]
    latout     = lat1[iilatmin:iilatmax+1]
    ndpout     = len(depthout)
    nlatout    = len(latout)
    nlonout    = len(lonout)
    
for mm,member in enumerate(ls_members):
    if(mm==0):
        # Initialize time array output #
        timingfull = np.array([])
    for kk in range(nExp):
        nFl = len(lsFilesThetaExps[mm][kk])
        for ff in range(nFl):
            # Load theta file and variable #
            locthetao = DirThetaoNC+lsFilesThetaExps[mm][kk][ff]
            ds = nc.Dataset(locthetao)
            fftime = np.array(ds.variables['time'])
            # Convert time to decimal time #
            if(SelModel=='MIROCES2L'):
                fftiming = datedecMIROCES2L(fftime)
            if(SelModel=='IPSLCM6A'):
                fftiming = datedecIPSLCM6A(fftime,datesRef[kk])
            # Find time indices of interest #
            if(kk==0):
                indstime = fftiming<=DatesCut[0]
            else:
                indstime = np.logical_and(fftiming>DatesCut[kk-1],fftiming<=DatesCut[kk])
            if(mm==0):
                timingfull = np.append(timingfull,fftiming[indstime])
            # Load thetao #
            if(dim2d==True):
                theta0 = np.array(ds.variables['thetao'][indstime,iloadz,min(jj1):max(jj1)+1,min(jj2):max(jj2)+1])
            elif(dim2d==False):
                theta0 = np.array(ds.variables['thetao'][indstime,iloadz,iilatmin:iilatmax+1,sellon])
            ds.close()
            if(kk==0 and ff==0):
                thetafull = np.copy(theta0)
            else:
                thetafull = np.row_stack((thetafull,theta0))
            # Load salinity file and variable #
            locso = DirSoNC+lsFilesSoExps[mm][kk][ff]
            ds = nc.Dataset(locso)
            if(dim2d==True):
                salt0 = np.array(ds.variables['so'][indstime,iloadz,min(jj1):max(jj1)+1,min(jj2):max(jj2)+1])
            elif(dim2d==False):
                salt0 = np.array(ds.variables['so'][indstime,iloadz,iilatmin:iilatmax+1,sellon])
            ds.close()
            # Freezing point profile #
            frpoint  = freezingPoint(salt0,depthout.reshape(1,len(depthout),1,1))
            # Thermal Forcing profile #
            fftf     = theta0-frpoint
            # Mask out grounded values #
            imaskout       = theta0>1e10
            fftf[imaskout] = 1.1e20
            if(kk==0 and ff==0):
                tffull = np.copy(fftf)
            else:
                tffull = np.row_stack((tffull,fftf))
            del(theta0,salt0,frpoint,fftf)
        
    if(saveBoxGreenlandNC):
        if(To2015hist):
            nameout = f'thetao_tf_{SelModel}hist_{member}.nc'
        elif(To2100histssp585):
            nameout = f'thetao_tf_{SelModel}hist2100ssp585_{member}.nc'
        elif(To2100histssp126):
            nameout = f'thetao_tf_{SelModel}hist2100ssp126_{member}.nc'
        ### Open netcdf ###
        outnc        = nc.Dataset(DirSaveNC+nameout,'w',format='NETCDF4')
        timedim      = outnc.createDimension('timeDim',size=len(timingfull)) 
        zdim         = outnc.createDimension('depthDim',size=len(depthout))
        if(dim2d):
            latdim   = outnc.createDimension('latDim',np.shape(latout)[0]) 
            londim   = outnc.createDimension('lonDim',np.shape(latout)[1]) 
        elif(dim2d==False):
            latdim   = outnc.createDimension('latDim',len(latout)) 
            londim   = outnc.createDimension('lonDim',len(lonout))
        
        time_nc      = outnc.createVariable('time','f4',('timeDim',))
        depth_nc     = outnc.createVariable('depth','f4',('depthDim',))
        if(dim2d):
            lat_nc   = outnc.createVariable('lat','f4',('latDim','lonDim',))
            lon_nc   = outnc.createVariable('lon','f4',('latDim','lonDim',))
        elif(dim2d==False):
            lat_nc   = outnc.createVariable('lat','f4',('latDim',))
            lon_nc   = outnc.createVariable('lon','f4',('lonDim',))            
        thetao_nc    = outnc.createVariable('thetao','f4',('timeDim','depthDim','latDim','lonDim',))
        tf_nc        = outnc.createVariable('thermalforcing','f4',('timeDim','depthDim','latDim','lonDim',))
            
        time_nc[:]          = timingfull
        depth_nc[:]         = depthout
        if(dim2d):
            lat_nc[:,:]     = latout
            lon_nc[:,:]     = lonout
        elif(dim2d==False):
            lat_nc[:]       = latout
            lon_nc[:]       = lonout
        thetao_nc[:,:,:,:]  = thetafull
        tf_nc[:,:,:,:]      = tffull
        
        depth_nc.units     = 'meter'
        time_nc.units      = 'yr'
        thetao_nc.units    = 'degC'
        tf_nc.units        = 'degC'
        outnc.close()
     
print('End of python job')
os._exit(0) 
    




