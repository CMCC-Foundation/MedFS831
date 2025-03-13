#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 08:20:33 2024

@author: marioadani
"""
import sys
sys.path.insert(1, '../UTILS')
from   utils import  readnc
import numpy as np
import matplotlib.pyplot as plt
from   netCDF4 import Dataset

#####################   
# Input file
#####################   
inputfile  = '../../medsea-nemo42/tools/DOMAINcfg/mesh_mask.nc'
outputfile = '../../medsea-nemo42/tools/DOMAINcfg/coordinates.bdy.nc'

#####################   
# read file
#####################   
src = readnc(inputfile)

#####################   
#create map of index 
#####################   
src['nbit'] = np.zeros( src['e1t'].shape ,dtype=np.int64)
src['nbjt'] = np.zeros( src['e1t'].shape ,dtype=np.int64)
src['nbiu'] = np.zeros( src['e1u'].shape ,dtype=np.int64)
src['nbju'] = np.zeros( src['e1u'].shape ,dtype=np.int64)
src['nbiv'] = np.zeros( src['e1v'].shape ,dtype=np.int64)
src['nbjv'] = np.zeros( src['e1v'].shape ,dtype=np.int64)

for i in range(src['nbit'].shape[2]):
    src['nbit'][:,:,i] = i+1
    src['nbiu'][:,:,i] = i+1
    src['nbiv'][:,:,i] = i+1
for j in range(src['nbit'].shape[1]):
    src['nbjt'][:,j,:] = j+1
    src['nbju'][:,j,:] = j+1
    src['nbjv'][:,j,:] = j+1
plt.pcolor( src['nbit'].squeeze());plt.colorbar()

#####################   
# List of variables to be processed for T,U,V grid respectively
#####################   
vlistt = ['nbit','nbjt','e1t','e2t','glamt','gphit']
vlistu = ['nbiu','nbju','e1u','e2u','glamu','gphiu']
vlistv = ['nbiv','nbjv','e1v','e2v','glamv','gphiv']

#####################   
# Extract boundaries
#####################   
out={}
for v in vlistt:
    #East
    out[v] = src[v][:,:,1].squeeze()                                                    # 1 = 2i
    out['tmask'] = src['tmask'][:,0,:,1].squeeze()
    #North
    out[v] = np.concatenate( (out[v],src[v][:,-2,:].squeeze()) )                        # -2 = jmt-1
    out['tmask'] = np.concatenate( (out['tmask'],src['tmask'][:,0,-2,:].squeeze()) )
    #West
    out[v] = np.concatenate( (out[v],src[v][:,:,-2].squeeze()) )                        # -2 = imt-1
    out['tmask'] = np.concatenate( (out['tmask'],src['tmask'][:,0,:,-2].squeeze()) )
    #South
    out[v] = np.concatenate( (out[v],src[v][:,1,:].squeeze()) )                         # 1 = 2 j  
    out['tmask'] = np.concatenate( (out['tmask'],src['tmask'][:,0,1,:].squeeze()) )
# U grid
for v in vlistu:
    #East
    out[v] = src[v][:,:,1].squeeze()                                                    # 1 = 2
    out['umask'] = src['umask'][:,0,:,1].squeeze()
    #North
    out[v] = np.concatenate( (out[v],src[v][:,-2,:].squeeze()) )                        # -2 =jmt-1
    out['umask'] = np.concatenate( (out['umask'],src['umask'][:,0,-2,:].squeeze()) )
    #West
    out[v] = np.concatenate( (out[v],src[v][:,:,-3].squeeze()) )                        # -2 =i,imt-1
    out['umask'] = np.concatenate( (out['umask'],src['umask'][:,0,:,-3].squeeze()) )
    #South
    out[v] = np.concatenate( (out[v],src[v][:,1,:-1].squeeze()) )                       # 1 = 2j
    out['umask'] = np.concatenate( (out['umask'],src['umask'][:,0,1,:].squeeze()) )
# V grid
for v in vlistv:
    #East
    out[v] = src[v][:,:,1].squeeze()                                                    # 1 = 2i
    out['vmask'] = src['vmask'][:,0,:,1].squeeze()
    #North
    out[v] = np.concatenate( (out[v],src[v][:,-3,:].squeeze()) )                        # -3 = jmt-2
    out['vmask'] = np.concatenate( (out['vmask'],src['vmask'][:,0,-3,:].squeeze()) )
    #West
    out[v] = np.concatenate( (out[v],src[v][:,:,-2].squeeze()) )                        # -2 = imt-1
    out['vmask'] = np.concatenate( (out['vmask'],src['vmask'][:,0,:,-2].squeeze()) )
    #South
    out[v] = np.concatenate( (out[v],src[v][:,1,:].squeeze()) )                         # 1 = 2j
    out['vmask'] = np.concatenate( (out['vmask'],src['vmask'][:,0,1,:].squeeze()) )

#####################   
# Remove land points    
#####################   
# T grid
seapoint=np.where(out['tmask']==1)
for v in vlistt:
    out[v] = out[v][seapoint]
# U grid
seapoint=np.where(out['umask']==1)
for v in vlistu:
    out[v] = out[v][seapoint]
# V grid
seapoint=np.where(out['vmask']==1)
for v in vlistv:
    out[v] = out[v][seapoint]


#####################   
# Remove duplicate  
#####################    
# T grid
coord = np.array((out['nbit'],out['nbjt'])).T
_, idx = np.unique(coord, axis=0, return_index=True)
for v in vlistt:
    out[v] = out[v][idx]
# U grid
coord = np.array((out['nbiu'],out['nbju'])).T
_, idx = np.unique(coord, axis=0, return_index=True)
for v in vlistu:
    out[v] = out[v][idx]
# V grid
coord = np.array((out['nbiv'],out['nbjv'])).T
_, idx = np.unique(coord, axis=0, return_index=True)
for v in vlistv:
    out[v] = out[v][idx]

#####################   
# Write coordinates
##################### 
yb = 1
xbT = out['nbjt'].shape[0]
xbU = out['nbju'].shape[0]
xbV = out['nbjv'].shape[0]

# Open file
nc=Dataset(outputfile,mode='w',format='NETCDF4')

# Define Dimensions
ybt_dim = nc.createDimension('yb', yb)
xbT_dim = nc.createDimension('xbT', xbT)
xbU_dim = nc.createDimension('xbU', xbU)
xbV_dim = nc.createDimension('xbV', xbV)

# Title    
nc.title='MED TEST CASE  Coordinates for BDY'
nc.author='M. Adani'

# Create Variables & fill in variables
for v in vlistt[0:2]:
    dummy = nc.createVariable(v, np.int32, ('yb','xbT'),zlib=True)
    dummy[:,:] = out[v]
dummy= nc.createVariable('nbrt', np.int32, ('yb','xbT'),zlib=True)
dummy[:,:] = 1
for v in vlistu[0:2]:
    dummy = nc.createVariable(v, np.int32, ('yb','xbU'),zlib=True)
    dummy[:,:] = out[v]
dummy = nc.createVariable('nbru', np.int32, ('yb','xbU'),zlib=True)
dummy[:,:] = 1
for v in vlistv[0:2]:
    dummy = nc.createVariable(v, np.int32, ('yb','xbV'),zlib=True)
    dummy[:,:] = out[v]
dummy = nc.createVariable('nbrv', np.int32, ('yb','xbV'),zlib=True)
dummy[:,:] = 1

for v in vlistt[3:]:
    dummy = nc.createVariable(v, np.float64, ('yb','xbT'),zlib=True)
    dummy[:,:] = out[v]
for v in vlistu[3:]:
    dummy = nc.createVariable(v, np.float64, ('yb','xbU'),zlib=True)
    dummy[:,:] = out[v]
for v in vlistv[3:]:
    dummy = nc.createVariable(v, np.float64, ('yb','xbV'),zlib=True)
    dummy[:,:] = out[v]
nc.close()

