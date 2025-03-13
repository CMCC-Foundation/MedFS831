#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 16:34:17 2025

@author: marioadani
"""

import sys
sys.path.insert(1, '../UTILS')
from utils import readnc, interp2d,SeaOverLand
from netCDF4 import Dataset
import numpy as np

filegrd   = '../../medsea-nemo42/tools/DOMAINcfg/domain_cfg.nc'
filemesh  = '../../medsea-nemo42/tools/DOMAINcfg/mesh_mask.nc'
filecoat  = '../../DATA/coastal_distance.nc'
fileout   = '../../DATA/grid_medfs831.nc'

grd = readnc(filegrd)
msh = readnc(filemesh)
cst = readnc(filecoat)


out={} #lon,lat,dep,tmsk,topo
out['lon']  = grd['nav_lon'].squeeze() 
out['lat']  = grd['nav_lat'].squeeze() 
out['dep']  = grd['nav_lev'].squeeze() 
out['dx']   = grd['e1t'].squeeze() 
out['dy']   = grd['e2t'].squeeze() 
out['dz']   = grd['e3t_1d'].squeeze()   
out['tmsk'] = msh['tmask'].squeeze() 
out['topo'] = grd['bathy_metry'].squeeze() 
out['dcoast'] = cst['coast'].squeeze() 
out['dcoast2d'] = cst['coast2d'].squeeze() 

#write file
imt = out['lon'].shape[1]
jmt = out['lat'].shape[0]
kmt = out['dep'].shape[0]
nc=Dataset(fileout,mode='w',format='NETCDF4')
# Define Dimensions
jmt_dim = nc.createDimension('jm', jmt)
imt_dim = nc.createDimension('im', imt)
kmt_dim = nc.createDimension('km', kmt)

# Title    
nc.title='MedFS831 test case'
nc.author='M. Adani'
# Create Variables
lonvar = nc.createVariable('lon', np.single, ('jm','im'),zlib=True)
lonvar.standard_name='longitude'
lonvar.long_name='Longitude'
lonvar.units ='degrees_east'

latvar = nc.createVariable('lat', np.single, ('jm','im'),zlib=True)
latvar.standard_name='latitude'
latvar.long_name='Latitude'
latvar.units='degrees_north'

depvar = nc.createVariable('dep', np.single, ('km'),zlib=True)
depvar.standard_name='depth'
depvar.long_name='depth'
depvar.units='meters'

dxvar = nc.createVariable('dx', np.single, ('jm','im'),zlib=True)
dxvar.coordinates='lat lon'
dxvar.standard_name='grid distance in longitude direction'
dxvar.long_name='grid distance in longitude direction'
dxvar.units='meters'

dyvar = nc.createVariable('dy', np.single, ('jm','im'),zlib=True)
dyvar.coordinates='lat lon'
dyvar.standard_name='grid distance in latitude direction'
dyvar.long_name='grid distance in latitude direction'
dyvar.units='meters'

dzvar = nc.createVariable('dz', np.single, ('km'),zlib=True)
dzvar.coordinates='depth'
dzvar.standard_name='grid distance in vertical direction'
dzvar.long_name='grid distance in vertical direction'
dzvar.units='meters'

topovar = nc.createVariable('topo', np.single, ('jm','im'),zlib=True)
topovar.coordinates='lat lon'
topovar.standard_name='model bathymetry'
topovar.long_name='model bathymetry'
topovar.units='meters'

cst2dvar = nc.createVariable('dcoast2d', np.single, ('jm','im'),zlib=True)
cst2dvar.coordinates='lat lon'
cst2dvar.standard_name='distance from coast'
cst2dvar.long_name='distance from coast'
cst2dvar.units='meters'

cst3dvar = nc.createVariable('dcoast', np.single, ('km','jm','im'),zlib=True)
cst3dvar.coordinates='depth lat lon'
cst3dvar.standard_name='distance from coast'
cst3dvar.long_name='distance from coast'
cst3dvar.units='meters'


mskvar = nc.createVariable('tmsk', np.single, ('km','jm','im'),zlib=True)
mskvar.coordinates='depth lat lon'
mskvar.standard_name='land mask on T point'
mskvar.long_name='land mask on T point'
mskvar.units='#'

# Fill in variables
lonvar[:,:] = out['lon'] 
latvar[:,:] = out['lat'] 
depvar[:]   = out['dep'] 
dxvar[:,:]  = out['dx'] 
dyvar[:,:]  = out['dy'] 
dzvar[:]    = out['dz'] 
topovar[:,:]= out['topo'] 
mskvar[:,:,:]= out['tmsk'] 
cst3dvar[:,:,:]= out['dcoast'] 
cst2dvar[:,:]= out['dcoast2d'] 

nc.close()

