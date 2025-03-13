#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 20 09:52:34 2024

@author: marioadani
"""

import sys
from netCDF4 import Dataset
import numpy as np
import glob
from cftime import  num2date

def read_nc(filename):
    vardict={}
    dimdict={}
    nc=Dataset(filename,mode='r')
    var_list=list(nc.variables.keys())
    dim_list=list(nc.dimensions.keys())
    for var in var_list:
        x=nc.variables[var]
        dim_name=x.dimensions
        val=x[:]
        vardict[var]=[val,dim_name]
    for dim in dim_list:
        x=nc.dimensions[dim].size
        dimdict[dim]=x
    nc.close()
    return vardict,dimdict

def write_nc(dirou,obs_arc,dims_arc,filetype):
    date=num2date(obs_arc['TIM'][0],calendar='standard',units='day since 1950-01-01 00:00')
    sdate=date[0].strftime('%Y%m%d')
    filename=dirou+'/'+filetype+sdate+'.nc'
    nc=Dataset(filename,'w')
    for keys in dims_arc.keys():
        dummy = nc.createDimension(keys, dims_arc[keys])
    for keys in  obs_arc.keys():
        dummy = nc.createVariable(keys, obs_arc[keys][0].dtype, obs_arc[keys][1],zlib=True)
        dummy[:] =  obs_arc[keys][0]
    nc.close()
    return



#---------------------------------------------------------
# Program starts here
dirin = sys.argv[1]
dirou =  sys.argv[2]

#dirin = '/Users/marioadani/SCRATCH/2020/'
#dirou = '/Users/marioadani/SCRATCH/2020'


filetype_list=['INSMIS_','SLAMIS_']

for filetype in filetype_list:
    print('Preparing: '+filetype)
    dirtype=dirin+'/'+filetype
    file_list=sorted(glob.glob(dirtype+"*.NC"))
    first_time=True
    # aggregate all the observation in a domain
    for file_name in file_list:
        print(file_name)
        obs,dims=read_nc(file_name)
        if first_time:
            dims_arc=dims
            key_list=obs.keys()
            obs_arc=obs
            first_time=False
        else:
            for key in key_list:
                if obs[key][1][0]=='OBS':
                    obs_arc[key][0]=np.concatenate((obs_arc[key][0].data, obs[key][0].data),axis=0)
                else:
                    obs_arc[key][0]=np.concatenate((obs_arc[key][0].data, obs[key][0].data),axis=1)
            dims_arc['OBS']=dims_arc['OBS']+dims['OBS']
    #  Initialize obs error array 
    obs_arc['ERR']=[np.zeros(obs_arc['LON'][0].shape),\
                    obs_arc['LON'][1]]
    # Exclude Nan... don't know why are present
    for keys in obs_arc.keys():
        try:
            idx = np.argwhere(np.isnan(obs_arc[keys][0]))
            print('keys : '+keys+' Nan found: '+str(len(idx)))
            obs_arc['FLC'][0][idx] = 0
            idx = np.argwhere(obs_arc[keys][0]>1e20)
            print('keys : '+keys+' Value gt: '+str(len(idx)))
            obs_arc['FLC'][0][idx] = 0
        except TypeError:
            pass
    # finally write file...
    write_nc(dirou,obs_arc,dims_arc,filetype)
