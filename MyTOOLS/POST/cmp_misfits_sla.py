#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 09:08:52 2024

@author: marioadani
"""
import sys
import scipy.interpolate as interp
import datetime as dt
from dateutil.relativedelta import relativedelta
import numpy as np
import copy
import glob
from cftime import num2pydate, date2num
sys.path.insert(1, '../UTILS')
from utils import setup_tides, readnc, find_nearest, \
    unbias, sel_idx, shuffle_obs, qc_mis, writenc

#------------------------------------------------------------------------------
# Edit by the user:
# List of experiments
exp_list=['EXP-A','EXP-S']
# Root input dir model data
rootdir = '../../'
# input dir obs data
dirobs='/Users/marioadani/SCRATCH/MFS831/SLAOBS/'
# time period
sdate = dt.datetime(2021,1,1)
edate = dt.datetime(2021,1,31)
correction = False
increments = False
observations = False
 
# Manage Tides
#filetides='/Users/marioadani/DATA/amppha2D_0_sossheig_20210701_20211231_mod_EAS7.nc'
filetides='../../DATA/TIDE/amppha2D_medfs831.nc'
med_tide=setup_tides(filetides)
minlon=min(med_tide['grid_longitudes'].flatten())
minlat=min(med_tide['grid_latitudes'].flatten())
maxlon=max(med_tide['grid_longitudes'].flatten())
maxlat=max(med_tide['grid_latitudes'].flatten())
x,y=np.meshgrid(med_tide['grid_longitudes'],med_tide['grid_latitudes'])
    
sdatestr = sdate.strftime('%Y%m%d')
edatestr = edate.strftime('%Y%m%d')

# for each experiment
for exp_name in exp_list:
    print('----------------------'+exp_name+'----------------------------------')
        
    ndate=sdate
    while ndate <= edate:
        ndatestr = ndate.strftime('%Y%m%d')

        print(ndatestr)
        filename=rootdir+'/'+exp_name+'/POST/SLAMIS_'+ndatestr+'.nc'
        # READ MISFITS
        sla=readnc(filename)
        # ..TO LOWER CASE
        sla={k.lower(): v for k, v in sla.items()}
        sla['tim']=num2pydate(sla['tim'],units='days since 1950-1-1',calendar='standard')
        sla=shuffle_obs(sla)
        # Put flag to 0 if Nans and values > 1e20 are found    
        sla = qc_mis(sla)
        # Remove value with flag not equal to 1
        idx = np.argwhere(sla['flc']==1)
        sla = sel_idx(sla,idx)
        # READ OBSERVATIONS
        sla['obs_ot']=np.zeros((len(sla['val'])))
        if observations:
            file_list=sorted(glob.glob(dirobs+'/*_europe_*'+ndatestr+'*.nc'))            
            for fileslaobs in file_list:
                obs=readnc(fileslaobs)
                for i in range(len(sla['tim'])):
                    idx=np.where(sla['tim'][i]==obs['time'])[0]
                    if idx.size>0:
                        sla['obs_ot'][i]=obs['ocean_tide'][idx]
        # READ INCREMENTS
        if increments:
            sla['incr']=np.zeros((len(sla['val'])))
            fileincr=rootdir+'/'+exp_name+'/INCREMENTS/'+ndatestr+'_corr_eta.nc'
            incr=readnc(fileincr)
            for i in range(len(sla['tim'])):
                sla['incr'][i]=np.sum(incr['eta'][sla['jb'][:,i],sla['ib'][:,i]]*sla['pq'][:,i])                
        if correction:
            # compute correction 
            cor = sla['res']-sla['val']+sla['bac']
            # sum up the correction
            sla['bac'] = sla['bac']+ cor
        #Compute tide
        sla['mod_ot']=np.zeros(len(sla['bac']))
        for i in range(len(sla['tim'])):
            x,val = find_nearest(med_tide['grid_longitudes'], sla['lon'][i])
            y,val = find_nearest(med_tide['grid_latitudes'], sla['lat'][i])
            if med_tide['msk'][y,x]==1: 
                med_tide['tidecon']=np.array([ [med_tide['M2_Amp'][y,x],0,med_tide['M2_Pha'][y,x],0],\
                                               [med_tide['K1_Amp'][y,x],0,med_tide['K1_Pha'][y,x],0],\
                                               [med_tide['O1_Amp'][y,x],0,med_tide['O1_Pha'][y,x],0],\
                                               [med_tide['S2_Amp'][y,x],0,med_tide['S2_Pha'][y,x],0],\
                                               [med_tide['P1_Amp'][y,x],0,med_tide['P1_Pha'][y,x],0],\
                                               [med_tide['N2_Amp'][y,x],0,med_tide['N2_Pha'][y,x],0],\
                                               [med_tide['Q1_Amp'][y,x],0,med_tide['Q1_Pha'][y,x],0],\
                                               [med_tide['K2_Amp'][y,x],0,med_tide['K2_Pha'][y,x],0] ],\
                                           dtype='>f8')
                med_tide['lat']=med_tide['grid_latitudes'][y]
                sla['mod_ot'][i]=med_tide(np.array([sla['tim'][i]])) 
                
        # remove tide
        sla['bac']=sla['bac']-sla['mod_ot']
        sla['res'] = sla['val']- sla['bac']
        # unbias residuals
        sla = unbias(sla,40)
        # Remove value with flag not equal to 1 that comes from unbias
        idx = np.argwhere(sla['flc']==1)
        sla = sel_idx(sla,idx)
        if ndate == sdate:
            arc_sla=copy.deepcopy(sla)
        else:
            for key in sla.keys():
                if key=='ib' or key=='jb' or key=='pq' or key=='tb' or key=='sb' or key=='bcp':
                    arc_sla[key]=np.append(arc_sla[key],sla[key],axis=1)
                else:
                    arc_sla[key]=np.append(arc_sla[key],sla[key],axis=0)
                                
        ndate=ndate+dt.timedelta(days=1)
    fileou=rootdir+'/DATA/sla_'+exp_name+'.nc'
    arc_sla['tim']=date2num(arc_sla['tim'],units='days since 1950-1-1',calendar='standard')
    writenc(arc_sla,fileou)
