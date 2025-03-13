#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 15:21:21 2024

@author: marioadani
"""
import datetime as dt
import sys
sys.path.insert(1, '../UTILS')
from utils import readnc,writenc,qc_mis ,sel_idx,find_nearest, insitu
import numpy as np
import copy

# Edit by the user:
# List of experiments
exp_list=['EXP-S','EXP-A']
# Root input dir model data
rootdir = '../../'
# time period
sdate = dt.datetime(2021,1,1)
edate = dt.datetime(2021,1,31)
sdatestr = sdate.strftime('%Y%m%d')
edatestr = edate.strftime('%Y%m%d')
grd=readnc('../../DATA/grid_medfs831.nc')
max_level,dummy=find_nearest(grd['dep'], 2100)

regs=readnc('../../DATA/STATIC/MedFS_831_regs.nc')
for exp_name in exp_list:
    insitup=insitu()
    print('----------------------'+exp_name+'----------------------------------')
    ndate=sdate
    while ndate<=edate:
        ndatestr = ndate.strftime('%Y%m%d')
        print(ndatestr)
        filename=rootdir+'/'+exp_name+'/POST/INSMIS_'+ndatestr+'.nc'
        ins=readnc(filename)
        # ..TO LOWER CASE
        ins={k.lower(): v for k, v in ins.items()}
        # only the one with valid background
        ins = qc_mis(ins)
        # only the best flag
        idx = np.argwhere(ins['flc']==1)
        ins = sel_idx(ins,idx)
        # READ INCREMENTS
        ins['incr']=np.zeros((len(ins['val'])))
        # for each day and profile avarage all observations within model layers
        profiles=list(set(ins['prof']))
        for pp in profiles:
            idx=pp==ins['prof']
            val=ins['val'][idx]
            bac=ins['bac'][idx]
            tim=ins['tim'][idx]
            lat=ins['lat'][idx]
            lon=ins['lon'][idx]
            dpt=ins['dpt'][idx]
            par=ins['par'][idx]
            prof=ins['prof'][idx]
            x=ins['ib'][0,idx]
            y=ins['jb'][0,idx]
            reg=regs['regs'][y,x]
            incr=ins['incr'][idx]
            # for each parameter
            for param in range(1,3):
                #for each model level
                for lvl in range(max_level):
                    idx=np.where( (dpt>=grd['dep'][lvl]) & (dpt<grd['dep'][lvl+1]) & (param==par) )[0]
                    # if there are values
                    if idx.shape[0]!=0:
                        for ll in insitup.list:
                            eval('insitup.'+ll+'.append(np.mean('+ll+'[idx]))')
        ndate=ndate+dt.timedelta(days=1)
    insitup.to_array()                
    insitup.par=np.int32(insitup.par)                    
    insitup.prof=np.int32(insitup.prof) 
    insitup.reg=np.int32(insitup.reg)
    insitup=insitup.to_dict()
    fileou=rootdir+'/DATA/ins_'+exp_name+'.nc'
    writenc(insitup,fileou) 
