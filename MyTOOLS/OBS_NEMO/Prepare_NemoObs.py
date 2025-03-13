#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 14:11:47 2024

@author: marioadani
"""
import os
import sys
import glob
import numpy as np
import datetime as dt
import sys
sys.path.insert(1, '../UTILS')
from prep_utils import readnc, insitu, get_familycode, get_rank,    \
    get_cell, read_layout, defdom, get_weight_2d, writenc, get_sat, \
    get_kty, dep_to_p, potemp, SeaLevelAnomaly, sel_idx, sorting
from cftime import date2num, num2date    

def prep_insitu(dirin):
    #
    # DEFINE GLOBAL LON/LAT
    #
    longrd_glob = grd['lon']['val'].data
    latgrd_glob = grd['lat']['val'].data
    depgrd      = grd['dep']['val'].data
    #
    # LIST OF INPUT FILE
    #
    list_file = glob.glob(dirin+'/*.nc');list_file.sort()
    #
    # DEFINE STRUCTURE OF INSITU OBS
    #
    ins = insitu()
    #
    # INITIALIZE NUMBER OF PROFILE
    #
    kprof     = 0           # initialize number of profile
    start_obs = 1           # start number of obs
    kmask     = 12          # Mask problems
    for file in list_file:
        print(file)
        kprof = kprof + 1 # assume one profile per file
        #
        # READ OBS FILE
        #
        varin = readnc(file)
        #
        # DEFINE SOM OBS VARIABLE
        #
        lonobs = varin['LONGITUDE']['val'].data
        latobs = varin['LATITUDE'] ['val'].data
        dptobs = varin['DEPH']     ['val'].data
        an_jultime = date2num(an_time,varin['TIME']['units'])
        tdist = varin['TIME']['val'] - an_jultime
        #
        # FIND PROCESSOR (RANK,IPROC,JPROC)
        #
        rank = get_rank(layout,dom,longrd_glob,latgrd_glob,lonobs,latobs)
        #
        # OBS NOT INSIDE OF THE MODEL DOMAIN
        #
        if rank == -999:
            continue
        area  = rank + 1
        iproc = layout['ii'][rank]
        jproc = layout['ij'][rank]
        #
        #DEFINE LOCAL LON/LAT
        #
        # -1 to take into account python start from 0
        # xend , yend do not have -1 because python does not consider the last point
        longrd_loc = longrd_glob[dom['ystart'][rank]-1:dom['yend'][rank],dom['xstart'][rank]-1:dom['xend'][rank]]
        latgrd_loc = latgrd_glob[dom['ystart'][rank]-1:dom['yend'][rank],dom['xstart'][rank]-1:dom['xend'][rank]]
        #
        # DEFINE LOCAL MASK
        #
        maskgrd_loc = grd['tmsk']['val'][: , dom['ystart'][rank]-1:dom['yend'][rank],dom['xstart'][rank]-1:dom['xend'][rank]]
        #
        # FIND CELL WITHIN PROCESSOR
        #
        igrd_loc_py, jgrd_loc_py, lon, lat = get_cell(longrd_loc,latgrd_loc,lonobs, latobs)
        #
        # FIND LOCAL INDEX i,j OF THE GRID 
        #
        # +1 to take into account fortran start from 1
        igrd_loc_fortran = igrd_loc_py+1
        jgrd_loc_fortran = jgrd_loc_py+1
        #
        # FIND GLOBAL INDEXES
        #
        igrd_glob_py = igrd_loc_py + dom['xstart'][rank]-1
        jgrd_glob_py = jgrd_loc_py + dom['ystart'][rank]-1
        igrd_glob_fortran = igrd_glob_py+1
        jgrd_glob_fortran = jgrd_glob_py+1
        #
        # FIND INTERPOLATION 2D WEIGHTS 
        #
        weighth = get_weight_2d(lon,lat,lonobs,latobs)
        #
        # DEFINE NUMBER OF OBS IN PROFILE
        #
        nobs   = len(varin['TEMP']['val'])
        #
        # FIND K LEVEL
        #
        kgrd_fortran = np.zeros((nobs),dtype=np.int8)
        weightv      = np.zeros((nobs))
        pq           = np.zeros((nobs,8))
        flc          = np.ones ((nobs))
        eve          = np.zeros((nobs))
        for iobs in range(nobs):
            if ( dptobs[iobs] < depgrd[0] ):
                k_py = 0
            else:
                k_py = np.where(depgrd<=dptobs[iobs])[0][-1]
            weightv[iobs] = max(1e-12,abs(dptobs[iobs]-depgrd[k_py]))/ \
                            max(1e-12,abs(depgrd[k_py+1]-depgrd[k_py]))
            # +1 to take into account fortran start from 1
            kgrd_fortran[iobs] = int(k_py + 1)
            kp1_py             = int(kgrd_fortran[iobs])
            r1 = weightv[iobs]
            r2 = 1 - weightv[iobs]
            summask = np.sum(maskgrd_loc[k_py:kp1_py+1,jgrd_loc_py:jgrd_loc_py+2,igrd_loc_py:igrd_loc_py+2].squeeze())
            if  ( summask < 1 ): #at least 1 points
                flc[iobs] = 0
                eve[iobs] = kmask
            else :
                pres = dep_to_p(dptobs[iobs],latobs)
                varin['TEMP']['val'][iobs] = potemp(varin['PSAL']['val'][iobs], varin['TEMP']['val'][iobs], pres, 0)
                #
                # FIND INTERPOLATION 3D WEIGHTS 
                #
                pq[iobs,0] = r2 * weighth[0] * maskgrd_loc[k_py  ,jgrd_loc_py  ,igrd_loc_py  ]
                pq[iobs,1] = r2 * weighth[1] * maskgrd_loc[k_py  ,jgrd_loc_py  ,igrd_loc_py+1]
                pq[iobs,2] = r2 * weighth[2] * maskgrd_loc[k_py  ,jgrd_loc_py+1,igrd_loc_py+1]
                pq[iobs,3] = r2 * weighth[3] * maskgrd_loc[k_py  ,jgrd_loc_py+1,igrd_loc_py  ]
                pq[iobs,4] = r1 * weighth[0] * maskgrd_loc[kp1_py,jgrd_loc_py  ,igrd_loc_py  ]
                pq[iobs,5] = r1 * weighth[1] * maskgrd_loc[kp1_py,jgrd_loc_py  ,igrd_loc_py+1]
                pq[iobs,6] = r1 * weighth[2] * maskgrd_loc[kp1_py,jgrd_loc_py+1,igrd_loc_py+1]
                pq[iobs,7] = r1 * weighth[3] * maskgrd_loc[kp1_py,jgrd_loc_py+1,igrd_loc_py  ]  
                pq[iobs,:] =  pq[iobs,:] / np.sum(pq[iobs,:])
        #
        # ARCHIVING
        #
        fmcode = get_familycode(varin['global']['family_code'])
        kty    = get_kty(fmcode)
        for i in range (1,3):
            if i == 2: # TEMPERATURE
                ins.val.extend(varin['TEMP']['val'].data.tolist()                )
            if i == 1: # SALINITY
                ins.val.extend(varin['PSAL']['val'].data.tolist()                )
            ins.par.extend  ( [ i                                       ] * nobs )
            ins.tim.extend  ( [ varin['TIME']['val'].data.tolist()      ] * nobs )
            ins.lat.extend  ( [ varin['LATITUDE']['val'].data.tolist()  ] * nobs )
            ins.lon.extend  ( [ varin['LONGITUDE']['val'].data.tolist() ] * nobs ) 
            ins.dpt.extend  (   varin['DEPH']['val'].data.tolist()               ) 
            ins.otype.extend( [ fmcode                                  ] * nobs )
            ins.plno.extend ( [ varin['global']['wmo_platform_code']    ] * nobs )
            ins.prof.extend ( [ kprof                                   ] * nobs )
            ins.flc.extend  (   flc.tolist()                                     )
            ins.tdist.extend( [ tdist                                   ] * nobs )
            ins.nind.extend ( [*range(start_obs,nobs+start_obs)         ]        )
            ins.ino.extend  ( [ 0                                       ] * nobs )
            ins.eve.extend  (   eve.tolist()                                     )
            ins.inst.extend ( [ 0                                       ] * nobs )
            ins.bia.extend  ( [ 0.0                                     ] * nobs )
            ins.kb.extend   (   kgrd_fortran.tolist()                            )
            ins.sd1.extend  ( [ area                                    ] * nobs )
            ins.sd2.extend  ( [ area                                    ] * nobs )
            ins.pq.extend   ( pq.tolist()                                        )
            ins.rb.extend   ( weightv.tolist()                                   )
            ins.kty.extend  ( [ kty                                     ] * nobs )
            ins.prind.extend( [ 0                                       ] * nobs )
            ins.ib.extend   ( [[ igrd_glob_fortran, igrd_glob_fortran+1, igrd_glob_fortran+1, igrd_glob_fortran   ]] * nobs )
            ins.jb.extend   ( [[ jgrd_glob_fortran, jgrd_glob_fortran  , jgrd_glob_fortran+1, jgrd_glob_fortran+1 ]] * nobs )
            ins.moi.extend  ( [[ 0                , 0                  , 0                  , 0                   ]] * nobs )
            ins.moj.extend  ( [[ 0                , 0                  , 0                  , 0                   ]] * nobs )

        start_obs=start_obs + nobs 
        
        #
        # DELETE SOME IMPORTANT VARIABLES FOR NEXT COMPUTATION
        #
        del rank,   igrd_glob_py,  jgrd_glob_py, igrd_loc_py, jgrd_loc_py,  k_py, kp1_py
    # 
    # DEFINE THE REMAINING DIMENSIONS
    #    
    ins.obs = len(ins.par)
    ins.levs= depgrd.shape[0]
    #
    # CONVERT LIST TO ARRAY
    #
    ins.list2array()
    #
    # SORTING THE ARRAY IN TEMPORAL ORDER
    #
    idx = np.argsort(ins.tim)
    ins = sorting(ins,idx,'ins')

    rank_list = list(set(ins.sd1))
    obj_list  = list(ins.__dict__.keys())
    #
    # COPY local array
    #
    ins_local = insitu()
    ins_local.list2array()
    #
    #WRITE NETCDF PER RANK
    #
    for area in rank_list:
        idx = np.where(ins.sd1==area)[0]
        for var in obj_list:
            varu = var.upper()
            if isinstance(eval('ins.'+var),int):
                if varu!='OBS':
                    exec('ins_local.'+var+' = ins.'+var )
                else:
                    ins_local.obs = idx.shape[0]
            else:
                if varu=='PLNO':
                    exec('ins_local.'+var+'=np.array(ins.'+var+'[idx,:])')
                elif varu =='PQ' or varu =='MOI'  or varu =='MOJ' or varu =='IB' or varu =='JB' :
                    exec('ins_local.'+var+'=np.array(ins.'+var+'[:,idx])')
                else:    
                    exec('ins_local.'+var+'=np.array(ins.'+var+'[idx])')
        fileou=dirou+'/INSOBS_'+str(area-1).zfill(4)+'.NC'
        writenc(fileou,ins_local)
    return ins

def prep_sla(dirin):
    #
    # DEFINE GLOBAL LON/LAT
    #
    longrd_glob = grd['lon']['val'].data
    latgrd_glob = grd['lat']['val'].data
    depgrd      = grd['dep']['val'].data
    kmask       = 12          # Mask problems
    #
    # LIST OF INPUT FILE
    #
    list_file = glob.glob(dirin+'/*');list_file.sort()
    #
    # DEFINE STRUCTURE OF INSITU OBS
    #
    sla = SeaLevelAnomaly()
    kobs   = 1
    for file in list_file:
        print(file)
        #
        # READ OBS FILE
        #
        varin = readnc(file)
        #
        # REMOVE OBS OUTSIDE TIME WINDOS 
        #
        timeobs = num2date(varin['time']['val'],varin['time']['units'])
        idx = np.where((timeobs >= st_time) & (timeobs <= an_time) )
        varin = sel_idx(varin,idx)
        #
        # DEFINE SOME OBS VARIABLE
        #
        lonobs = varin['longitude']['val'].data;lonobs[lonobs>180]=lonobs[lonobs>180]-360
        latobs = varin['latitude'] ['val'].data
        an_jultime = date2num(an_time,varin['time']['units'])
        tdist = varin['time']['val'] - an_jultime
        nobs   = len(varin['longitude']['val'])
        for iobs in range(nobs):
            pq = np.zeros((4))
            #
            # FIND PROCESSOR (RANK,IPROC,JPROC)
            #
            rank = get_rank(layout,dom,longrd_glob,latgrd_glob,lonobs[iobs],latobs[iobs])
            #
            # OBS NOT INSIDE OF THE MODEL DOMAIN
            #
            if rank == -999:
                continue
            area  = rank + 1
            iproc = layout['ii'][rank]
            jproc = layout['ij'][rank]
            #
            #DEFINE LOCAL LON/LAT
            #
            # -1 to take into account python start from 0
            # xend , yend do not have -1 because python does not consider the last point
            longrd_loc = longrd_glob[dom['ystart'][rank]-1:dom['yend'][rank],dom['xstart'][rank]-1:dom['xend'][rank]]
            latgrd_loc = latgrd_glob[dom['ystart'][rank]-1:dom['yend'][rank],dom['xstart'][rank]-1:dom['xend'][rank]]
            #
            # DEFINE LOCAL MASK
            #
            maskgrd_loc = grd['tmsk']['val'][0 , dom['ystart'][rank]-1:dom['yend'][rank],dom['xstart'][rank]-1:dom['xend'][rank]].squeeze()
            #
            # FIND CELL WITHIN PROCESSOR
            #
            igrd_loc_py, jgrd_loc_py, lon, lat = get_cell(longrd_loc,latgrd_loc,lonobs[iobs], latobs[iobs])
            #
            # FIND LOCAL INDEX i,j OF THE GRID 
            #
            # +1 to take into account fortran start from 1
            igrd_loc_fortran = igrd_loc_py+1
            jgrd_loc_fortran = jgrd_loc_py+1
            #
            # FIND GLOBAL INDEXES
            #
            igrd_glob_py = igrd_loc_py + dom['xstart'][rank]-1
            jgrd_glob_py = jgrd_loc_py + dom['ystart'][rank]-1
            igrd_glob_fortran = igrd_glob_py+1
            jgrd_glob_fortran = jgrd_glob_py+1
            flc = 1
            eve = 0
            #
            # EXCLUDE DATA ON THE ATLANTIC BOX
            #
            atlbox = regs['regs']['val'].mask[jgrd_glob_py,igrd_glob_py]==True
            if (atlbox):
                continue
            #
            # EXCLUDE DATA ON LAND
            #
            summask = np.sum(maskgrd_loc[jgrd_loc_py:jgrd_loc_py+2,igrd_loc_py:igrd_loc_py+2].squeeze())
            if  ( summask < 3 ): #at least 3 points
                continue
            #
            # FIND INTERPOLATION 2D WEIGHTS 
            #
            weighth = get_weight_2d(lon,lat,lonobs[iobs],latobs[iobs])
            pq[0] = weighth[0] * maskgrd_loc[jgrd_loc_py  ,igrd_loc_py  ]
            pq[1] = weighth[1] * maskgrd_loc[jgrd_loc_py  ,igrd_loc_py+1]
            pq[2] = weighth[2] * maskgrd_loc[jgrd_loc_py+1,igrd_loc_py+1]
            pq[3] = weighth[3] * maskgrd_loc[jgrd_loc_py+1,igrd_loc_py  ]
            pq[:] = pq[:] / np.sum(pq)
            dpt  = np.min(grd['topo']['val'][jgrd_glob_py:jgrd_glob_py+2,igrd_glob_py:igrd_glob_py+2])
            #
            # EXCLUDE DATA TOO SHALLOW
            #
            if ( dpt < min_sla_dpt ):
                continue
            platform = os.path.basename(file);platform = platform[:-2]
            ksat = get_sat(platform)
            #
            # EXCLUDE NOT REALISTIC DATA (+-10m)
            #
            if ( abs(varin['SLApDACpOT']['val'].data[iobs]) > 10 ):
                continue
            #
            # ARCHIVING
            #
            sla.tim.extend  ( [ varin['time']['val'].data[iobs]      ])
            sla.lon.extend  ( [ varin['longitude']['val'].data[iobs] ])
            sla.lat.extend  ( [ varin['latitude']['val'].data[iobs]  ])
            sla.track.extend( [ varin['track']['val'].data[iobs]     ])
            sla.tdist.extend( [ tdist[iobs]                          ])
            sla.nind.extend ( [ kobs                                 ])
            sla.pq.append   (  pq.tolist()                            )
            sla.ino.extend  ( [ 0                                    ])
            sla.dpt.extend  ( [ dpt                                  ])
            sla.sd1.extend  ( [ area                                 ])
            sla.sd2.extend  ( [ area                                 ])
            sla.eve.extend  ( [ eve                                  ])
            sla.flc.extend  ( [ flc                                  ])
            sla.ksat.extend ( [ ksat                                 ])
            sla.val.extend  ( [ varin['SLApDACpOT']['val'].data[iobs]   ])
            sla.bia.extend  ( [ 0.0                                  ])
            sla.bot.extend  ( [ -1                                   ])
            sla.ib.append   ( [ igrd_glob_fortran, igrd_glob_fortran+1, igrd_glob_fortran+1, igrd_glob_fortran   ])
            sla.jb.append   ( [ jgrd_glob_fortran, jgrd_glob_fortran  , jgrd_glob_fortran+1, jgrd_glob_fortran+1 ])
            sla.moi.append  ( [ 0                , 0                  , 0                  , 0                   ])
            sla.moj.append  ( [ 0                , 0                  , 0                  , 0                   ])
            
            kobs += 1
    # 
    # DEFINE THE REMAINING DIMENSIONS
    #    
    sla.obs = len(sla.val)        
    sla.levs= depgrd.shape[0]
    #
    # CONVERT LIST TO ARRAY
    #    
    sla.list2array()
    #
    # SORTING THE ARRAY IN TEMPORAL ORDER
    #
    idx = np.argsort(sla.tim)
    sla = sorting(sla,idx,'sla')
    #
    # COPY local array
    #
    sla_local = SeaLevelAnomaly()
    sla_local.list2array()
    rank_list = list(set(sla.sd1))
    obj_list  = list(sla.__dict__.keys())
    #
    #WRITE NETCDF PER RANK
    #
    for area in rank_list:
        idx = np.where(sla.sd1==area)[0]
        for var in obj_list:
            varu = var.upper()
            if isinstance(eval('sla.'+var),int):
                if varu!='OBS':
                    exec('sla_local.'+var+' = sla.'+var )
                else:
                    sla_local.obs = idx.shape[0]
            else:
                if varu =='PQ' or varu =='MOI'  or varu =='MOJ' or varu =='IB' or varu =='JB' :
                    exec('sla_local.'+var+'=np.array(sla.'+var+'[:,idx])')
                else:    
                    exec('sla_local.'+var+'=np.array(sla.'+var+'[idx])')
        fileou=dirou+'/SLAOBS_'+str(area-1).zfill(4)+'.NC'
        writenc(fileou,sla_local)
    return sla

#--------------------------------------------------------------------------------------
#--------------- MAIN STARTS FROM HERE ------------------------------------------------
#--------------------------------------------------------------------------------------

#
# INPUT DIRECTORY
#
dirin= sys.argv[1] 
dirin_ins=sys.argv[2]
dirin_sla=sys.argv[3]
dirou=sys.argv[4]
sdate=sys.argv[5]
#
# FILEOUTPUT
#
fileout_ins = dirou+'/obs_ins.nc'
fileout_sla = dirou+'/obs_sla.nc'
#
# MIN SLA DPT
#
min_sla_dpt = 0.
#
# TIME WINDOW IN HOURS
#
time_window = 24
#
# COMPUTE START 
#
year  = int(sdate[0:4])
month = int(sdate[4:6])
day   = int(sdate[6:8])
st_time = dt.datetime(year,month,day,0,0)
#
# ANALYSIS TIME 
#
an_time = st_time + dt.timedelta(hours=time_window)
#
# GRID FILE
#
#MFS24 filegrd=dirin+'/GRID_BS.nc'
filegrd=dirin+'/grid_medfs831.nc'
#
# GRID FILE
#
#MFS24 filereg=dirin+'/MFS_24_mdt_mask_regs.nc'
filereg=dirin+'/MedFS_831_regs.nc'
#
# LAYOUT FILE
#
#MFS24 layoutfile=dirin+'/layout.dat'
layoutfile=dirin+'/layout.dat'
#
# NUMBER OF HALO POINTS
#
nn_hls = 1 # number of halo point define in NEMO namelist
#
# READ GRID
#
grd = readnc(filegrd)
#
# READ MASK OF THE REGIONS 
#
regs = readnc(filereg); 
regs['regs']['val'].mask[regs['regs']['val'].data==17]=True
#
# READ LAYOUT
#
layout = read_layout(layoutfile)
#
# DEFINE DOMAIN
#
dom = defdom(layout,nn_hls)
#
# PREP INSITU

ins = prep_insitu(dirin_ins)

# PREP SLA
#
sla = prep_sla(dirin_sla)
#
# OPEN ASCII FILE
#
filetab = dirou+'/OBS_TAB.DAT'
f = open(filetab, "w")
#
# FILL IN ASCII FILE
#   
# area, INS , SLA, SST, SSS
for proc in layout['rank']:
    area = proc + 1
    idx_ins = np.where(ins.sd1==area)[0]
    idx_sla = np.where(sla.sd1==area)[0]
    f.writelines(['{0:12d}'.format(area)            ,\
                  '{0:12d}'.format(idx_ins.shape[0]),\
                  '{0:12d}'.format(idx_sla.shape[0]),\
                  '{0:12d}'.format(0)               ,\
                  '{0:12d}'.format(0)+'\n'])
f.close()
#
# WRITE COMPLETE NC FILE
#    
writenc(fileout_ins,ins)
writenc(fileout_sla,sla)

