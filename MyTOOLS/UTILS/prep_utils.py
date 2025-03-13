#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 14:16:09 2024

@author: marioadani
"""
from netCDF4 import Dataset, chartostring
import numpy as np
import re


def sorting(invar,idx,typ):
    if typ == 'sla':
        outvar = SeaLevelAnomaly()
    if typ == 'ins':
        outvar = insitu()
    outvar.list2array()
    obj_list  = list(invar.__dict__.keys())
    for var in obj_list:
        if isinstance(eval('invar.'+var),int):
            exec('outvar.'+var+'=invar.'+var)
        else:
            if eval('invar.'+var+'.ndim')==1:
                exec('outvar.'+var+'=invar.'+var+'[idx]')
            elif eval('invar.'+var+'.ndim')==2:
                if var=='plno':
                   exec('outvar.'+var+'=invar.'+var+'[idx,:]')
                else:
                    exec('outvar.'+var+'=invar.'+var+'[:,idx]')
    return outvar

def get_weight_2d(lon,lat,lonobs,latobs):
    weighth = np.zeros((4))
    n = 0
    for lo,la in zip(lon,lat):
            dist = haversine(lonobs,latobs,lo,la)
            # weighth[n] = dist
            weighth[n] = 1.0/max(1e-12,dist)
            n += 1
    weighth = weighth/sum(weighth)
    return weighth

def get_rank(layout,dom,longrd_glob,latgrd_glob,lonobs,latobs):
    rankout = -999
    # domain of the tile 
    #   P4(lon[3],lat[3])   P3(lon[2],lat[2])
    #
    #   P1(lon[0],lat[0])   P2(lon[1],lat[1])
    lon = np.zeros((4))
    lat = np.zeros((4))
    for rank in layout['rank']:
        # -1 to take into account python start from 0
        lon[0] = longrd_glob[dom['ystart'][rank]-1, dom['xstart'][rank]-1]
        lon[1] = longrd_glob[dom['ystart'][rank]-1, dom['xend']  [rank]-1]
        lon[2] = longrd_glob[dom['yend']  [rank]-1, dom['xend']  [rank]-1]
        lon[3] = longrd_glob[dom['yend']  [rank]-1, dom['xstart'][rank]-1]
        lat[0] = latgrd_glob[dom['ystart'][rank]-1, dom['xstart'][rank]-1]
        lat[1] = latgrd_glob[dom['ystart'][rank]-1, dom['xend']  [rank]-1]
        lat[2] = latgrd_glob[dom['yend']  [rank]-1, dom['xend']  [rank]-1]
        lat[3] = latgrd_glob[dom['yend']  [rank]-1, dom['xstart'][rank]-1]
        isin = linquad(lonobs, latobs, lon, lat)
        if (isin):
            rankout = rank
            break
    return rankout

def get_cell(longrd_loc,latgrd_loc,lonobs, latobs):
    lon = np.ones((4),dtype=np.float64) * -999
    lat = np.ones((4),dtype=np.float64) * -999
    igrd_loc_py = -999
    jgrd_loc_py = -999
    jmt,imt = longrd_loc.shape
    for j in range(jmt-1):
        for i in range(imt-1):
            lon[0] = longrd_loc[j  , i  ]
            lon[1] = longrd_loc[j  , i+1]
            lon[2] = longrd_loc[j+1, i+1]
            lon[3] = longrd_loc[j+1, i  ]
            lat[0] = latgrd_loc[j  , i  ]
            lat[1] = latgrd_loc[j  , i+1]
            lat[2] = latgrd_loc[j+1, i+1]
            lat[3] = latgrd_loc[j+1, i  ]
            isin = linquad(lonobs, latobs, lon, lat)
            if (isin):
                igrd_loc_py = i
                jgrd_loc_py = j
    #FIND LON,LAT OF THE GRID 
    lon[0] = longrd_loc[jgrd_loc_py  , igrd_loc_py  ]
    lon[1] = longrd_loc[jgrd_loc_py  , igrd_loc_py+1]
    lon[2] = longrd_loc[jgrd_loc_py+1, igrd_loc_py+1]
    lon[3] = longrd_loc[jgrd_loc_py+1, igrd_loc_py  ]
    lat[0] = latgrd_loc[jgrd_loc_py  , igrd_loc_py  ]
    lat[1] = latgrd_loc[jgrd_loc_py  , igrd_loc_py+1]
    lat[2] = latgrd_loc[jgrd_loc_py+1, igrd_loc_py+1]
    lat[3] = latgrd_loc[jgrd_loc_py+1, igrd_loc_py  ]
    return igrd_loc_py,jgrd_loc_py,lon,lat

def sel_idx(varin,idx):
    varout = varin.copy()
    list_keys = [x for x in list(varin.keys()) if ((x != "global") and (x != "tpa_correction"))]
    for key in list_keys:
        varout[key]['val']=varin[key]['val'][idx]
    return varout

def dep_to_p(p_dep, p_lat):
  """
  Converts depth to pressure.

  Args:
    p_dep: Depth in meters
    p_lat: Latitude in degrees

  Returns:
    Pressure in Pascals
  """

  z_x = np.sin(p_lat / 57.29578) ** 2
  z_c1 = (5.92 + 5.25 * z_x) * 1e-3
  z_c2 = 2.21e-6
  z_d = (z_c1 - 1) ** 2 - 4 * z_c2 * p_dep
  dep_to_p = ((1 - z_c1) - np.sqrt(z_d)) / (2 * z_c2)

  return dep_to_p

def potemp(ps, pt, pp, ppr):
  """
  Calculates the potential temperature.

  """

  a1 = 1.067610e-05
  a2 = -1.434297e-06
  a3 = -7.566349e-09
  a4 = -8.535585e-06
  a5 = 3.074672e-08
  a6 = 1.918639e-08
  a7 = 1.788718e-10

  zpol = a1 + a2 * ps + a3 * (pp + ppr) + a4 * pt + a5 * ps * pt + a6 * pt * pt + a7 * pt * (pp + ppr)
  potemp = pt + (pp - ppr) * zpol

  return potemp


def linquad(px, py, pxv, pyv):
  """
  Determines if a point (px, py) is within a quadrilateral defined by four points (pxv, pyv).

  Args:
    px: Longitude of the point.
    py: Latitude of the point.
    pxv: Array of longitudes of the quadrilateral's vertices.
    pyv: Array of latitudes of the quadrilateral's vertices.

  Returns:
    True if the point is inside the quadrilateral, False otherwise.
  """

  zst1 = (px - pxv[0]) * (py - pyv[3]) - (py - pyv[0]) * (px - pxv[3])
  if zst1 <= 0.0:
    zst2 = (px - pxv[3]) * (py - pyv[2]) - (py - pyv[3]) * (px - pxv[2])
    if zst2 <= 0.0:
      zst3 = (px - pxv[2]) * (py - pyv[1]) - (py - pyv[2]) * (px - pxv[1])
      if zst3 <= 0.0:
        zst4 = (px - pxv[1]) * (py - pyv[0]) - (py - pyv[1]) * (px - pxv[0])
        if zst4 <= 0.0:
          return True
  return False    


def haversine(lo1,la1,lo2,la2):
#    aversine formula to compute distance (km)
#   between points in latlon coordinates

    za1 = 6378.1370
    zb1 = 6356.7523
    deg2rad =  ( 2 * np.arcsin(1) ) / 180
    lo1 = lo1 * deg2rad
    la1 = la1 * deg2rad
    lo2 = lo2 * deg2rad
    la2 = la2 * deg2rad
    
    dlon=abs(lo1-lo2)
    dlat=abs(la1-la2)

    za = np.sin( dlat/2 )**2 + np.cos(la1)*np.cos(la2)* \
         np.sin(dlon/2)**2
    zc = 2 * np.arctan2(np.sqrt(za),np.sqrt(1-za))
    zr = np.sqrt( za1**4*np.cos(la1)**2 + zb1**4*np.sin(la1)**2 ) / \
         np.sqrt( za1**2*np.cos(la1)**2 + zb1**2*np.sin(la1)**2 )
    out = zc*zr
    return out

def find_sw(array,value):
    dif=(value-array)
    valid_idx = np.where(dif >= 0)[0]
    out = valid_idx[dif[valid_idx].argmin()]
    return out

def find_nearest(array, value):
    #it needs numpy module
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx, array[idx]


def readnc(filename):
    out={}
    nc=Dataset(filename)
    var_list = list(nc.variables.keys())
    for var in var_list:
        out[var]={}
        out[var]['val'] = nc.variables[var][:].squeeze()
        for attr in nc.variables[var].ncattrs():
            out[var][attr] = nc.variables[var].getncattr(attr)
    out['global']={}
    for attr in nc.ncattrs():
        out['global'][attr] = nc.getncattr(attr)
    nc.close()
    return out

def writenc(filename,ins):
    obj_list = list(ins.__dict__.keys())
    obj_list_c = [x.upper() for x in obj_list ]
    nc=Dataset(filename,mode='w',format='NETCDF4')
    for var in obj_list:
        varu = var.upper()
        # Define Dimensions
        if  isinstance(eval('ins.'+var),int):
            if varu=='PREDS':
                dimdim = nc.createDimension(varu, None)
            else:
                dimdim = nc.createDimension(varu, eval('ins.'+var))
        # Define variable
        else:
            if varu=='PQ':
                varvar = nc.createVariable(varu, np.double, ('NPQ2','OBS'),zlib=True)
                varvar[:,:] = eval('ins.'+var) 
            elif varu=='MOI' or varu=='MOJ' or varu=='IB' or varu=='JB':
                varvar = nc.createVariable(varu, np.int32, ('NPQ','OBS'),zlib=True)
                varvar[:,:] = eval('ins.'+var) 
            elif varu=='PLNO':
                varvar = nc.createVariable(varu, 'S1', ('OBS','CHARLEN1'))
                varvar[:,0:7] = eval('ins.'+var) 
            else:
                # if isinstance(eval('ins.'+var+'[0]'),np.int32) or isinstance(eval('ins.'+var+'[0]'),np.int64):
                if isinstance(eval('ins.'+var+'[0]'),int)      or isinstance(eval('ins.'+var+'[0]'),np.int32) \
                or isinstance(eval('ins.'+var+'[0]'),np.int64) or isinstance(eval('ins.'+var+'[0]'),np.int16):
                   varvar = nc.createVariable(varu, np.int32, ('OBS'),zlib=True)
                   varvar[:] = eval('ins.'+var) 
                elif isinstance(eval('ins.'+var+'[0]'),float) or isinstance(eval('ins.'+var+'[0]'),np.float32):
                    varvar = nc.createVariable(varu, np.double, ('OBS'),zlib=True)
                    varvar[:] = eval('ins.'+var) 
    nc.close()
    return 

def read_layout(filename):
    str_tbs='rank    ii    ij   jpi   jpj nimpp njmpp mpiwe mpiea mpiso mpino mpisw mpise mpinw mpine'
    found = False
    first = False
    out={}
    f = open(filename, 'r')
    for line in f:
        found = bool(re.search(str_tbs, line))
        columns = line.split()
        if ( found ):   
            for item in columns:
                out[item] = []
            first = True
            continue
        if (first):
            try:
                for key,item in zip(out.keys(),columns):
                    out[key].append(int(item))
            except ValueError:
                break
    for key in out.keys():
        out[key] = np.array(out[key])
    f.close()
    return out

def defdom(layout,nn_halo):
    out = {}
    layout['mpiwe'][layout['mpiwe']>=0]=1 #west
    layout['mpiea'][layout['mpiea']>=0]=1 #east
    layout['mpiso'][layout['mpiso']>=0]=1 #south
    layout['mpino'][layout['mpino']>=0]=1 #north
    
    layout['mpiwe'][layout['mpiwe']<0]=0 #west
    layout['mpiea'][layout['mpiea']<0]=0 #east
    layout['mpiso'][layout['mpiso']<0]=0 #south
    layout['mpino'][layout['mpino']<0]=0 #north   
    dimx = layout['jpi'] - 2 * nn_halo #- layout['mpiea'] - layout['mpiwe']
    dimy = layout['jpj'] - 2 * nn_halo #- layout['mpino'] - layout['mpiso']
    out['xstart'] = layout['nimpp'] #+ layout['mpiwe']
    out['xend'] = out['xstart'] + dimx - 1
    out['ystart'] = layout['njmpp'] #+ layout['mpiso']
    out['yend'] = out['ystart'] + dimy - 1
    return out

def get_kty(invar):
    outvar = {401:1, 741:2, 831:3, 820:4,"KSAIR2":5}
    return outvar[invar]

def get_sat(platform):
    
    outvar = { "ERS1" :         1,
               "ERS2" :         2,
               "ENVISAT" :      3,
               "GFO" :          4,
               "JASON1" :       5,
               "JASON2" :       6,
               "TP" :           7,
               "CRYOSAT2" :     8,
               "GEOSAT" :       9,
               "ALTIKA" :      10,
               "JASON2N" :     11,
               "JASON3" :      12,
               "SENTINEL3A" :  13,
               "JASON2G" :     14,
               "SENTINEL3B" :  15,
               "HY-2A" :       16,
               "HY-2B" :       17,
               "SENTINEL6A" :  18,
             }
    return outvar[platform]

def get_familycode(family_code):
    KKXBT    = 401        # XBT/MBT
    KKTESAC  = 741        # TESAC
    KKARGO   = 831        # ARGO
    KKBUOY   = 820        # BUOYS
    
    outvar = {"PF": KKARGO,     
           "BA": KKXBT,   
           "XB": KKXBT,   
           "OC": KKXBT,   
           "TE": KKTESAC, 
           "CT": KKTESAC, 
           "HF": KKBUOY,  
           "MO": KKBUOY,  
           "GL": KKBUOY
           }
    return outvar[family_code]


class insitu:
    def __init__(self):
        # DIMENSIONS
        self.obs      = 0;
        self.levs     = 0;
        self.preds    = 0;
        self.npq      = 4;
        self.npq2     = 8;
        self.charlen1 = 8;
        self.charlen2 = 16;
        # VARIABLES
        self.val      = [];     # value
        self.par      = [];     # parameter for T/S
        self.lon      = [];     # longitude
        self.lat      = [];     # latitude
        self.tim      = [];     # time
        self.dpt      = [];     # depth
        self.otype    = [];     # observational type 
        self.plno     = [];     # wmo_platform_code
        self.prof     = [];     # profile number
        self.flc      = [];     # flag
        self.nind     = [];     # index of obs  numeri consecutivi
        self.ino      = [];     # ???           tutti 0
        self.eve      = [];     # eventi        tutti 0
        self.inst     = [];     # instrument    tutti 0
        self.tdist    = [];     # time distance from analysis?
        self.ib       = [];     
        self.jb       = [];     
        self.kb       = [];     
        self.bia      = [];     # bias?
        self.sd1      = [];     # ???
        self.sd2      = [];     # ???
        self.moi      = [];     # ??? tutti 0
        self.moj      = [];     # ??? tutti 0
        self.kty      = [];     #???
        self.prind    = [];     #???
        self.rb       = [];
        self.pq       = [];
        
    def list2array(self):
       obj_list = list(self.__dict__.keys())
       for var in obj_list:
           if not isinstance(eval('self.'+var),int):
               if var=='plno':
                   self.plno=np.array(self.plno,dtype='c')
               else:
                   exec('self.'+var+'=np.array(self.'+var+').T')
                   
class SeaLevelAnomaly:
    def __init__(self):
        # DIMENSIONS
        self.obs      = 0;
        self.levs     = 0;
        self.preds    = 0;
        self.npq      = 4;
        self.npq2     = 4;
        self.charlen1 = 8;
        self.charlen2 = 16;
        # VARIABLES
        self.nind     = [];
        self.ino      = [];
        self.eve      = [];
        self.tdist    = [];
        self.lon      = [];
        self.lat      = [];
        self.dpt      = [];
        self.tim      = [];
        self.val      = [];
        self.bia      = [];
        self.ib       = [];
        self.jb       = [];
        self.pq       = [];
        self.track    = [];
        self.ksat     = [];
        self.bot      = [];
        self.flc      = [];
        self.sd1      = [];
        self.sd2      = [];
        self.moi      = [];
        self.moj      = [];

    def list2array(self):
       obj_list = list(self.__dict__.keys())
       for var in obj_list:
           if not isinstance(eval('self.'+var),int):
              exec('self.'+var+'=np.array(self.'+var+').T')

