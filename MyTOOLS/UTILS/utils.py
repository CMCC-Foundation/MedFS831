#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 08:21:13 2024

@author: marioadani
"""
from netCDF4 import Dataset, chartostring
import datetime as dt 
from cftime import num2pydate
import numpy as np
from ttide import t_getconsts as tgc
import ttide as tt
from ttide import time as tm
import copy
import scipy

#------------------------------------------------------------------------------
# MIXED FUNCTIONS
def sel_idx(varin,idx):
    varout = varin.copy()
    for key in list(varin.keys()):
        if key=='ib' or key=='jb' or key=='pq' \
        or key=='tb' or key=='sb' or key=='bcp':
            varout[key]=varin[key][:,idx].squeeze()
        else:
            varout[key]=varin[key][idx].squeeze()
    return varout

def qc_mis(varin):
    for key in list(varin.keys()):
        try:
            idx = np.argwhere(np.isnan(varin[key]))
            varin['flc'][idx] = 0
            idx = np.argwhere(abs(varin[key]>1e10))
            varin['flc'][idx] = 0
        except TypeError:
            pass        
    return varin

def find_nearest(array, value):
    #it needs numpy module
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx, array[idx]

def interp2d(x, y, xp, yp, zp):
    
    """
    Bilinearly interpolate over regular 2D grid.

    `xp` and `yp` are 1D arrays defining grid coordinates of sizes :math:`N_x`
    and :math:`N_y` respectively, and `zp` is the 2D array, shape
    :math:`(N_x, N_y)`, containing the gridded data points which are being
    interpolated from. Note that the coordinate grid should be regular, i.e.
    uniform grid spacing. `x` and `y` are either scalars or 1D arrays giving
    the coordinates of the points at which to interpolate. If these are outside
    the boundaries of the coordinate grid, the resulting interpolated values
    are evaluated at the boundary.

    Parameters
    ----------
    x : 1D array or scalar
        x-coordinates of interpolating point(s).
    y : 1D array or scalar
        y-coordinates of interpolating point(s).
    xp : 1D array, shape M
        x-coordinates of data points zp. Note that this should be a *regular*
        grid, i.e. uniform spacing.
    yp : 1D array, shape N
        y-coordinates of data points zp. Note that this should be a *regular*
        grid, i.e. uniform spacing.
    zp : 2D array, shape (M, N)
        Data points on grid from which to interpolate.

    Returns
    -------
    z : 1D array or scalar
        Interpolated values at given point(s).
        
    downloaded from https://github.com/aneeshnaik/interp2d with MIT licence
    """
    # if scalar, turn into array
    scalar = False
    if not isinstance(x, (list, np.ndarray)):
        scalar = True
        x = np.array([x])
        y = np.array([y])

    # grid spacings and sizes
    hx = xp[1] - xp[0]
    hy = yp[1] - yp[0]
    Nx = xp.size
    Ny = yp.size

    # snap beyond-boundary points to boundary
    x[x < xp[0]] = xp[0]
    y[y < yp[0]] = yp[0]
    x[x > xp[-1]] = xp[-1]
    y[y > yp[-1]] = yp[-1]

    # find indices of surrounding points
    i1 = np.floor((x - xp[0]) / hx).astype(int)
    i1[i1 == Nx - 1] = Nx - 2
    j1 = np.floor((y - yp[0]) / hy).astype(int)
    j1[j1 == Ny - 1] = Ny - 2
    i2 = i1 + 1
    j2 = j1 + 1

    # get coords and func at surrounding points
    x1 = xp[i1]
    x2 = xp[i2]
    y1 = yp[j1]
    y2 = yp[j2]
    z11 = zp[i1, j1]
    z21 = zp[i2, j1]
    z12 = zp[i1, j2]
    z22 = zp[i2, j2]

    # interpolate
    t11 = z11 * (x2 - x) * (y2 - y)
    t21 = z21 * (x - x1) * (y2 - y)
    t12 = z12 * (x2 - x) * (y - y1)
    t22 = z22 * (x - x1) * (y - y1)
    z = (t11 + t21 + t12 + t22) / (hx * hy)
    if scalar:
        z = z[0]
    return z

def SeaOverLand(varin,npts):
    varout = copy.deepcopy(varin)
    varout[np.where(varout==0)]=np.nan
    idx_notnan=np.where(~np.isnan(varout))
    jmt,imt=varin.shape
    for n in range(npts):
        dummy=np.zeros((4,jmt,imt))*np.nan
        dummy[0,:,:]=np.roll(varout,-1,axis=0)
        dummy[1,:,:]=np.roll(varout,1,axis=0)
        dummy[2,:,:]=np.roll(varout,-1,axis=1)
        dummy[3,:,:]=np.roll(varout,1,axis=1)
        dummy=np.nanmean(dummy,axis=0)
        dummy[idx_notnan]=varin[idx_notnan]
        varout = dummy
    return varout

#------------------------------------------------------------------------------
# READINIG/WRITING GENERAL FUNCIONS FUNCTIONS
def readnc(filename):
    out={}
    nc=Dataset(filename)
    var_list = list(nc.variables.keys())
    for var in var_list:
        if var == 'time':
            dummy1 =  nc.variables[var]
            dummy2 =  nc.variables[var][:].compressed()
            out[var] = np.array(num2pydate(dummy2, units=dummy1.units,calendar=dummy1.calendar)) 
        else:
            out[var] = nc.variables[var][:]
    nc.close()
    return out 

def writenc(varin,filename):
    nc=Dataset(filename,'w')
    for keys in varin.keys():
        for d in range(1,varin[keys].ndim+1):
            dummy = nc.createDimension('dim'+str(d)+'_'+keys, varin[keys].shape[d-1])
    for keys in  varin.keys():
        if varin[keys].ndim==1:
            dummy = nc.createVariable(keys, varin[keys].dtype, ('dim1_'+keys),zlib=True)
        elif varin[keys].ndim==2:
            dummy = nc.createVariable(keys, varin[keys].dtype, ('dim1_'+keys,'dim2_'+keys),zlib=True)
        elif  varin[keys].ndim==3:
            dummy = nc.createVariable(keys, varin[keys].dtype, ('dim1_'+keys,'dim2_'+keys,\
                                                                'dim3_'+keys),zlib=True)
        elif  varin[keys].ndim==4:
            dummy = nc.createVariable(keys, varin[keys].dtype, ('dim1_'+keys,'dim2_'+keys,\
                                                                'dim3_'+keys,'dim4_'+keys),zlib=True)
        dummy[:] =  varin[keys]
    nc.close()
    return

def readobs(filename):
    argo = {}
    sla  = {}
    fId = open(filename, "r")
    #Number of sla obs
    dummy=np.fromfile(fId,np.int32,count=1)
    sla['nobs']=np.fromfile(fId,np.int64,count=1)[0]
    dummy=np.fromfile(fId,np.int32,count=1)
    #sla obs
    dummy=np.fromfile(fId,np.int32,count=1)
    sla['ino']=np.fromfile(fId,np.int64,count=sla['nobs'])
    sla['flg']=np.fromfile(fId,np.int64,count=sla['nobs'])
    sla['lon']=np.fromfile(fId,np.float64,count=sla['nobs'])
    sla['lat']=np.fromfile(fId,np.float64,count=sla['nobs'])
    sla['tim']=np.fromfile(fId,np.float64,count=sla['nobs'])
    sla['val']=np.fromfile(fId,np.float64,count=sla['nobs'])
    sla['bac']=np.fromfile(fId,np.float64,count=sla['nobs'])
    sla['err']=np.fromfile(fId,np.float64,count=sla['nobs'])
    sla['res']=np.fromfile(fId,np.float64,count=sla['nobs'])
    sla['bia']=np.fromfile(fId,np.float64,count=sla['nobs'])
    sla['inc']=np.fromfile(fId,np.float64,count=sla['nobs'])
    sla['b_a']=np.fromfile(fId,np.float64,count=sla['nobs'])
    sla['dpt']=np.fromfile(fId,np.float64,count=sla['nobs'])
    sla['dtm']=np.fromfile(fId,np.float64,count=sla['nobs'])
    dummy=np.fromfile(fId,np.int32,count=1)
    
    #Number of arg obs
    dummy=np.fromfile(fId,np.int32,count=1)
    argo['nobs']=np.fromfile(fId,np.int64,count=1)[0]
    dummy=np.fromfile(fId,np.int32,count=1)
    
    #arg obs
    dummy=np.fromfile(fId,np.int32,count=1)
    argo['ino']=np.fromfile(fId,np.int64,count=argo['nobs'])
    argo['flg']=np.fromfile(fId,np.int64,count=argo['nobs'])
    argo['par']=np.fromfile(fId,np.int64,count=argo['nobs'])
    argo['lon']=np.fromfile(fId,np.float64,count=argo['nobs'])
    argo['lat']=np.fromfile(fId,np.float64,count=argo['nobs'])
    argo['dpt']=np.fromfile(fId,np.float64,count=argo['nobs'])
    argo['tim']=np.fromfile(fId,np.float64,count=argo['nobs'])
    argo['val']=np.fromfile(fId,np.float64,count=argo['nobs'])
    argo['bac']=np.fromfile(fId,np.float64,count=argo['nobs'])
    argo['err']=np.fromfile(fId,np.float64,count=argo['nobs'])
    argo['res']=np.fromfile(fId,np.float64,count=argo['nobs'])
    argo['bia']=np.fromfile(fId,np.float64,count=argo['nobs'])
    argo['inc']=np.fromfile(fId,np.float64,count=argo['nobs'])
    argo['b_a']=np.fromfile(fId,np.float64,count=argo['nobs'])
    dummy=np.fromfile(fId,np.int32,count=1)
    fId.close()
    return sla, argo        


def read_iterate(filein):

    out={}
    out['it']=[]
    out['nf']=[]
    out['nint']=[]
    out['nact']=[]
    out['sub']=[]
    out['itls']=[]
    out['stepl']=[]
    out['tstep']=[]
    out['projg']=[]
    out['fun']=[]
    
    with open(filein) as f:
        lines_after_20 = f.readlines()[20:]
        for nl in range(len(lines_after_20)):
            line=lines_after_20[nl].split()
            line = [line.replace('D', 'E') for line in line]
            out['it'].append(np.float64(line[0]))
            out['nf'].append(np.float64(line[1]))
            try:
                out['nint'].append(np.float64(line[2]))
            except ValueError:
                out['nint'].append(np.nan)
            try:   
                out['nact'].append(np.float64(line[3]))
            except ValueError:
                out['nact'].append(np.nan)
            out['sub'].append(line[4])
            try:
                out['itls'].append(np.float64(line[5]))
            except ValueError:
                out['itls'].append(np.nan)
            try:
                out['stepl'].append(np.float64(line[6]))
            except ValueError:
                 out['stepl'].append(np.nan)
            try:           
                out['tstep'].append(np.float64(line[7]))
            except ValueError:
                 out['tstep'].append(np.nan) 
            out['projg'].append(np.float64(line[8]))
            out['fun'].append(np.float64(line[9]))
    return out
#------------------------------------------------------------------------------
#TIDES && SLA
def read_tides(filename):
    nc=Dataset(filename)
    out={}
    out['M2_Amp']=nc.variables['M2_Amp'][:]
    out['K1_Amp']=nc.variables['K1_Amp'][:]
    out['O1_Amp']=nc.variables['O1_Amp'][:]
    out['S2_Amp']=nc.variables['S2_Amp'][:]
    out['P1_Amp']=nc.variables['P1_Amp'][:]
    out['N2_Amp']=nc.variables['N2_Amp'][:]
    out['Q1_Amp']=nc.variables['Q1_Amp'][:]
    out['K2_Amp']=nc.variables['K2_Amp'][:]
    
    out['M2_Pha']=nc.variables['M2_Pha'][:]
    out['K1_Pha']=nc.variables['K1_Pha'][:]
    out['O1_Pha']=nc.variables['O1_Pha'][:]
    out['S2_Pha']=nc.variables['S2_Pha'][:]
    out['P1_Pha']=nc.variables['P1_Pha'][:]
    out['N2_Pha']=nc.variables['N2_Pha'][:]
    out['Q1_Pha']=nc.variables['Q1_Pha'][:]
    out['K2_Pha']=nc.variables['K2_Pha'][:]
    try :
      msk=nc.variables['sossheig'][0,:,:].squeeze()
      msk[~msk.mask]=1
      out['msk']=msk
    except KeyError:
      out['msk']=nc.variables['tmask'][:,:]
    out['grid_longitudes']=nc.variables['lon'][:]
    out['grid_latitudes']=nc.variables['lat'][:]

    nc.close()
    return out

def setup_tides(filetides):
    # Read tides constituents
    tide   = read_tides(filetides)
    # Model tides
    const, sat, shallow = tgc.t_getconsts(np.empty(0))
    med_tide=tt.TTideCon()
    med_tide['nameu']=np.array(['M2  ','K1  ','O1  ','S2  ',\
                                'P1  ','N2  ','Q1  ','K2  '],dtype='S4') 
    nconst=len(med_tide['nameu'])
    med_tide['fu']=np.zeros(nconst)    
    for i in range(0,nconst):
        idx=np.where(med_tide['nameu'][i]==const['name'])[0][0]
        med_tide['fu'][i]=const['freq'][idx]
    med_tide['ltype']='nodal'
    med_tide['synth']=0
    med_tide['ref_time']=tm.datetime(1950,1,1)
    for keys in tide.keys():
        med_tide[keys]=tide[keys]
    return med_tide

def shuffle_obs(sla):
    sla1=copy.deepcopy(sla)
    slano=len(sla['flc'])
    # for keys in sla.keys():
    #     sla1[keys]=np.zeros((slano))
    for i in range(1,slano+1):
        idx=np.where(i==sla['nind'])[0][0]
        for keys in sla.keys():
            if sla[keys].size == slano:
                sla1[keys][i-1]=sla[keys][idx]
            else:
                sla1[keys][:,i-1]=sla[keys][:,idx]
    return sla1

def unbias(sla,dsm):
    sla['biao'] = np.zeros(len(sla['res']))
    sla['biam'] = np.zeros(len(sla['res']))
    minobspt=4
    for iter in range(1):
        #bias
        satino=sla['track'][0]
        satk=sla['ksat'][0]
        #ADD new not in 3dVAR
        # sattim=sla['tim'][0]
        # dsm = 100.
        i1 = 0
        slano=len(sla['flc'])
        for k in range(1,slano):
          dxx = 6371.0*3.14/180.0 * (sla['lon'][k]-sla['lon'][k-1]) * np.cos(sla['lat'][k]*3.14/180.0)
          dyy = 6371.0*3.14/180.0 * (sla['lat'][k]-sla['lat'][k-1])
          if( ((sla['track'][k] != satino) or (np.sqrt(dxx**2+dyy**2)>dsm) )  & (k != slano-1) ):
            sumt = 0.0
            sumto = 0.0
            sumi = 0.0
            sumtm = 0.0
            for i in range(i1,k):
                if(sla['flc'][i] == 1):
                    sumt =  sumt + sla['res'][i]
                    sumto = sumto + sla['val'][i]
                    sumtm = sumtm + sla['bac'][i]
                    sumi =  sumi + 1.0
            if (sumi > 0.0):
                sumt = sumt/sumi
                sumto = sumto/sumi
                sumtm = sumtm/sumi
            if (sumi > minobspt):
                for i in range(i1,k):
                    sla['res'][i] = sla['res'][i] - sumt
                    sla['val'][i] = sla['val'][i] - sumto
                    sla['bac'][i] = sla['bac'][i] - sumtm
                    sla['biao'][i] = sumto
                    sla['biam'][i] = sumtm
                    sla['bia'][i] = sumt  
            else: # Track too short
                for i in range(i1,k):
                    sla['res'][i] = sla['res'][i] - sumt
                    sla['val'][i] = sla['val'][i] - sumto
                    sla['bac'][i] = sla['bac'][i] - sumtm
                    sla['biao'][i] = sumto
                    sla['biam'][i] = sumtm
                    sla['bia'][i] = sumt  
                    sla['flc'][i] = 0
            satino=sla['track'][k]
            satk=sla['ksat'][k]
            i1 = k
          elif (k == slano-1):         
            sumt = 0.0 
            sumtm = 0.0
            sumto = 0.0                                 
            sumi = 0.0                                  
            for i in range(i1,k+1):                            
                if(sla['flc'][i]==1):             
                    sumt = sumt + sla['res'][i]    
                    sumto = sumto + sla['val'][i]
                    sumtm = sumtm + sla['bac'][i]
                    sumi = sumi + 1.0               
            if(sumi>0.0): 
                sumt = sumt/sumi  
                sumto = sumto/sumi
                sumtm = sumtm/sumi
            if (sumi > minobspt):
                for i in range (i1,k+1):                           
                    sla['res'][i] = sla['res'][i] - sumt
                    sla['val'][i] = sla['val'][i] - sumto    
                    sla['bac'][i] = sla['bac'][i] - sumtm
                    sla['biao'][i] = sumto
                    sla['biam'][i] = sumtm
                    sla['bia'][i] = sumt
            else: # Track too short
                for i in range(i1,k+1):
                    sla['res'][i] = sla['res'][i] - sumt
                    sla['val'][i] = sla['val'][i] - sumto
                    sla['bac'][i] = sla['bac'][i] - sumtm
                    sla['biao'][i] = sumto
                    sla['biam'][i] = sumtm
                    sla['bia'][i] = sumt  
                    sla['flc'][i] = 0                  
    return sla
#------------------------------------------------------------------------------
#STATISTICS
class statistics:
    def __init__(self):
        self.list =['bias','rmse','anb','mnb','mae','mne','corr','stdo','stdm',\
                    'crmse','nmsd','mfb','mfe','fb','skvar','cv','nmse','fac2',\
                    'stdoa','ioa','obs','mod','nobs']
        self.bias = []
        self.rmse = []
        self.anb = []
        self.mnb = []
        self.mae = []
        self.mne=[]
        self.corr=[]
        self.stdo=[]
        self.stdm=[]
        self.stdoa=[]
        self.crmse=[]
        self.nmsd=[]
        self.mfb=[]
        self.mfe=[]
        self.fb=[]
        self.skvar=[]
        self.cv=[]
        self.nmse=[]
        self.fac2=[]
        self.ioa=[]
        self.obs=[]
        self.mod=[]
        self.nobs=[]
        self.time=[]
        self.level=[]
        self.param=[]

    def compute(self,obs,mod):
        self.bias.append(mk_bias(obs,mod,1))
        self.rmse.append(mk_rmse(obs,mod,1))
        self.anb.append(mk_anb(obs,mod,1)) 
        self.mnb.append(mk_mnb(obs,mod,1))
        self.mae.append(mk_mae(obs,mod,1))
        self.mne.append(mk_mne(obs,mod,1))
        self.corr.append(mk_corr(obs,mod,2))
        self.stdm.append(mk_std(obs,1))
        self.stdo.append(mk_std(mod,1))
        self.stdoa.append(mk_std(obs-mod,1))
        self.crmse.append(mk_crmse(obs,mod,1))
        self.nmsd.append(mk_nmsd(obs,mod,1))
        self.mfb.append(mk_mfb(obs,mod,1))
        self.mfe.append(mk_mfe(obs,mod,1))
        self.fb.append(mk_fb(obs,mod,1))
        self.skvar.append(mk_skvar(obs,mod,1))
        self.cv.append(mk_cv(obs,mod,1))
        self.nmse.append(mk_nmse(obs,mod,1))
        self.fac2.append(mk_fac2(obs,mod,1))
        self.ioa.append(mk_ioa(obs,mod,1))
        self.obs.append(np.mean(obs))
        self.mod.append(np.mean(mod))
        self.nobs.append(len(mod))
    
    def to_array(self):
        for ll in self.list:
            exec('self.'+ll+'=np.array(self.'+ll+')')
        

def mk_bias(x,y,nval):
#BIAS
    if x.size < nval:
       out=np.nan
    else:
       out=(x-y).mean()
    return out

def mk_rmse(x,y,nval):
#RMSE
    if x.size < nval:
       out=np.nan
    else:
       out=np.sqrt(((x-y)**2).mean())
    return out

def mk_anb(x,y,nval):
#Average normalised absolute Bias
    if x.size < nval:
       out=np.nan
    else:
       XM=x.mean()
       YM=y.mean()
       out=(XM-YM)/YM
    return out

def mk_mnb(x,y,nval):
#Mean Normalised BIAS
    if x.size < nval:
       out=np.nan
    else:
       out=((x-y)/y).mean()
    return out

def mk_mae(x,y,nval):
#Mean Absolute Gross Error
    if x.size < nval:
       out=np.nan
    else:
       out=abs(x-y).mean()
    return out

def mk_mne(x,y,nval):
#Mean Normalised Gross Error
    if x.size < nval:
       out=np.nan
    else:
       out=(abs(x-y)/y).mean()
    return out

def mk_corr(x,y,nval):
#Correlation
    if x.size < nval:
       out=np.nan
    else:
       out=scipy.stats.pearsonr(x.squeeze(), y.squeeze())[0]
    return out
#Alternative way to compute correlation. Save results but slower function
#def mk_corr(x,y):
#Correlation
#    if x.size == 0:
#       out=np.nan
#    else:
#       out=np.corrcoef(x, y)[0][1]
#    return out

def mk_std(x,nval):
#Standard deviation
    if x.size < nval:
       out=np.nan
    else:
       out=np.std(x)
    return out


def mk_crmse(x,y,nval):
#Centered Root Mean Square Error / Standard deviation of error
    if x.size < nval:
       out=np.nan
    else:
       xm=x.mean()
       ym=y.mean()
       out=np.sqrt((((x-xm)-(y-ym))**2).mean())
    return out

def mk_nmsd(x,y,nval):
#Normalised Mean Standard Deviation
    if x.size < nval:
       out=np.nan
    else:
       xm=np.std(x)
       ym=np.std(y)
       out=(xm-ym)/ym
    return out

def mk_mfb(x,y,nval):
#Mean Fractional Bias
    if x.size < nval:
       out=np.nan
    else:
       out=((x-y)/((x+y)*0.5)).mean()
    return out

def mk_mfe(x,y,nval):
#Mean Fractional Error
    if x.size < nval:
       out=np.nan
    else:
       out=(abs(x-y)/((x+y)*0.5)).mean()
    return out

def mk_fb(x,y,nval):
#Fractional Bias
    if x.size < nval:
       out=np.nan
    else:
       xm=x.mean()
       ym=y.mean()
       out=(xm-ym)/(0.5*(xm+ym))
    return out

def mk_skvar(x,y,nval):
#Skill Variance
    if x.size < nval:
       out=np.nan
    else:
       xm=np.std(x)
       ym=np.std(y)
       out=xm/ym
    return out

def mk_cv(x,y,nval):
#Coefficient of Variation
    if x.size < nval:
       out=np.nan
    else:
       xm=x.mean()
       ym=y.mean()
       out=np.sqrt((((x-xm)-(y-ym))**2).mean())/ym
    return out

def mk_nmse(x,y,nval):
#Normalised Mean Square Error
    if x.size < nval:
       out=np.nan
    else:
       xm=x.mean()
       ym=y.mean()
       out=((x-y)**2).mean()/(xm*ym)
    return out

def mk_fac2(x,y,nval):
#Fraction of Predictions within a factor of two of observations
    if x.size < nval:
       out=np.nan
    else:
       a=abs(x/y)
       index1=np.where(np.logical_and(a>=0.5, a<=2))
       index0=np.where(np.logical_or(a<0.5, a>2))
       a[index1]=1
       a[index0]=0
       out=a.mean()
    return out

def mk_ioa(x,y,nval):
#Index of Agreement
    if x.size < nval:
       out=np.nan
    else:
       ym=y.mean()
       out=1-( sum((x-y)**2)/sum((abs(x-ym)+abs(y-ym))**2) )
    return out

#-----------------------------------------------------------------------
class insitu:
    def __init__(self):
        self.list=['bac','val','tim','lon','lat','dpt','prof','par','reg','incr']
        for ll in self.list:
            exec('self.'+ll+'=[]')
        return
    
    def to_array(self):
        for ll in self.list:
            exec('self.'+ll+'=np.array(self.'+ll+',dtype=np.float32)')
        return

    def to_dict(self):
        self1={}
        for ll in self.list:
            exec('self1["'+ll+'"]=self.'+ll)
        return self1

#-----------------------------------------------------------------------
class obsstat_insitu:
    def __init__(self):
        self.list=['bkgT','valT','incT','resT',\
                   'bkgS','valS','incS','resS',\
                   'tim','lon','lat','depth']
        for ll in self.list:
            exec('self.'+ll+'=[]')
        return
    
    def to_array(self):
        for ll in self.list:
            exec('self.'+ll+'=np.array(self.'+ll+',dtype=np.float32)')
        return

    def to_dict(self):
        self1={}
        for ll in self.list:
            exec('self1["'+ll+'"]=self.'+ll)
        return self1
