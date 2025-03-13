#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 09:38:54 2024

@author: marioadani
"""
import sys
sys.path.insert(1, '../UTILS')
from utils import readnc
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import cartopy.feature as cfea
from matplotlib.colors import BoundaryNorm
import datetime as dt
import numpy as np
import  matplotlib as mpl
#--------------------------------
from urllib.request import urlopen
import ssl
try:
   _create_unverified_https_context = ssl._create_unverified_context
except AttributeError:
    # Legacy Python that doesn't verify HTTPS certificates by default
    pass
else:
    # Handle target environment that doesn't support HTTPS verification
    ssl._create_default_https_context = _create_unverified_https_context
#--------------------------------


# Use custom colormap function from Earle
def custom_div_cmap(numcolors, name='custom_div_cmap',
                    mincol='blue', midcol='white', maxcol='red'):
    """ Create a custom diverging colormap with three colors
    
    Default is blue to white to red with 11 colors.  Colors can be specified
    in any way understandable by matplotlib.colors.ColorConverter.to_rgb()
    """

    from matplotlib.colors import LinearSegmentedColormap 
    
    cmap = LinearSegmentedColormap.from_list(name=name, 
                                             colors =[mincol, midcol, maxcol],
                                             N=numcolors)
    return cmap


#------------------------------------------------------------------------------
# LOAD INS DATASET
dirobs='../..//EXP-A/POST/'
sdate=dt.datetime(2021,1,1)
edate=dt.datetime(2021,1,31)
ndate=sdate
lon=[]
lat=[]
otype=[]
while ndate <=edate:
    ndatestr=ndate.strftime('%Y%m%d')
    print(ndatestr)
    filename=dirobs+'/INSMIS_'+ndatestr+'.nc'
    obs=readnc(filename)
    prof_list=list(set(obs['PROF']))
    for pp in prof_list:
        idx=np.where(pp==obs['PROF'])[0][0]
        lon.append(obs['LON'][idx])
        lat.append(obs['LAT'][idx])
        otype.append(obs['OTYPE'][idx])
    ndate=ndate+dt.timedelta(days=1)
lon_ins=np.array(lon)
lat_ins=np.array(lat)
otype=np.array(otype)
idxbt=np.where(otype==401)
idarg=np.where(otype==831)
#------------------------------------------------------------------------------
#LOAD SLA DATASET
dirobs='../../EXP-A/POST/'
sdate=dt.datetime(2021,1,1)
edate=dt.datetime(2021,1,31)
number_track_red=483
number_track_green=723
date_track_red=dt.datetime(2021,8,4)
ndate=sdate
lon=[]
lat=[]
while ndate <=edate:
    ndatestr=ndate.strftime('%Y%m%d')
    print(ndatestr)
    filename=dirobs+'/SLAMIS_'+ndatestr+'.nc'
    obs=readnc(filename)
    if ndate==dt.datetime(2021,8,4):
        idx=obs['TRACK']==number_track_red
        lon_red=obs['LON'][idx]
        lat_red=obs['LAT'][idx]
    for keys in obs.keys():
        try:
            idx = np.argwhere(np.isnan(obs[keys][0]))
            print('keys : '+keys+' Nan found: '+str(len(idx)))
            obs['FLC'][0][idx] = 0
            idx = np.argwhere(abs(obs[keys][0])>1e20)
            print('keys : '+keys+' Value gt: '+str(len(idx)))
            obs['FLC'][0][idx] = 0
        except TypeError:
            pass
    lon.extend(list(obs['LON']))
    lat.extend(list(obs['LAT']))
    ndate=ndate+dt.timedelta(days=1)
lon_sla=np.array(lon)
lat_sla=np.array(lat)
#------------------------------------------------------------------------------
#LOAD BATHY
grd=readnc('../../DATA/grid_medfs831.nc')
grd['topo']=-grd['topo']
grd['topo'][grd['topo']==0]=0.7
#------------------------------------------------------------------------------
# GEOMETRY FOR PROJECTION
minlon=min(grd['lon'].flatten())
minlat=min(grd['lat'].flatten())
maxlon=max(grd['lon'].flatten())
maxlat=max(grd['lat'].flatten())
avelat=np.mean(grd['lat'].flatten())
avelon=np.mean(grd['lon'].flatten())
projection=ccrs.LambertConformal(central_longitude=avelon, central_latitude=avelat)
bbox=[minlon,maxlon,minlat,maxlat]
#------------------------------------------------------------------------------
# COLORS 
#Setup contour intervals and colormaps for PANEL 1
blevels_p1 = [-5000,  -4750, -4500, -4250, -4000, -3750, -3500, -3250, -3000, \
           -2750, -2500, -2250, -2000, -1750, -1500, -1250, -1000, -750, -500, -350 ,-250, -150, 0]
N = len(blevels_p1)-1
cmap2_p1 = custom_div_cmap(N, mincol='DarkBlue', midcol='CornflowerBlue' ,maxcol='w')
cmap2_p1.set_extremes(over='white',under='indigo')
bnorm_p1 = BoundaryNorm(blevels_p1, ncolors=N, clip=False)
#Setup contour intervals and colormaps for PANEL 3
blevels_p3 = np.linspace(0,17,18)
N = len(blevels_p3)-1
cmap2_p3 = mpl.colormaps['tab20c']
bnorm_p3 = BoundaryNorm(blevels_p3, ncolors=N, clip=False)
#------------------------------------------------------------------------------
# FIGURE
plt.rcParams.update({
#    "text.usetex": True,
    "text.usetex": False,
    "font.family": "sans-serif",
    "font.sans-serif": "Helvetica",
})
fig=plt.figure(figsize=(11.69,8.27))
#------------------------------------------------------------------------------
#PANEL 1
ax=fig.add_subplot(2,1,1,projection=projection)

ax._autoscaleXon = False
ax._autoscaleYon = False

gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=2, color='gray', alpha=0.5, linestyle='--',x_inline =False, y_inline =False)
gl.bottom_labels = False
gl.right_labels = False
gl.xlabel_style = {'size': 12, 'color': 'k', 'rotation':0 }
gl.ylabel_style = {'size': 12, 'color': 'k', 'rotation':0 }
ax.set_extent(bbox,ccrs.PlateCarree())
ax.add_feature(cfea.COASTLINE,lw=.5)
#plot bathy
#pc0=ax.contourf(grd['lon'],grd['lat'],grd['topo'],transform=ccrs.PlateCarree(),\
#                levels=blevels_p1,cmap=cmap2_p1,vmin=-5000, vmax=0,norm=bnorm_p1,extend='both')
#plot contour 1000 350 150
#pc1=ax.contour(pc0,transform=ccrs.PlateCarree(),\
#                levels=[-1000,-350,-150],colors=('navy','blue','royalblue'),linewidths=(1,1,1))
pc1=ax.contour(grd['lon'],grd['lat'],grd['topo'],transform=ccrs.PlateCarree(),\
                levels=[-1000,-350,-150],colors=('navy','blue','royalblue'),linewidths=(1,1,1))
# 
ax.set_position([0, 0.55, 1., 0.4])
#plot border of the model domain
#EAST
pc31=plt.plot(grd['lon'][:,0],grd['lat'][:,0],transform=ccrs.PlateCarree(),\
             linewidth=3,linestyle='-',color='k')
#WEST
pc32=plt.plot(grd['lon'][:,-1],grd['lat'][:,-1],transform=ccrs.PlateCarree(),\
             linewidth=3,linestyle='-',color='k')
#SOUTH
pc33=plt.plot(grd['lon'][0,:],grd['lat'][0,:],transform=ccrs.PlateCarree(),\
             linewidth=3,linestyle='-',color='k')
#NORTH
pc34=plt.plot(grd['lon'][-1,:],grd['lat'][-1,:],transform=ccrs.PlateCarree(),\
             linewidth=3,linestyle='-',color='k')
#plot colorbar
#pos=ax.get_position()
#cbax1=fig.add_axes([pos.x0, 0.5, (pos.x1-pos.x0), 0.025])
#fig.colorbar(pc0,cax=cbax1)
#cbarn=fig.colorbar(pc0,cax=cbax1,location='bottom',extend='both',\
#              spacing='uniform', extendfrac=0.02, drawedges=True)
#cbarn.set_ticks(ticks=blevels_p1,labels=['$'+str(-int(blvl))+'$' for blvl in blevels_p1])
#cbarn.set_label(label='$[m]$',labelpad=-15,x=1.01,fontsize=15)
#cbarn.ax.tick_params(labelsize=15,rotation=45)

ax.plot(lon_ins[idarg],lat_ins[idarg],'.',color='lime',markersize=1,label='$argo$',transform=ccrs.PlateCarree())
ax.plot(lon_ins[idxbt],lat_ins[idxbt],'o',color='yellow',markersize=2,label='$xbt$',transform=ccrs.PlateCarree())
frame = ax.legend(loc='upper right',markerscale=5.,fontsize=15)
frame = frame.get_frame()
frame.set_facecolor('lightgray')
#------------------------------------------------------------------------------
#PANEL 2
ax=fig.add_subplot(2,1,2,projection=projection)
ax._autoscaleXon = False
ax._autoscaleYon = False
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=2, color='gray', alpha=0.5, linestyle='--',x_inline =False, y_inline =False)
gl.top_labels = False
gl.right_labels = False
gl.xlabel_style = {'size': 12, 'color': 'k', 'rotation':0}
gl.ylabel_style = {'size': 12, 'color': 'k', 'rotation':0}
ax.set_extent(bbox,ccrs.PlateCarree())
ax.add_feature(cfea.COASTLINE,lw=.5)
# Plot sla data
pc1=ax.plot(lon_sla,lat_sla,'.',color='gray',markersize=0.5,label='sla',transform=ccrs.PlateCarree())
# Plot border
#EAST
pc31=ax.plot(grd['lon'][:,0],grd['lat'][:,0],transform=ccrs.PlateCarree(),\
             linewidth=3,linestyle='-',color='k')
#WEST
pc32=ax.plot(grd['lon'][:,-1],grd['lat'][:,-1],transform=ccrs.PlateCarree(),\
             linewidth=3,linestyle='-',color='k')
#SOUTH
pc33=ax.plot(grd['lon'][0,:],grd['lat'][0,:],transform=ccrs.PlateCarree(),\
             linewidth=3,linestyle='-',color='k')
#NORTH
pc34=ax.plot(grd['lon'][-1,:],grd['lat'][-1,:],transform=ccrs.PlateCarree(),\
             linewidth=3,linestyle='-',color='k')
pc2=ax.contour(grd['lon'],grd['lat'],grd['topo'],transform=ccrs.PlateCarree(),\
                levels=[-1000,-350,-150],colors=('navy','blue','royalblue'),linewidths=(1,1,1))
# Plot red track referred to figure 2
pos=ax.get_position()
ax.set_position([0, 0.02, 1., 0.4])

#------------------------------------------------------------------------------
plt.savefig('./figure01.png')
#plt.show()
plt.close()

