#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 08:07:22 2024

@author: marioadani
"""

import sys
sys.path.insert(1, '../UTILS')
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from utils import readnc, statistics,find_nearest,sel_idx
import matplotlib.colors               as cl
from cftime import num2pydate
import datetime as dt


exp_list=['EXP-S','EXP-A']
st_list=['corr','bias','crmse','stdo','stdm']
rootdir = '../../DATA/'
grd=readnc('../../DATA/grid_medfs831.nc')
grd['dep']=grd['dep'][0::1]
max_lvl,dummy=find_nearest(grd['dep'], 2100.)
max_lvl=max_lvl+1
grd['dep']=grd['dep'][0:max_lvl]
dpt_list=[int(-x) for x in list(grd['dep'])] 
#Find the mean of the depth where the statistics are computed 
dpt_list = [ int(sum(dpt_list[x : x + 2]) / 2) for x in range(0, len(dpt_list)-1)]
arc_exp={}

for exp_name in exp_list:
    if  exp_name!='EXP-S':
        arc_exp[exp_name]={}
        arc_exp[exp_name+'A']={}
    else:
        arc_exp[exp_name]={}
    for st in st_list:
        for param in range(1,3):
            if exp_name!='EXP-S':
                exec('arc_exp[exp_name]["'+st+str(param)+'"]=np.zeros((len(dpt_list)))')
                exec('arc_exp[exp_name+"A"]["'+st+str(param)+'"]=np.zeros((len(dpt_list)))')
            else:
                exec('arc_exp[exp_name]["'+st+str(param)+'"]=np.zeros((len(dpt_list)))')
                

for exp_name in exp_list:
    print('------------'+exp_name+'------------')
    filename=rootdir+'/ins_'+exp_name+'.nc'
    ins=readnc(filename)
    # cycle in dates parameters and depth
    for param in range(1,3):
            # for lvl in range(grd['dep'].shape[0]-1):
            for lvl in range(len(dpt_list)):
                idx=np.where( (ins['dpt']>=grd['dep'][lvl]) & (ins['dpt']<grd['dep'][lvl+1]) &\
                              (ins['par']==param          ) )
                stat=statistics()
                stat.compute(ins['val'][idx], ins['bac'][idx])
                if exp_name!='EXP-S':
                    stata=statistics()
                    stata.compute(ins['val'][idx], ins['bac'][idx]+ins['incr'][idx])
                for st in st_list:
                    if exp_name!='EXP-S':
                        exec('arc_exp[exp_name][st+str(param)][lvl]=stat.'+st+'[0]')
                        exec('arc_exp[exp_name+"A"][st+str(param)][lvl]=stata.'+st+'[0]')
                    else:
                        exec('arc_exp[exp_name][st+str(param)][lvl]=stat.'+st+'[0]')

depthmin=dpt_list[0]
depthmax=dpt_list[-1]
depthud=-400  
dpt_list=abs(np.array(dpt_list))
height=0.2
width=0.2
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
ax1=fig.add_subplot(4,4,1)
ax1.plot(arc_exp['EXP-S']['crmse1'],dpt_list,color='blue')
ax1.plot(arc_exp['EXP-A']['crmse1'],dpt_list,color='g')
ax1.set_yticks(np.linspace(0,500,6)[:-1])
ax1.grid()
ax1.set_ylim([0,500])
ax1.set_xlim([0,0.35])
ax1.invert_yaxis()
ax1.set_xticklabels('')
ax1.tick_params(labelsize=15)
ax1.set_ylabel('$depth\ [m]$',fontsize=15,loc='bottom')
ax1.set_position([0.07,0.75,width,height])
ax1.set_title('$uRMSE$',fontsize=15)

ax2=fig.add_subplot(4,4,5)
ax2.plot(arc_exp['EXP-S']['crmse1'],dpt_list,color='blue',label='$EXP-S$')
ax2.plot(arc_exp['EXP-A']['crmse1'],dpt_list,color='g',label='$EXP{-}A$')
ax2.grid()
ax2.set_xlim([0,0.35])
ax2.set_ylim([500,dpt_list[-1]])
ax2.invert_yaxis()
ax2.tick_params(labelsize=15)
ax2.set_xlabel('$[psu]$',fontsize=15)
ax2.legend(loc='lower right')
ax2.set_position([0.07,0.55,width,height])
#------------------------------------------------------------------------------
ax1=fig.add_subplot(4,4,2)
ax1.plot(-arc_exp['EXP-S']['bias1'],dpt_list,color='blue')
ax1.plot(-arc_exp['EXP-A']['bias1'],dpt_list,color='g')
# ax1.plot(-arc_exp['EXP-AA']['bias1'],dpt_list,color='g',linestyle='--')
ax1.set_yticks(np.linspace(0,500,6)[:-1])
ax1.grid()
ax1.set_ylim([0,500])
ax1.set_xlim([-0.15,0.15])
ax1.invert_yaxis()
ax1.set_xticklabels('')
ax1.set_yticklabels('')
ax1.set_position([0.28,0.75,width,height])
ax1.set_title('$MB$',fontsize=15)

ax2=fig.add_subplot(4,4,6)
ax2.plot(-arc_exp['EXP-S']['bias1'],dpt_list,color='blue')
ax2.plot(-arc_exp['EXP-A']['bias1'],dpt_list,color='g')
# ax2.plot(-arc_exp['EXP-AA']['bias1'],dpt_list,color='g',linestyle='--')
ax2.grid()
ax2.set_yticklabels('')
ax2.set_ylim([500,dpt_list[-1]])
ax2.set_xlim([-0.15,0.15])
ax2.invert_yaxis()
ax2.tick_params(labelsize=15)
ax2.set_xlabel('$[psu]$',fontsize=15)
ax2.set_position([0.28,0.55,width,height])
#------------------------------------------------------------------------------
ax1=fig.add_subplot(4,4,3)
ax1.plot(arc_exp['EXP-S']['corr1'],dpt_list,color='blue',)
ax1.plot(arc_exp['EXP-A']['corr1'],dpt_list,color='g')
ax1.set_yticks(np.linspace(0,500,6)[:-1])
ax1.grid()
ax1.set_ylim([0,500])
ax1.set_xlim([0.95,1])
ax1.invert_yaxis()
ax1.set_xticklabels('')
ax1.set_yticklabels('')
ax1.set_position([0.49,0.75,width,height])
ax1.set_title('$CC$',fontsize=15)

ax2=fig.add_subplot(4,4,7)
ax2.plot(arc_exp['EXP-S']['corr1'],dpt_list,color='blue',)
ax2.plot(arc_exp['EXP-A']['corr1'],dpt_list,color='g')
# ax2.plot(arc_exp['EXP-AA']['corr1'],dpt_list,color='g',linestyle='--')
ax2.grid()
ax2.set_yticklabels('')
ax2.set_ylim([500,dpt_list[-1]])
ax2.set_xlim([0.95,1])
ax2.invert_yaxis()
ax2.tick_params(labelsize=15)
ax2.set_xlabel('$[\#]$',fontsize=15)
ax2.set_position([0.49,0.55,width,height])
#------------------------------------------------------------------------------
ax1=fig.add_subplot(4,4,4)
ax1.plot(arc_exp['EXP-S']['stdm1']-arc_exp['EXP-S']['stdo1'],dpt_list,color='blue')
ax1.plot(arc_exp['EXP-A']['stdm1']-arc_exp['EXP-A']['stdo1'],dpt_list,'g')
ax1.set_yticks(np.linspace(0,500,6)[:-1])
ax1.set_xticks(np.linspace(-0.025,0.025,3))
ax1.grid()
ax1.set_ylim([0,500])
ax1.set_xlim([-0.05,0.05])
ax1.invert_yaxis()
ax1.set_xticklabels('')
ax1.set_yticklabels('')
ax1.set_position([0.70,0.75,width,height])
ax1.set_title('$SDE$',fontsize=15)

ax2=fig.add_subplot(4,4,8)
ax2.plot(arc_exp['EXP-S']['stdm1']-arc_exp['EXP-S']['stdo1'],dpt_list,color='blue')
ax2.plot(arc_exp['EXP-A']['stdm1']-arc_exp['EXP-A']['stdo1'],dpt_list,'g')
ax2.grid()
ax2.set_yticklabels('')
ax2.set_ylim([500,dpt_list[-1]])
ax2.set_xlim([-0.05,0.05])
ax2.set_xticks(np.linspace(-0.025,0.025,3))
ax2.invert_yaxis()
ax2.tick_params(labelsize=15)
ax2.set_xlabel('$[psu]$',fontsize=15)
ax2.set_position([0.70,0.55,width,height])



#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
ax1=fig.add_subplot(4,4,9)
ax1.plot(arc_exp['EXP-S']['crmse2'],dpt_list,color='blue')
ax1.plot(arc_exp['EXP-A']['crmse2'],dpt_list,'g')
ax1.set_yticks(np.linspace(0,500,6)[:-1])
ax1.grid()
ax1.set_ylim([0,500])
ax1.invert_yaxis()
ax1.set_xticklabels('')
ax1.tick_params(labelsize=15)
ax1.set_ylabel('$depth\ [m]$',fontsize=15,loc='bottom')
ax1.set_position([0.07,0.27,width,height])

ax2=fig.add_subplot(4,4,13)
ax2.plot(arc_exp['EXP-S']['crmse2'],dpt_list,color='blue')
ax2.plot(arc_exp['EXP-A']['crmse2'],dpt_list,'g')
# ax2.plot(arc_exp['EXP-AA']['crmse2'],dpt_list,'g',linestyle='--')
ax2.grid()
ax2.set_ylim([500,dpt_list[-1]])
ax2.invert_yaxis()
ax2.tick_params(labelsize=15)
ax2.set_xlabel('$[{^o}C]$',fontsize=15)
ax2.set_position([0.07,0.07,width,height])
#------------------------------------------------------------------------------
ax1=fig.add_subplot(4,4,10)
ax1.plot(-arc_exp['EXP-S']['bias2'],dpt_list,color='blue')
ax1.plot(-arc_exp['EXP-A']['bias2'],dpt_list,'g')
ax1.set_yticks(np.linspace(0,500,6)[:-1])
ax1.grid()
ax1.set_ylim([0,500])
ax1.set_xlim([-0.3,0.3])
ax1.invert_yaxis()
ax1.set_xticklabels('')
ax1.set_yticklabels('')
ax1.set_position([0.28,0.27,width,height])

ax2=fig.add_subplot(4,4,14)
ax2.plot(-arc_exp['EXP-S']['bias2'],dpt_list,color='blue')
ax2.plot(-arc_exp['EXP-A']['bias2'],dpt_list,'g')
ax2.grid()
ax2.set_yticklabels('')
ax2.set_ylim([500,dpt_list[-1]])
ax2.set_xlim([-0.3,0.3])
ax2.invert_yaxis()
ax2.tick_params(labelsize=15)
ax2.set_xlabel('$[{^o}C]$',fontsize=15)
ax2.set_position([0.28,0.07,width,height])
#------------------------------------------------------------------------------
ax1=fig.add_subplot(4,4,11)
ax1.plot(arc_exp['EXP-S']['corr2'],dpt_list,color='blue')
ax1.plot(arc_exp['EXP-A']['corr2'],dpt_list,'g')
ax1.set_yticks(np.linspace(0,500,6)[:-1])
ax1.grid()
ax1.set_ylim([0,500])
ax1.invert_yaxis()
ax1.set_xticklabels('')
ax1.set_yticklabels('')
ax1.set_position([0.49,0.27,width,height])

ax2=fig.add_subplot(4,4,15)
ax2.plot(arc_exp['EXP-S']['corr2'],dpt_list,color='blue')
ax2.plot(arc_exp['EXP-A']['corr2'],dpt_list,'g')
ax2.grid()
ax2.set_yticklabels('')
ax2.set_ylim([500,dpt_list[-1]])
ax2.invert_yaxis()
ax2.tick_params(labelsize=15)
ax2.set_xlabel('$[\#]$',fontsize=15)
ax2.set_position([0.49,0.07,width,height])
#------------------------------------------------------------------------------
ax1=fig.add_subplot(4,4,12)
ax1.plot(arc_exp['EXP-S']['stdm2']-arc_exp['EXP-S']['stdo2'],dpt_list,color='blue')
ax1.plot(arc_exp['EXP-A']['stdm2']-arc_exp['EXP-A']['stdo2'],dpt_list,'g')
ax1.set_yticks(np.linspace(0,500,6)[:-1])
ax1.grid()
ax1.set_ylim([0,500])
ax1.set_xlim([-0.2,0.2])
ax1.set_xticks(np.linspace(-0.1,0.1,3))
ax1.invert_yaxis()
ax1.set_xticklabels('')
ax1.set_yticklabels('')
ax1.set_position([0.70,0.27,width,height])

ax2=fig.add_subplot(4,4,16)
ax2.plot(arc_exp['EXP-S']['stdm2']-arc_exp['EXP-S']['stdo2'],dpt_list,color='blue')
ax2.plot(arc_exp['EXP-A']['stdm2']-arc_exp['EXP-A']['stdo2'],dpt_list,'g')
ax2.grid()
ax2.set_yticklabels('')
ax2.set_ylim([500,dpt_list[-1]])
ax2.set_xlim([-0.2,0.2])
ax2.set_xticks(np.linspace(-0.1,0.1,3))
ax2.invert_yaxis()
ax2.tick_params(labelsize=15)
ax2.set_xlabel('$[{^o}C]$',fontsize=15)
ax2.set_position([0.70,0.07,width,height])
plt.savefig('./figure06.png')







