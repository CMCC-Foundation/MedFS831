#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 26 17:35:30 2024

@author: marioadani
"""

import sys
sys.path.insert(1, '../UTILS')
from utils import readnc,statistics, sel_idx
import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
from cftime import num2pydate
from dateutil.relativedelta import relativedelta
import numpy as np


#------------------------------------------------------------------------------
# 0) User Param
exp_list=['EXP-S','EXP-A']
rootdir = '../../DATA/'
st_list=['crmse','corr','skvar']
dpt_list=np.linspace(125,4000,32)
dpt_list=np.linspace(0,4200,13) # 350mm
dpt_list=np.linspace(0,4000,9)  # 500
dpt_list=np.array([0,150,350,500,1000,1500,2000,2500,3000,3500,4000])
sdate=dt.datetime(2021,1,1)

tstep=1
days=31
numdays=np.int32(days/tstep)+1
date_list=date_list = [sdate + relativedelta(days=tstep*x) for x in range(numdays)]

# 

#------------------------------------------------------------------------------
# 1) Read input dataset
#------------------------------------------------------------------------------
arc_inp={}
for  exp_name in exp_list:
    print('Reading input data: '+exp_name)
    filename=rootdir+'/sla_'+exp_name+'.nc'
    arc_inp[exp_name]=readnc(filename)
    arc_inp[exp_name]['tim']=num2pydate(arc_inp[exp_name]['tim'],units='days since 1950-1-1',calendar='standard')

#------------------------------------------------------------------------------
# 2) Compute statitstic based on time
#------------------------------------------------------------------------------
arc_exp={}    
for  exp_name in exp_list:
    print('Compute statistics based time for input data: '+exp_name)
    stat=statistics()
    for ll in range(len(date_list)-1):
        idx=( (arc_inp[exp_name]['tim']>=date_list[ll]) & (arc_inp[exp_name]['tim']<date_list[ll+1]) )
        stat.compute(arc_inp[exp_name]['val'][idx],arc_inp[exp_name]['bac'][idx])
    stat.to_array()
    arc_exp[exp_name]=stat
          
#------------------------------------------------------------------------------    
# 3) Plot time series
plt.rcParams.update({
#    "text.usetex": True,
    "text.usetex": False,
    "font.family": "sans-serif",
    "font.sans-serif": "Helvetica",
})
fig=plt.figure(figsize=(11.69,8.27))
gs = fig.add_gridspec(3, 1)  
# PANEL 1
ax1 = fig.add_subplot(gs[0, :])
ax1.plot(date_list[1:],arc_exp['EXP-S'].crmse*100,'-',color='b',linewidth=3,\
        label='$EXP-S$')
ax1.plot(date_list[1:],arc_exp['EXP-A'].crmse*100,'-',color='g',linewidth=3,\
        label='$EXP{-}A$')
ax1.set_xticklabels('')
ax1.set_ylabel('$uRMSE\ [cm]$',fontsize=15)
ax1.grid()
# ax1.legend(fontsize=15)
ax1.set_xlim(date_list[0],date_list[-1])
ax1.set_ylim(2,7)
ax1.tick_params(labelsize=15)
ax1.set_position([0.07,0.69,0.88,0.3])

# PANEL 2
ax2 = fig.add_subplot(gs[1, :])
ax2.plot(date_list[1:],arc_exp['EXP-S'].corr,'-',color='b',linewidth=3,\
        label='$EXP-S$')
ax2.plot(date_list[1:],arc_exp['EXP-A'].corr,'-',color='g',linewidth=3,\
        label='$EXP{-}A$')
ax2.set_xticklabels('')
ax2.set_ylabel('$CC\ [\#]$',fontsize=15)
ax2.grid()
# ax2.legend(fontsize=15)
ax2.set_xlim(date_list[0],date_list[-1])
ax2.set_ylim(0,0.9)
ax2.tick_params(labelsize=15)
ax2.set_position([0.07,0.38,0.88,0.3])

# PANEL 3
ax3 = fig.add_subplot(gs[2, :])
ax3.plot(date_list[1:],(arc_exp['EXP-S'].stdm-arc_exp['EXP-S'].stdo)*100,'-',color='b',linewidth=3,\
        label='$EXP-S$')
ax3.plot(date_list[1:],(arc_exp['EXP-A'].stdm-arc_exp['EXP-A'].stdo)*100,'-',color='g',linewidth=3,\
        label='$EXP{-}A$')
# ax3.plot(date_list[1:],(arc_exp['EXP-AA'].stdm-arc_exp['EXP-AA'].stdo)*100,'--',color='g',linewidth=3,\
#         label='$EXP{-}1_{r}$')
# ax1.set_xticks(date_list)
# ax1.set_xticklabels('')
ax3.set_ylabel('$SDE\ [cm]$',fontsize=15)
ax3.set_xlabel('$time$',fontsize=15)
ax3.grid()
ax3.legend(ncol=3,fontsize=15)
ax3.tick_params(labelsize=15)
ax3.set_xlim(date_list[0],date_list[-1])
ax3.set_position([0.07,0.07,0.88,0.3])
plt.savefig('figure04.png')


# #------------------------------------------------------------------------------    
# # 3) Compute ste statitstic based on depth
# #------------------------------------------------------------------------------
# arc_exp={}
# for exp_name in exp_list:
#    print('Compute statistics based on water column depth for input data: '+exp_name)
#    stat=statistics()
#    for ll in range(len(dpt_list)-1):
#        idx=( (arc_inp[exp_name]['dpt']>=dpt_list[ll]) & (arc_inp[exp_name]['dpt']<dpt_list[ll+1]) )
#        stat.compute(arc_inp[exp_name]['val'][idx],arc_inp[exp_name]['bac'][idx])
#    stat.to_array()
#    arc_exp[exp_name]=stat
# 
#------------------------------------------------------------------------------    
# # 3) Plot function of depth series
# size = range(len(dpt_list[0:-1]))
# plt.rcParams.update({
#     "text.usetex": True,
#     "font.family": "sans-serif",
#     "font.sans-serif": "Helvetica",
# })
# fig=plt.figure(figsize=(11.69,8.27))
# gs = fig.add_gridspec(3, 1)  
# # PANEL 1
# ax1 = fig.add_subplot(gs[0, :])
# ax11 = ax1.twinx() # Create another axes that shares the same x-axis as ax.
# p1=ax1.bar([pos  for pos in size],arc_exp['EXP-S'].crmse*100,color='b',edgecolor='black',\
#         width=0.2, align='center', label='$EXP-S$')
# p2=ax1.bar([pos + 0.2 for pos in size],arc_exp['EXP-A'].crmse*100,color='g',edgecolor='black',\
#         width=0.2, align='center',label='$EXP{-}A$')
# ax1.set_xticks([pos - 0.4 for pos in size], [str(int(dpt)) for dpt in dpt_list[0:-1]])
# ax1.grid()
# ax1.tick_params(labelsize=15)
# ax1.set_ylabel('$uRMSE\ [cm]$',fontsize = 15)
# ax1.set_ylim([0,8])
# ax1.legend(ncol=3,loc='upper left',fontsize = 15)  
# ax1.set_xticklabels('')
# p4=ax11.bar([pos + 0.4 for pos in size],np.array(arc_exp['EXP-A'].nobs)/1000,color='k',hatch='x',edgecolor='black',\
#           width=0.2, align='center', label='$NOBS$')
# ax11.legend(loc='upper right',fontsize = 15)
# ax11.tick_params(labelsize=15)
# ax11.set_ylabel('$N.\ of\ SLA\ obs\ *1000$',fontsize = 15)
# ax11.grid(linestyle='--')
# ax1.set_position([0.07,0.69,0.86,0.28])
# # PANEL 2
# ax2 = fig.add_subplot(gs[1, :])
# ax2.bar([pos  for pos in size],arc_exp['EXP-S'].corr,color='b',edgecolor='black',\
#         width=0.2, align='center', label='$EXP-S$')
# ax2.bar([pos + 0.2 for pos in size],arc_exp['EXP-A'].corr,color='g',edgecolor='black',\
#         width=0.2, align='center', label='$EXP{-}A$')
# ax2.set_xticks([pos - 0.4 for pos in size], [str(int(dpt)) for dpt in dpt_list[0:-1]])
# ax2.grid()
# ax2.set_xticklabels('')
# ax2.tick_params(labelsize=15)
# ax2.set_ylabel('$CC\ [\#]$',fontsize = 15)
# ax2.set_position([0.07,0.38,0.86,0.28])
# # PANEL 3
# ax3 = fig.add_subplot(gs[2, :])
# ax3.bar([pos  for pos in size],(arc_exp['EXP-S'].stdm-arc_exp['EXP-S'].stdo)*100,color='b',edgecolor='black',\
#         width=0.2, align='center', label='SIM')
# ax3.bar([pos + 0.2 for pos in size],(arc_exp['EXP-A'].stdm-arc_exp['EXP-A'].stdo)*100,color='g',edgecolor='black',\
#         width=0.2, align='center', label='$EXP{-}A$')
# ax3.set_xticks([pos - 0.4 for pos in size], ['$'+str(int(dpt))+'$' for dpt in dpt_list[0:-1]])
# ax3.grid()
# ax3.tick_params(labelsize=15)
# ax3.set_ylabel('$SDE\ [cm]$',fontsize = 15)
# ax3.set_xlabel('$depth\ [m]$',fontsize = 15)
# ax3.set_position([0.07,0.07,0.86,0.28])
# plt.savefig('./figure05.png')











   
