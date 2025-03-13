#!/usr/bin/env python

import warnings  # Manage python warnings
import numpy as np  # Maths module
import xarray as xr  # Xarray (netCDF++)
# import datetime as dtime          # Date module
import os, sys, glob, shutil  # System modules
import gsw  # Sea water formulas
import csv
import subprocess

warnings.filterwarnings("ignore", message='The Panel class is removed from pandas.')
warnings.filterwarnings("ignore", message='the \'box\' keyword is deprecated and will be removed in a future version.')
warnings.filterwarnings("ignore", message='Mean of empty slice')

print('----------------------------------------------------------------')
print('                        Preprocessing of in-situ                ')
print('----------------------------------------------------------------')

# %% Parameters
# dataset = "Merge_v4"
INPUT_DIR = sys.argv[1]
WORK_DIR = sys.argv[2]
OUT_DIR = sys.argv[3]
#dataset = sys.argv[4]
REGION = sys.argv[4]

# %% Directories
if str(REGION) == 'med':
    datasets_all = [
       'CMEMS_MED_his_old',
       'CMEMS_MED_his_new',
       'CMEMS_GLO_his',
       'SDN_v2',
       'Merge_v2',
       'Merge_v3',
       'Merge_v4',
    ]
    geo = [-18.125, 36.25, 30.1878, 45.9375] # med
elif str(REGION) == 'blk':
    datasets_all = [
       'insitu_BS_NRT_CMEMS',
       'Merge_v1'
    ]
    geo = [27.37,41.96,40.86,46.80]   # blk

# INPUT_DIR = '/work/opa/am09320/Reanalysis/Inputs/OBS/MERGED/'+dataset+'.prof/'
# WORK_DIR = '/work/opa/am09320/Reanalysis/Inputs/OBS/PREPROC/'+dataset+'.prep/'
print('Preprocessing dataset: \n' + INPUT_DIR)
# geo = [-6, 36, 30, 46]
#geo = [-6, 36.25, 30.1878, 45.9375]
keep_name = False
no_TS = True

os.makedirs(WORK_DIR, exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)

# N_days = (dtime.datetime(year+1,1,1)-dtime.datetime(year,1,1)).days
# date_vec = [dtime.datetime(year,1,1)+dtime.timedelta(iday) for iday in range(N_days)]

# List of possible observations
obs_type_list = {'PF': 'profiling floats',
                 'GL': 'gliders',
                 'MO': 'fixed buoys, mooring time series, fixed observations',
                 'RF': 'River flow',
                 'HF': 'HF radar',
                 'DB': 'drifting buoys',
                 'DC': 'drifting buoy reporting calculated sea water current',
                 'CT': 'oceanographic CTD profiles',
                 'BO': 'bottle data',
                 'FB': 'ferrybox',
                 'SF': 'towed CTD data (ex: scanfish)',
                 'TG': 'Tide gauge station',
                 'TS': 'thermosalinograph data',
                 'XB': 'XBT or XCTD profiles',
                 'ML': 'mini logger for fishery observing system',
                 'SM': 'Sea mammals data',
                 'XX': 'Not yet identified data type',
                 'BA': 'XBT profiles',
                 }
obs_type_used = ['PF', 'BA','GL','XB']


# obs_type_used = ['PF']

# %% Functions
def is_gap(x, gapmin=40, depthmax=300):
    x_nonan = x[~np.isnan(x)]
    x_nonan = x_nonan[x_nonan < depthmax]
    return np.any(np.diff(x_nonan) > gapmin)


def check_gap(a, dim, gapmin=40, depthmax=300):
    return xr.apply_ufunc(is_gap, a, input_core_dims=[[dim]],
                          kwargs={'gapmin': gapmin, 'depthmax': depthmax},
                          vectorize=True)


def is_inMed(lon, lat):
    # Region
    geo = [-18.125, 36.4, 30.2, 45.9]
    mask_Med = np.logical_and(np.logical_and(lon > geo[0], lon < geo[1]), np.logical_and(lat > geo[2], lat < geo[3]))
    # Remove Bay of Biscay
#    geo = [-7, -.5, 42.9, 46]
#    mask_BoB = np.logical_and(np.logical_and(lon > geo[0], lon < geo[1]), np.logical_and(lat > geo[2], lat < geo[3]))
#    mask_Med = np.logical_and(mask_Med, ~mask_BoB)
    # Remove Black Sea
    geo = [26.9, 44, 40.1, 48]
    mask_BS = np.logical_and(np.logical_and(lon > geo[0], lon < geo[1]), np.logical_and(lat > geo[2], lat < geo[3]))
    mask_Med = np.logical_and(mask_Med, ~mask_BS)
    # Return masks
    return mask_Med

def is_inBS(lon,lat):
    # Region
    mask_BS = np.logical_and(np.logical_and(lon>geo[0],lon<geo[1]),np.logical_and(lat>geo[2],lat<geo[3]))
    return mask_BS

def prepare_for_save(nc_in, vars_obs, vars_insitu):
    # Remove useless dimensions
    nc_out = nc_in.dropna('DEPTH', how='all', subset=vars_obs)
    # No FillValue for coords
    for var in ['LATITUDE', 'LONGITUDE', 'TIME']:
        nc_out[var].attrs['_FillValue'] = False
    # FillValue is -999.0 for variables
    for var in vars_insitu:
        if var in nc_out.variables:
            nc_out[var].attrs['_FillValue'] = -999.0
            nc_out[var] = nc_out[var].fillna(-999.0)
    return nc_out


def extract_obs_type(filename):
    if "LATEST" in filename:
        return filename.split('_')[3]
    else:
        return filename.split('_')[2]

def write_list_rej_data_argo(nc_in,time,lon,lat,pres,temp,psal,wmocode,flag):
    nreja=0
    if (flag ==1) or (flag ==6):
        for ij in range(len(pres[0])):
            nreja=nreja+1
            with open(WORK_DIR + '/' + 'ARGO_LISTRJ.dat', 'a', newline='') as csvfile:
                spamwriter = csv.writer(csvfile, delimiter=',')
                spamwriter.writerow(['%8s'%(str(time[0])[0:10].replace("-", "")),'%10.5f'%lon,'%10.5f'%lat,'%10.5f'%pres[0][ij],
                                    '%10.5f'%temp[0][ij],'%10.5f'%psal[0][ij],
                                    '%7s'%wmocode,'%1d'%flag])
    else:
        for ij in range(len(nc_in['DEPH'].values[0])):
           if (np.isnan(nc_in['DEPH'].values[0][ij]) or
               np.isnan(nc_in['TEMP'].values[0][ij]) or
               np.isnan(nc_in['PSAL'].values[0][ij])):
               nreja=nreja+1
               with open(WORK_DIR + '/' + 'ARGO_LISTRJ.dat', 'a', newline='') as csvfile:
                   spamwriter = csv.writer(csvfile, delimiter=',')
                   spamwriter.writerow(['%8s'%(str(time[0])[0:10].replace("-", "")),'%10.5f'%lon,'%10.5f'%lat,'%10.5f'%pres[0][ij],
                                       '%10.5f'%temp[0][ij],'%10.5f'%psal[0][ij],
                                       '%7s'%wmocode,'%1d'%flag])
    return(nreja)

def write_list_rej_data_xbt(nc_in,time,lon,lat,pres,temp,wmocode,flag):
    nrejx=0
    if (flag ==1):
        for ij in range(len(pres[0])):
            nrejx=nrejx+1
            with open(WORK_DIR + '/' + 'XBT_LISTRJ.dat', 'a', newline='') as csvfile:
                    spamwriter = csv.writer(csvfile, delimiter=',')
                    spamwriter.writerow(['%8s'%(str(time[0])[0:10].replace("-", "")),'%10.5f'%lon,'%10.5f'%lat,'%10.5f'%pres[0][ij],
                                        '%10.5f'%temp[0][ij],'%8s'%wmocode,'%1d'%flag])
    else:
        for ij in range(len(nc_in['DEPH'].values[0])):
            if (np.isnan(nc_in['DEPH'].values[0][ij]) or
                np.isnan(nc_in['TEMP'].values[0][ij])):
                nrejx=nrejx+1
                with open(WORK_DIR + '/' + 'XBT_LISTRJ.dat', 'a', newline='') as csvfile:
                    spamwriter = csv.writer(csvfile, delimiter=',')
                    spamwriter.writerow(['%8s'%(str(time[0])[0:10].replace("-", "")),'%10.5f'%lon,'%10.5f'%lat,'%10.5f'%pres[0][ij],
                                        '%10.5f'%temp[0][ij],'%8s'%wmocode,'%1d'%flag])
    return(nrejx)

def replace_word(A,B):
    cmd="sed -it 's/"+A+"/"+B+"/g' "+WORK_DIR+ "/" + "pprejobs.dat"
    print(cmd)
    repword=subprocess.getstatusoutput(cmd)
    return repword

# Output text file
file_output = WORK_DIR + '/' + 'preprocessed_files.log'
files_done = []
if os.path.exists(file_output):
    fid = open(file_output, 'r+')
    # Get the observation files that were already done
    lines = fid.readlines()
    for line in lines:
        line_split = line.split(' : ')
        files_done.append(line_split[0])
else:
    fid = open(file_output, 'w')

# %% Preprocessing files
# Number of data rejected by QC of D/T/S for ARGO 
noqcaTOT=0 
# Number of data rejected by QC of D/T for XBT
noqcxTOT=0
# Number of data rejected by check on T/S in range for ARGO
nolmaTOT=0
# Number of data rejected by QC for time or position for ARGO
nodpaTOT=0
# Number of data rejected by QC for time or position for XBT
nodpxTOT=0
# Number of data rejected by check on thermocline  ARGO
nohlaTOT=0
# Number of data with depth < 2 for ARGO
ndueaTOT = 0
# Number of data with depth < 2 for XBT
nduexTOT = 0
# Number of entire profiles rejected by QC for ARGO
nprofa=0
# Number of entire profiles rejected by QC for XBT
nprofx=0
# Number of good data for ARGO
nvalsaTOT=0
# Number of good data for XBT
nvalsxTOT=0

files_insitu  = sorted(glob.glob(INPUT_DIR + '/*.nc'))
for file in files_insitu:

    # Check if file is already done
    print(file)
    if file in files_done:
        print('File already preprocessed')
        continue

        # Characteristics from filename
    filebase = os.path.basename(file)
    obs_type = extract_obs_type(filebase)
    is_seadatanet = (filebase.split('_')[0] == 'SD')

    if not (obs_type in obs_type_used):
        print('Platform type not wanted! Not preprocessed')
        fid.write(file + ' : unwanted platform!\n')
        fid.flush()
        continue

    # For ARGO, only take the first profile each day
    if (obs_type == 'PF'):
        pprof = filebase.split('_')[-2]
        if (pprof[0] == 'p') and (pprof[1:].isdigit()):
            if int(pprof[1:]) > 1:
                print('ARGO float: not the first profile! Not preprocessed')
                fid.write(file + ' : not first ARGO profile!\n')
                fid.flush()
                continue

    # Open yearly file
    nc_insitu = xr.open_dataset(file)

    # Clean dataset
    myvars = []
    for var in ['TIME', 'LATITUDE', 'LONGITUDE', 'TIME_QC', 'POSITION_QC', 'DC_REFERENCE', 'DIRECTION', \
                'DEPH', 'DEPH_QC', 'PRES', 'PRES_QC', 'TEMP', 'TEMP_QC', 'PSAL', 'PSAL_QC']:
        if var in nc_insitu.variables:
            myvars.append(var)
    nc_insitu = nc_insitu[myvars]

    if str(REGION) == 'blk':
        if ('POTENTIAL_TEMP' in nc_insitu.variables):
            nc_insitu = nc_insitu.rename({'POTENTIAL_TEMP':'TEMP','POTENTIAL_TEMP_QC':'TEMP_QC'})

    # Characteristics of the file
    is_temp = ('TEMP' in nc_insitu.variables)
    is_salt = ('PSAL' in nc_insitu.variables)
    is_pres = ('PRES' in nc_insitu.variables)
    is_dept = ('DEPH' in nc_insitu.variables)
    is_ts = (nc_insitu.dims['DEPTH'] == 1)  # file is a timeseries
    is_nodepth = (nc_insitu.dims['DEPTH'] == 0)  # when TEMP or PSAL contain only nan

    if not (obs_type in obs_type_list.keys()):
        print('Observation type not supported: ' + obs_type + '! Applying standard preprocessing.')
    else:
        print('Processing: ' + obs_type_list[obs_type])

    # Only process file if TEMP is present
    if not ('TEMP' in nc_insitu.variables):
        print('File does not contain TEMP! Not preprocessed')
        fid.write(file + ' : no TEMP!\n')
        fid.flush()
        continue

    # Check if it is timeseries
    if no_TS and is_ts:
        print('File has only one depth! No timeseries is on...')
        fid.write(file + ' : timeseries not accepted!\n')
        fid.flush()
        continue

    if is_nodepth:
        print('Corrupted file with no depth! Not preprocessed')
        fid.write(file + ' : not accepted!\n')
        fid.flush()
        continue


    vars_obs = []
    vars_insitu = ['DEPH', 'DEPH_QC']
    if is_temp:
        vars_obs = vars_obs + ['TEMP']
        vars_insitu = vars_insitu + ['TEMP', 'TEMP_QC']
    if is_salt:
        vars_obs = vars_obs + ['PSAL']
        vars_insitu = vars_insitu + ['PSAL', 'PSAL_QC']

    ###############################################################################
    # Preprocessing

    # Put lon, lat as variables
    nc_insitu['TIME_TMP'] = xr.DataArray(nc_insitu.TIME.values, dims=['PROF'])
    time_val = nc_insitu.TIME.values
    time_attrs = nc_insitu.TIME.attrs
    nc_insitu = nc_insitu.drop('TIME').rename({'PROF': 'TIME'})
    lat_attrs = nc_insitu.LATITUDE.attrs
    lon_attrs = nc_insitu.LONGITUDE.attrs
    lon_val = nc_insitu.LONGITUDE.values
    lat_val = nc_insitu.LATITUDE.values
    if is_pres:
        pres_val = nc_insitu.PRES.values
    else:
        pres_val = nc_insitu.DEPH.values
    temp_val = nc_insitu.TEMP.values
    if is_salt:
        psal_val = nc_insitu.PSAL.values
    if (nc_insitu.LATITUDE.shape[0] == 1):
        lat_val = np.ones(nc_insitu.dims['TIME']) * lat_val
        lon_val = np.ones(nc_insitu.dims['TIME']) * lon_val
        pos_qc_val = np.ones(nc_insitu.dims['TIME']) * nc_insitu.POSITION_QC.values
        nc_insitu = nc_insitu.drop(['POSITION_QC'])
        nc_insitu['POSITION_QC'] = xr.DataArray(pos_qc_val, dims=['TIME'])
    else:
        nc_insitu['POSITION_QC'] = nc_insitu.POSITION_QC.rename({'POSITION': 'TIME'})
    nc_insitu['LAT'] = xr.DataArray(lat_val, dims=['TIME'])
    nc_insitu['LON'] = xr.DataArray(lon_val, dims=['TIME'])
    nc_insitu = nc_insitu.drop(['LATITUDE', 'LONGITUDE'])
    # Out of domain
    if str(REGION) == 'med':
        nc_insitu = nc_insitu.where(is_inMed(nc_insitu.LON, nc_insitu.LAT), drop=True)
    elif str(REGION) == 'blk':
        nc_insitu = nc_insitu.where(is_inBS(nc_insitu.LON, nc_insitu.LAT), drop=True)

    # Exit if no data
    if nc_insitu.dims['TIME'] == 0:
        print('File has no data on region!')
        fid.write(file + ' : no data in region!\n')
        fid.flush()
        continue
# Mettere qui l'append a FILEOBS del nome del file per GLO
    # if TIME.dims > 1 check "QC flag on depth":
    # if we have n dim of TIME and DEPTH, a single dim of DEPTH (with all nan) can drop all TIME dims
    if nc_insitu.dims['TIME'] > 1:
        raise Exception("File with TIME dimension > 1 can't pass quality control")

    # QC for time and position
    # print('TIME/POS QC: '+str(nc_insitu.TIME_QC.values[0])+'/'+str(nc_insitu.POSITION_QC.values[0]))
    if not is_seadatanet:
        nc_insitu = nc_insitu.where((nc_insitu.TIME_QC == 1) & (nc_insitu.POSITION_QC == 1), drop=True)

        # Exit if no data
        if nc_insitu.dims['TIME'] == 0:
            print('File has bad QC for time or position!')
            fid.write(file + ' : bad QC for time/position!\n')
            fid.flush()
#            continue
            fqfpos=1
            if (obs_type == 'PF'):
                nodpa=write_list_rej_data_argo(nc_insitu,time_val,lon_val[0],lat_val[0],pres_val,temp_val,psal_val,nc_insitu.attrs["wmo_platform_code"],fqfpos)
                nodpaTOT=nodpaTOT+(nodpa*2)
                nprofa=nprofa+1
            else:
                nodpx=write_list_rej_data_xbt(nc_insitu,time_val,lon_val[0],lat_val[0],pres_val,temp_val,nc_insitu.attrs["wmo_platform_code"],fqfpos)
                nodpxTOT=nodpxTOT+nodpx
                nprofx=nprofx+1
            continue

    # Direction must be "down"
    if ('DIRECTION' in nc_insitu.variables) and (obs_type == 'CT'):
        nc_insitu = nc_insitu.where(nc_insitu.DIRECTION == b'D', drop=True).drop(['DIRECTION'])

    # Exit if no data
    if nc_insitu.dims['TIME'] == 0:
        print('File has only ascending profiles!')
        fid.write(file + ' : Only ascending profiles!\n')
        fid.flush()
        continue

    # Check if XBT are with depths every 2 m (good) or every 0.7 m (raw)
    # If raw exit
    if (obs_type == 'BA') or (obs_type == 'XB'):
        if (pres_val[0][1]-pres_val[0][0]) < 2:
            print('File has only raw profiles!')
            fid.write(file + ' : Only raw profiles!\n')
            fid.flush()
            continue 

    # Pressure VS depth
    if is_pres:
        dept_pres = -gsw.geostrophy.z_from_p(nc_insitu.PRES.transpose(), nc_insitu.LAT).transpose()
        if is_dept:
            # Merge depth and depth_from_pres arrays
            dept_pres = np.nanmean(np.stack((nc_insitu.DEPH.values, dept_pres)), axis=0)
            nc_insitu['DEPH_QC'] = nc_insitu['DEPH_QC'].where(nc_insitu.PRES_QC != 1, other=1)
        else:
            nc_insitu['DEPH_QC'] = nc_insitu.PRES_QC
        nc_insitu['DEPH'] = xr.DataArray(dept_pres, dims=['TIME', 'DEPTH'])
        nc_insitu = nc_insitu.drop(['PRES', 'PRES_QC'])

    # Exit if no data
    if nc_insitu.dims['TIME'] == 0:
        print('Problem in Pressure/depth!')
        fid.write(file + ' : pb pressure/depth!\n')
        fid.flush()
        continue

    # Remove useless dimensions
    nc_insitu = nc_insitu.dropna('TIME', how='all', subset=vars_obs)

    # Separate 2D variables
    nc_TS = nc_insitu[vars_insitu]
    nc_not_TS = nc_insitu.drop(vars_insitu)

    # Exit if no data
    if nc_TS.dropna('TIME', how='all', subset=vars_obs).dims['TIME'] == 0:
        print('No data in file')
        fid.write(file + ' : no data in file!\n')
        fid.flush()
        if (obs_type == 'PF'):
            nprofa=nprofa+1
        else:
            nprofx=nprofx+1 
        continue

    # QC for depth: this check puts Nan in the levels where DEPH_QC != 1
    nc_TS: xr.Dataset = nc_TS.where((nc_TS.DEPH_QC == 1))
#    print(nc_TS)
    # Exit if no data 
    if nc_TS.dropna('TIME', how='all', subset=vars_obs).dims['TIME'] == 0:
        print('QC of DEPTH failed')
        fid.write(file + ' : depth QC fail!\n')
        fid.flush()
        fqfpts=4
        if (obs_type == 'PF'):
            noqca=write_list_rej_data_argo(nc_TS,time_val,lon_val[0],lat_val[0],pres_val,temp_val,psal_val,nc_insitu.attrs["wmo_platform_code"],fqfpts)
            noqcaTOT=noqcaTOT+(noqca*2)
            nprofa=nprofa+1
        else:
            noqcx=write_list_rej_data_xbt(nc_TS,time_val,lon_val[0],lat_val[0],pres_val,temp_val,nc_insitu.attrs["wmo_platform_code"],fqfpts)
            noqcxTOT=noqcxTOT+noqcx
            nprofx=nprofx+1
        continue

    # QC for temperature: this check delete levels where TEMP_QC != 1
    nc_TS['TEMP'] = nc_TS.TEMP.where((nc_TS.TEMP_QC == 1))
    # QC for salinity: this check delete levels where PSAL_QC != 1
    if is_salt:
        nc_TS['PSAL'] = nc_TS.PSAL.where((nc_TS.PSAL_QC == 1))
    fqfpts=4
    if (obs_type == 'PF'):
        noqca=write_list_rej_data_argo(nc_TS,time_val,lon_val[0],lat_val[0],pres_val,temp_val,psal_val,nc_insitu.attrs["wmo_platform_code"],fqfpts)
        noqcaTOT=noqcaTOT+(noqca*2)
    else:
        noqcx=write_list_rej_data_xbt(nc_TS,time_val,lon_val[0],lat_val[0],pres_val,temp_val,nc_insitu.attrs["wmo_platform_code"],fqfpts)
        noqcxTOT=noqcxTOT+noqcx
    # Exit if no data - QC flag on depth for T/S
    if nc_TS.dropna('DEPTH', how='any', subset=vars_obs).dims['DEPTH'] == 0:
        print('QC of T/S failed')
        fid.write(file + ' : T/S QC fail!\n')
        fid.flush()
        if (obs_type == 'PF'):
            nprofa=nprofa+1
        else:
            nprofx=nprofx+1
        continue

    # Correct range for temperature
    nc_TS['TEMP'] = nc_TS['TEMP'].where((nc_TS['TEMP'] < 35) & (nc_TS['TEMP'] > 0))
    # Correct range for salinity
    if is_salt:
        nc_TS['PSAL'] = nc_TS['PSAL'].where((nc_TS['PSAL'] < 45) & (nc_TS['PSAL'] > 0))

    # Exit if no data
    if nc_TS.dropna('TIME', how='all', subset=vars_obs).dims['TIME'] == 0:
        print('T/S not in range')
        fid.write(file + ' : T/S not within range!\n')
        fid.flush()
        fqflm=5
        if (obs_type == 'PF'):
            nolma=write_list_rej_data_argo(nc_TS,time_val,lon_val[0],lat_val[0],pres_val,temp_val,psal_val,nc_insitu.attrs["wmo_platform_code"],fqflm)
            nolmaTOT=nolmaTOT+(nolma*2)
            nprofa=nprofa+1
        continue
    # Check gap for thermocline (Only ARGO and gliders)
    if obs_type in ['PF', 'GL']:
        depth_nan = nc_TS.DEPH.where(~np.isnan(nc_TS.TEMP))
        nc_TS = nc_TS.where(~check_gap(depth_nan, dim='DEPTH', gapmin=40, depthmax=300))

    if ((nc_TS.dropna('TIME', how='all', subset=vars_obs).dims['TIME'] == 0) or
       (nc_TS.dropna('DEPTH', how='any', subset=vars_obs).DEPH.values[0][0] > 40 )):  
        print('Thermocline check failed')
        fid.write(file + ' : Thermocline check failed!\n')
        fid.flush()
        fqfhl=6
        if (obs_type == 'PF'):
            nohla=write_list_rej_data_argo(nc_insitu,time_val,lon_val[0],lat_val[0],pres_val,temp_val,psal_val,nc_insitu.attrs["wmo_platform_code"],fqfhl)
            nohlaTOT=nohlaTOT+(nohla*2)
            nprofa=nprofa+1
        continue

    # Minimum 2 meters depth (Only ARGO and XBT)
    if obs_type in ['PF', 'BA', 'XB','GL']:
        A=nc_TS.where((nc_TS.DEPH < 2))
        nc_TS = nc_TS.where((nc_TS.DEPH >= 2))
        if (obs_type == 'PF'):
            nvalsaTOT=nvalsaTOT+((nc_TS.dropna('DEPTH', how='any', subset=vars_obs).dims['DEPTH'])*2)
            ndueaTOT=ndueaTOT+((A.dropna('DEPTH', how='any', subset=vars_obs).dims['DEPTH'])*2)
        else:
            nvalsxTOT=nvalsxTOT+(nc_TS.dropna('DEPTH', how='any', subset=vars_obs).dims['DEPTH'])
            nduexTOT=nduexTOT+(A.dropna('DEPTH', how='any', subset=vars_obs).dims['DEPTH'])
    # Exit if no data
    if nc_TS.dropna('TIME', how='all', subset=vars_obs).dims['TIME'] == 0:
        print('No data below 2m')
        fid.write(file + ' : too shallow data!\n')
        fid.flush()
        if (obs_type == 'PF'):
            nprofa=nprofa+1
        else:
            nprofx=nprofx+1
        continue
    # Remerge the two datasets
    nc_insitu = xr.merge((nc_not_TS, nc_TS))
    nc_insitu.attrs = nc_not_TS.attrs

    # Remove useless dimensions
    nc_insitu = nc_insitu.dropna('TIME', how='all', subset=vars_obs)
    nc_insitu = nc_insitu.dropna('DEPTH', how='any', subset=vars_obs)

    # Change type of DC_REFERENCE
    if 'DC_REFERENCE' in nc_insitu.variables:
        attrs_tmp = nc_insitu.DC_REFERENCE.attrs
        nc_insitu['DC_REFERENCE'] = (nc_insitu.DC_REFERENCE).astype(str)

    # Write family code if not present
    if not ('family_code' in nc_insitu.attrs):
        nc_insitu.attrs['family_code'] = obs_type

    # Exit if no data
    if (nc_insitu.dims['TIME'] == 0):
        print('File has failed preprocessing!')
        fid.write(file + ' : failed preprocessing!\n')
        fid.flush()
        continue

    ###############################################################################
    # Separate in daily files
    date_start = max(np.datetime64('1985-01-01'), nc_insitu.TIME_TMP.values[0].astype('datetime64[D]'))
    date_end = nc_insitu.TIME_TMP.values[-1].astype('datetime64[D]')
    N_days = ((date_end - date_start) / np.timedelta64(1, 'D')).astype(int) + 1
    date_vec = [date_start + np.timedelta64(i, 'D') for i in range(N_days)]

    if (N_days == 0):
        print('No data in the period!')
        fid.write(file + ' : no data in period!\n')
        fid.flush()
        continue

    for dday_np64 in date_vec:
        dday = dday_np64.item()
        idx_day = (nc_insitu.TIME_TMP >= dday_np64) & (nc_insitu.TIME_TMP < (dday_np64 + np.timedelta64(1, 'D')))
        if np.sum(idx_day) > 0:
            # Create daily xarray
            nc_insitu_day = nc_insitu.where((idx_day), drop=True)
            # Put back lon and lat as coords
            nc_insitu_day['POSITION_QC'] = nc_insitu_day.POSITION_QC.rename({'TIME': 'POSITION'})
            nc_insitu_day['LATITUDE'] = xr.DataArray(nc_insitu_day.LAT.values, dims=['POSITION'], attrs=lat_attrs)
            nc_insitu_day['LONGITUDE'] = xr.DataArray(nc_insitu_day.LON.values, dims=['POSITION'], attrs=lon_attrs)
            # Add the time vector
            time_out = [(ddate - np.datetime64('1950-01-01')) / np.timedelta64(1, 'D') for ddate in
                        nc_insitu_day.TIME_TMP.values]
            nc_insitu_day['TIME'] = xr.DataArray(time_out, dims=['TIME'], attrs=time_attrs)
            nc_insitu_day['TIME'].attrs['units'] = 'days since 1950-1-1 00:00:00'
            nc_insitu_day = nc_insitu_day.drop(['LON', 'LAT', 'TIME_TMP'])
            # Preparation for writing netCDF
            nc_insitu_day = prepare_for_save(nc_insitu_day, vars_obs, vars_insitu)
            
            # Save in file
            if keep_name:
                file_day = f'{OUT_DIR}/{filebase}'.replace('GL_', 'MO_')
            else:
                filebase_no_latest = filebase.replace("_LATEST", "")
                file_split = filebase_no_latest.split('_')
                file_date = dday.strftime('%Y%m%d')
                file_float = file_split[3]
                file_obs_type = obs_type
                if not file_obs_type == 'PF':
                    file_obs_type = 'GL' if is_salt else 'BA'
                
                file_day = f'{OUT_DIR}/MO_PR_{file_obs_type}_{file_float}_{file_date}.nc'
            nc_insitu_day.to_netcdf(file_day)
            print(' - ' + file_day)
    fid.write(file + ' : preprocessing OK.\n')
    fid.flush()
# Number of input data for ARGO
naTOT=nvalsaTOT+nolmaTOT+noqcaTOT+nodpaTOT+nohlaTOT+ndueaTOT
# Number of input data for XBT
nxTOT=nvalsxTOT+noqcxTOT+nodpxTOT+nduexTOT
# Put ARGO statistics in the final tab

replace_word("obsa",str(nvalsaTOT))
replace_word("nodpa",str(nodpaTOT))
replace_word("noqca",str(noqcaTOT))
replace_word("nolma",str(nolmaTOT))
replace_word("nohla",str(nohlaTOT))
replace_word("nduea",str(ndueaTOT))
replace_word("pda",str(nprofa))  
replace_word("ta",str(naTOT)) 

# Put XBT statistics in the final tab
replace_word("obsx",str(nvalsxTOT))
replace_word("nodpx",str(nodpxTOT))
replace_word("noqcx",str(noqcxTOT))
replace_word("nduex",str(nduexTOT))
replace_word("pdx",str(nprofx))
replace_word("tx",str(nxTOT))
fid.close()
