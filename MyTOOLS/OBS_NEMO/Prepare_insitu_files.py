#!/usr/bin/env python

import warnings  # Manage python warnings
import numpy as np  # Maths module
import xarray as xr  # Xarray (netCDF++)
# import datetime as dtime          # Date module
import os, sys, glob, shutil  # System modules
import gsw  # Sea water formulas

warnings.filterwarnings("ignore", message='The Panel class is removed from pandas.')
warnings.filterwarnings("ignore", message='the \'box\' keyword is deprecated and will be removed in a future version.')
warnings.filterwarnings("ignore", message='Mean of empty slice')

print('----------------------------------------------------------------')
print('                        Preprocessing of in-situ                ')
print('----------------------------------------------------------------')

# dataset='CMEMS_GLO_his'
INPUT_DIR = sys.argv[1]
WORK_DIR = sys.argv[2]
OUT_DIR = sys.argv[3]
REGION = sys.argv[4]

# INPUT_DIR='/work/oda/re15918/data/Reanalysis/Inputs/OBS/RAW/'+dataset+'/'
# OUT_DIR='/work/opa/am09320/Reanalysis/Inputs/OBS/Profiles_raw/' + dataset + '.prof' + '/'
os.makedirs(WORK_DIR, exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)

# %% Parameters
# geo = [-6, 36, 30, 46]
if str(REGION) == 'med':
  geo = [-18.125, 36.25, 30.1878, 45.9375]
elif str(REGION) == 'blk':
  geo = [27.37, 41.96, 40.86, 46.80]
keep_name = True
no_TS = True

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
                 'BA': 'XBT profiles',}
obs_type_used = ['PF', 'BA', 'GL', 'XB']


# %% Functions
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


def get_var_adj(var, var_adj, dim):
    # Replace var with var_adj:
    # Only replace profiles when there is at least one data point in profile
    is_adj = (~np.isnan(var_adj)).any(dim=dim)
    return var_adj.where(is_adj, other=var.where(~is_adj))


def extract_obs_type(filename):
    if "LATEST" in filename:
        return filename.split('_')[3]
    else:
        return filename.split('_')[2]


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
i_debug = 0
files_insitu = glob.glob(INPUT_DIR + '/*.nc')
for file in files_insitu:

    # Check if file is already done
    print(file)
    if file in files_done:
        print('File already preprocessed')
        continue

        # Characteristics from filename
    filebase = os.path.basename(file)
    is_ts = '_TS_' in filebase
    obs_type = extract_obs_type(filebase)

    # Do not process unwanted files
    if not (obs_type in obs_type_used):
        print('Platform type not wanted! Not preprocessed')
        fid.write(file + ' : unwanted platform!\n')
        fid.flush()
        continue
    if no_TS and is_ts:
        print('This is a TS file! Not preprocessed')
        fid.write(file + ' : timeseries file!\n')
        fid.flush()
        continue

    # Open yearly file
    nc_insitu = xr.open_dataset(file)

    if str(REGION) == 'med':    # indentation below to be adjusted
        # Rename for SicilyChannel
        if 'POTENTIAL_TEMP' in nc_insitu.variables:
            if 'TEMP' in nc_insitu.variables:
                nc_insitu = nc_insitu.drop(['TEMP', 'TEMP_QC']).rename(
                    {'POTENTIAL_TEMP': 'TEMP', 'POTENTIAL_TEMP_QC': 'TEMP_QC'})
            else:
                nc_insitu = nc_insitu.rename({'POTENTIAL_TEMP': 'TEMP', 'POTENTIAL_TEMP_QC': 'TEMP_QC'})

        # Replace var by var_adjusted
        for var in ['TEMP', 'PSAL', 'PRES']:
            if var + '_ADJUSTED' in nc_insitu.variables:
                # Only profiles where var_adjusted is present
                nc_insitu[var] = get_var_adj(nc_insitu[var], nc_insitu[var + '_ADJUSTED'], 'DEPTH')
    
    elif str(REGION) == 'blk':  # indentation below to be adjusted
            if ('POTENTIAL_TEMP' in nc_insitu.variables):
                nc_insitu = nc_insitu.rename({'POTENTIAL_TEMP':'TEMP','POTENTIAL_TEMP_QC':'TEMP_QC'})

    # Clean dataset
    myvars = []
    for var in ['TIME', 'LATITUDE', 'LONGITUDE', 'TIME_QC', 'POSITION_QC', 'DC_REFERENCE', 'DIRECTION',
                'DEPH', 'DEPH_QC', 'PRES', 'PRES_QC', 'TEMP', 'TEMP_QC', 'PSAL', 'PSAL_QC']:
        if var in nc_insitu.variables:
            myvars.append(var)
    nc_insitu = nc_insitu[myvars]

    # Characteristics of the file
    is_temp = ('TEMP' in nc_insitu.variables)
    is_salt = ('PSAL' in nc_insitu.variables)
    is_pres = ('PRES' in nc_insitu.variables)
    is_dept = ('DEPH' in nc_insitu.variables)

    if not (obs_type in obs_type_list.keys()):
        print('Observation type not supported: ' + obs_type + '! Applying standard preprocessing.')
    else:
        print('Processing: ' + obs_type_list[obs_type])

    # Only process file if TEMP is present
    if not is_temp and not is_salt:
        print('File does not contain TEMP/PSAL! Not preprocessed')
        fid.write(file + ' : no TEMP/PSAL!\n')
        fid.flush()
        continue

    vars_obs = []
    vars_insitu = []
    if is_temp:
        vars_obs = vars_obs + ['TEMP']
        vars_insitu = vars_insitu + ['TEMP', 'TEMP_QC']
    if is_salt:
        vars_obs = vars_obs + ['PSAL']
        vars_insitu = vars_insitu + ['PSAL', 'PSAL_QC']
    if is_dept:
        vars_insitu = vars_insitu + ['DEPH', 'DEPH_QC']
    if is_pres:
        vars_insitu = vars_insitu + ['PRES', 'PRES_QC']

    ###############################################################################
    # Preprocessing

    # Put lon, lat as variables
    nc_insitu['TIME_TMP'] = xr.DataArray(nc_insitu.TIME.values, dims=['PROF'])
    time_attrs = nc_insitu.TIME.attrs
    nc_insitu = nc_insitu.drop('TIME').rename({'PROF': 'TIME'})
    lat_attrs = nc_insitu.LATITUDE.attrs
    lon_attrs = nc_insitu.LONGITUDE.attrs
    lon_val = nc_insitu.LONGITUDE.values
    lat_val = nc_insitu.LATITUDE.values
    print(nc_insitu.dims)
    if nc_insitu.dims['TIME'] == 1:
        lat_val = np.ones(nc_insitu.dims['TIME']) * lat_val
        lon_val = np.ones(nc_insitu.dims['TIME']) * lon_val
        pos_qc_val = np.ones(nc_insitu.dims['TIME']) * nc_insitu.POSITION_QC.values
        nc_insitu.drop(['POSITION_QC'])
        nc_insitu['POSITION_QC'] = xr.DataArray(pos_qc_val, dims=['TIME'])
    else:
        nc_insitu['POSITION_QC'] = nc_insitu.POSITION_QC.rename({'POSITION': 'TIME'})
    nc_insitu['LAT'] = xr.DataArray(lat_val, dims=['TIME'])
    nc_insitu['LON'] = xr.DataArray(lon_val, dims=['TIME'])
    nc_insitu = nc_insitu.drop(['LATITUDE', 'LONGITUDE'])
    
    nc_insitu.load() # to avoid error in nc_insitu.where
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

    # Write family code if not present
    if not ('family_code' in nc_insitu.attrs):
        nc_insitu.attrs['family_code'] = obs_type

    ###############################################################################
    # Separate into 1 file per profile
    date_start = max(np.datetime64('1985-01-01'), nc_insitu.TIME_TMP.values[0].astype('datetime64[D]'))
    date_end = nc_insitu.TIME_TMP.values[-1].astype('datetime64[D]')
    N_days = ((date_end - date_start) / np.timedelta64(1, 'D')).astype(int) + 1
    date_vec = [date_start + np.timedelta64(i, 'D') for i in range(N_days)]

    if N_days == 0:
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
            # Save in file
            
            prefix_file = OUT_DIR + "/" + filebase[:-3]
            # For timeseries we keep only one file per day otherwise it is per profile
            if is_ts:
                file_day = prefix_file + dday.strftime('_%Y%m%d.nc')
                nc_insitu_day = prepare_for_save(nc_insitu_day, vars_obs, vars_insitu)
                nc_insitu_day.to_netcdf(file_day)
                print(' - ' + file_day)
            else:
                N_prof = nc_insitu_day.dims['TIME']
                for i_prof in range(N_prof):
                    file_prof = prefix_file + '_p%02d' % (i_prof + 1) + dday.strftime('_%Y%m%d.nc')
                    nc_insitu_prf = nc_insitu_day.isel(TIME=[i_prof], POSITION=[i_prof])
                    nc_insitu_prf = prepare_for_save(nc_insitu_prf, vars_obs, vars_insitu)
                    nc_insitu_prf.to_netcdf(file_prof)
                    print(' - ' + file_prof)
    fid.write(file + ' : preprocessing OK.\n')
    fid.flush()
fid.close()
