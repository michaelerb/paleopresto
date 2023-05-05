#==============================================================================
# Make a standardized netcdf file for the LMR.
#    author: Michael P. Erb
#    date  : 5/4/2023
#==============================================================================

import sys
import numpy as np
import xarray as xr
import utils


#%% SETTINGS

version_txt = '2.0'
#version_txt = '2.1'
var_txt     = 'tas'
#var_txt     = 'precip'

#version_txt = sys.argv[1]
#var_txt     = sys.argv[2]


#%% LOAD DATA
# Note: data can be downloaded at: https://www.ncei.noaa.gov/access/paleo-search/study/27850

print('=== Processing LMR, v'+version_txt+' ===')
data_dir = '/projects/pd_lab/data/paleoclimate_reconstructions/LMR_final/'

# Load data
if var_txt == 'tas':
    #
    data_xarray = xr.open_dataset(data_dir+'air_MCruns_ensemble_mean_LMRv'+version_txt+'.nc')
    var_ens = data_xarray['air'].values
    time    = data_xarray['time'].values
    lat     = data_xarray['lat'].values
    lon     = data_xarray['lon'].values
    data_xarray.close()
    #
    data_xarray = xr.open_dataset(data_dir+'gmt_MCruns_ensemble_full_LMRv'+version_txt+'.nc')
    var_global = data_xarray['gmt'].values
    data_xarray.close()
    var_global = np.reshape(var_global,(2001,20*100))
    #
elif var_txt == 'precip':
    #
    data_xarray = xr.open_dataset(data_dir+'prate_MCruns_ensemble_mean_LMRv'+version_txt+'.nc')
    var_ens = data_xarray['prate'].values
    time    = data_xarray['time'].values
    lat     = data_xarray['lat'].values
    lon     = data_xarray['lon'].values
    var_ens = var_ens*60*60*24 # Convert precipitation units from kg/m2/s to mm/day
    #
    # Compute the global mean of precip
    lat_weights  = np.cos(np.deg2rad(data_xarray.lat))
    var_global = data_xarray.prate.weighted(lat_weights).mean(('lon','lat'))
    data_xarray.close()

# Get the years and ages
year = np.array([time_selected.year for time_selected in time])
age = 1950-year


#%% CALCULATIONS

# Calculate lat and lon bounds
lat_bounds,lon_bounds = utils.bounding_latlon(lat,lon)

# Process data
var_mean = np.mean(var_ens,axis=1)

# Format the variables
var_global = np.swapaxes(var_global,0,1)
var_ens    = np.swapaxes(var_ens,   0,1)
var_global = np.expand_dims(var_global,axis=0)
var_ens    = np.expand_dims(var_ens,   axis=0)
var_mean   = np.expand_dims(var_mean,  axis=0)

# Get other metadata
methods = ['LMR']
ens_spatial = np.arange(var_ens.shape[1])+1
ens_global  = np.arange(var_global.shape[1])+1

# If this data can't be reformatted to the standard format, add a note here 
notes = ['ens_are_not_available']


#%% FORMAT DATA

# Create new array
if   var_txt == 'tas':    units = 'degrees Celsius'
elif var_txt == 'precip': units = 'mm per day'
data_xarray_output = xr.Dataset(
    {
        var_txt+'_global':(['method','ens_global','age'],             var_global,{'units':units}),
        var_txt+'_mean':  (['method','age','lat','lon'],              var_mean,  {'units':units}),
        var_txt+'_ens':   (['method','ens_spatial','age','lat','lon'],var_ens,   {'units':units})
    },
    coords={
        'method':     (['method'],methods),
        'notes':      (['notes'],notes),
        'ens_global': (['ens_global'],ens_global,{'description':'ensemble members'}),
        'ens_spatial':(['ens_spatial'],ens_spatial,{'description':'means of the 20 Monte Carlo interations'}),
        'age':        (['age'],age,{'units':'yr BP'}),
        'lat':        (['lat'],lat,{'units':'degrees_north'}),
        'lon':        (['lon'],lon,{'units':'degrees_east'}),
        'lat_bounds': (['lat_bounds'],lat_bounds,{'units':'degrees_north'}),
        'lon_bounds': (['lon_bounds'],lon_bounds,{'units':'degrees_east'}),
    },
)


#%% SAVE DATA

# Save new array
output_dir = '/projects/pd_lab/data/paleoclimate_reconstructions/presto_format/'
output_name = 'lmr_v'+version_txt.replace('.','_')+'_0_'+var_txt+'_annual.nc'
data_xarray_output.to_netcdf(output_dir+output_name)
