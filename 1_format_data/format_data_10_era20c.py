#==============================================================================
# Make a standardized netcdf file for the ERA-20C reanalysis.
#    author: Michael P. Erb
#    date  : 9/12/2023
#=============================================================================

import sys
import numpy as np
import xarray as xr
import utils
import calendar


# SETTINGS

#quantity_txt = 'annual'
#quantity_txt = 'jja'
#quantity_txt = 'djf'
quantity_txt = sys.argv[1]


#%% LOAD DATA
# Note: data can be downloaded at: https://apps.ecmwf.int/datasets/data/era20c-moda/levtype=sfc/type=an/

print('=== Processing ERA-20C ===')

# Load data
data_dir = '/projects/pd_lab/data/modern_datasets/ERA20C/'
data_xarray = xr.open_dataset(data_dir+'t2m_monthly_190001_to_201012_era20c.nc')
var_spatial_monthly = data_xarray['2T_GDS4_SFC_S123'].values
lon = data_xarray['g4_lon_2'].values
lat = data_xarray['g4_lat_1'].values
#time_str = handle['initial_time0'].values

years = np.arange(1900,2011,1)
age = 1950-years

# Compute the global mean
lat_weights = np.cos(np.deg2rad(data_xarray.g4_lat_1))
var_global_monthly = data_xarray['2T_GDS4_SFC_S123'].weighted(lat_weights).mean(('g4_lon_2','g4_lat_1')).values
data_xarray.close()


#%% CALCULATIONS

# Reshape the data
nyears = int(var_spatial_monthly.shape[0]/12)
nlat   = var_spatial_monthly.shape[1]
nlon   = var_spatial_monthly.shape[2]
var_spatial_monthly_reshape = np.reshape(var_spatial_monthly,(nyears,12,nlat,nlon)) 
var_global_monthly_reshape  = np.reshape(var_global_monthly, (nyears,12)) 

# Compute the annual-means weighted by the correct number of days in each month.
var_spatial_mean = np.zeros((nyears,nlat,nlon)); var_spatial_mean[:] = np.nan
var_global_mean  = np.zeros((nyears));           var_global_mean[:]  = np.nan
i=0;year=years[i]
for i,year in enumerate(years):
    if calendar.isleap(year): days_in_months = np.array([31,29,31,30,31,30,31,31,30,31,30,31])
    else:                     days_in_months = np.array([31,28,31,30,31,30,31,31,30,31,30,31])
    if quantity_txt == 'annual':
        var_spatial_mean[i,:,:] = np.average(var_spatial_monthly_reshape[i,:,:,:],axis=0,weights=days_in_months)
        var_global_mean[i]      = np.average(var_global_monthly_reshape[i,:],     axis=0,weights=days_in_months)
    elif quantity_txt == 'jja':
        var_spatial_mean[i,:,:] = np.average(var_spatial_monthly_reshape[i,[5,6,7],:,:],axis=0,weights=days_in_months[[5,6,7]])
        var_global_mean[i]      = np.average(var_global_monthly_reshape[i,[5,6,7]],     axis=0,weights=days_in_months[[5,6,7]])
    elif quantity_txt == 'djf':
        #
        if calendar.isleap(year+1): weights_DJF = np.array([31,31,29])
        else:                       weights_DJF = np.array([31,31,28])
        #
        data_D = np.expand_dims(var_spatial_monthly_reshape[i,11,:,:],axis=0)
        try:    data_JF = var_spatial_monthly_reshape[i+1,[0,1],:,:]
        except: data_JF = np.zeros((2,nlat,nlon)); data_JF[:] = np.nan
        data_DJF = np.concatenate((data_D,data_JF),axis=0)
        var_spatial_mean[i,:,:] = np.average(data_DJF,axis=0,weights=weights_DJF)
        #
        data_global_D = np.expand_dims(var_global_monthly_reshape[i,11],axis=0)
        try:    data_global_JF = var_global_monthly_reshape[i+1,[0,1]]
        except: data_global_JF = np.zeros((2)); data_global_JF[:] = np.nan
        data_global_DJF = np.concatenate((data_global_D,data_global_JF),axis=0)
        var_global_mean[i] = np.average(data_global_DJF,axis=0,weights=weights_DJF)


#%% CALCULATIONS 2

# Calculate lat and lon bounds
lat_bounds,lon_bounds = utils.bounding_latlon(lat,lon)

# Format the variables
var_spatial_mean = np.expand_dims(var_spatial_mean,axis=0)
var_global_mean  = np.expand_dims(var_global_mean, axis=0)

# ERA-20C only has one ensemble member, so the mean and the ensemble member are the same. #TOOD: Is this true?
var_spatial_members = np.expand_dims(var_spatial_mean,axis=1)
var_global_members  = np.expand_dims(var_global_mean, axis=1)

# Get other metadata
methods = ['ERA-20C']
ens_spatial = np.array([1])
ens_global = ens_spatial

# If this data can't be reformatted to the standard format, add a note here 
notes = ['ERA-20C only has one ensemble member, so the mean and the ensemble member are the same.']

# Check the shape of the variables
print(var_spatial_members.shape)
print(var_spatial_mean.shape)
print(var_global_members.shape)
print(var_global_mean.shape)


#%% FORMAT DATA

# Create new array
data_xarray_output = xr.Dataset(
    {
        'tas_global_mean':    (['method','age'],                          var_global_mean,    {'units':'degrees Celsius'}),
        'tas_global_members': (['method','ens_global','age'],             var_global_members, {'units':'degrees Celsius'}),
        'tas_spatial_mean':   (['method','age','lat','lon'],              var_spatial_mean,   {'units':'degrees Celsius'}),
        'tas_spatial_members':(['method','ens_spatial','age','lat','lon'],var_spatial_members,{'units':'degrees Celsius'})
    },
    coords={
        'method':     (['method'],methods),
        'notes':      (['notes'],notes),
        'ens_global': (['ens_global'],ens_global,{'description':'ensemble members'}),
        'ens_spatial':(['ens_spatial'],ens_spatial,{'description':'ensemble members'}),
        'age':        (['age'],age,{'units':'yr BP'}),
        'lat':        (['lat'],lat,{'units':'degrees_north'}),
        'lon':        (['lon'],lon,{'units':'degrees_east'}),
        'lat_bounds': (['lat_bounds'],lat_bounds,{'units':'degrees_north'}),
        'lon_bounds': (['lon_bounds'],lon_bounds,{'units':'degrees_east'}),
    },
    attrs={
        'dataset_name':      'ERA-20C',
        'dataset_source_url':'https://apps.ecmwf.int/datasets/data/era20c-moda/levtype=sfc/type=an/',
    },
)


#%% SAVE DATA

# Save new array
output_dir = '/projects/pd_lab/data/paleoclimate_reconstructions/presto_format/'
data_xarray_output.to_netcdf(output_dir+'era20c_v1_0_0_tas_'+quantity_txt+'.nc')

