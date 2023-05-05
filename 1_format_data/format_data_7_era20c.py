#==============================================================================
# Make a standardized netcdf file for the ERA-20C reanalysis.
#    author: Michael P. Erb
#    date  : 5/5/2023
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
tas = data_xarray['2T_GDS4_SFC_S123'].values
lon = data_xarray['g4_lon_2'].values
lat = data_xarray['g4_lat_1'].values
#time_str = handle['initial_time0'].values

years = np.arange(1900,2011,1)
age = 1950-years

# Compute the global mean
lat_weights = np.cos(np.deg2rad(data_xarray.g4_lat_1))
tas_global = data_xarray['2T_GDS4_SFC_S123'].weighted(lat_weights).mean(('g4_lon_2','g4_lat_1')).values
data_xarray.close()


#%% CALCULATIONS

# Reshape the data
nyears = int(tas.shape[0]/12)
nlat   = tas.shape[1]
nlon   = tas.shape[2]
tas_reshape        = np.reshape(tas,(nyears,12,nlat,nlon)) 
tas_global_reshape = np.reshape(tas_global,(nyears,12)) 

# Compute the annual-means weighted by the correct number of days in each month.
tas_mean        = np.zeros((nyears,nlat,nlon)); tas_mean[:]        = np.nan
tas_global_mean = np.zeros((nyears));           tas_global_mean[:] = np.nan
i=0;year=years[i]
for i,year in enumerate(years):
    if calendar.isleap(year): days_in_months = np.array([31,29,31,30,31,30,31,31,30,31,30,31])
    else:                     days_in_months = np.array([31,28,31,30,31,30,31,31,30,31,30,31])
    if quantity_txt == 'annual':
        tas_mean[i,:,:]    = np.average(tas_reshape[i,:,:,:],   axis=0,weights=days_in_months)
        tas_global_mean[i] = np.average(tas_global_reshape[i,:],axis=0,weights=days_in_months)
    elif quantity_txt == 'jja':
        tas_mean[i,:,:]    = np.average(tas_reshape[i,[5,6,7],:,:],   axis=0,weights=days_in_months[[5,6,7]])
        tas_global_mean[i] = np.average(tas_global_reshape[i,[5,6,7]],axis=0,weights=days_in_months[[5,6,7]])
    elif quantity_txt == 'djf':
        #
        if calendar.isleap(year+1): weights_DJF = np.array([31,31,29])
        else:                       weights_DJF = np.array([31,31,28])
        #
        data_D = np.expand_dims(tas_reshape[i,11,:,:],axis=0)
        try:    data_JF = tas_reshape[i+1,[0,1],:,:]
        except: data_JF = np.zeros((2,nlat,nlon)); data_JF[:] = np.nan
        data_DJF = np.concatenate((data_D,data_JF),axis=0)
        tas_mean[i,:,:] = np.average(data_DJF,axis=0,weights=weights_DJF)
        #
        data_global_D = np.expand_dims(tas_global_reshape[i,11],axis=0)
        try:    data_global_JF = tas_global_reshape[i+1,[0,1]]
        except: data_global_JF = np.zeros((2)); data_global_JF[:] = np.nan
        data_global_DJF = np.concatenate((data_global_D,data_global_JF),axis=0)
        tas_global_mean[i] = np.average(data_global_DJF,axis=0,weights=weights_DJF)


#%% CALCULATIONS 2

# Calculate lat and lon bounds
lat_bounds,lon_bounds = utils.bounding_latlon(lat,lon)

# Format the variables
tas_ens    = np.expand_dims(np.expand_dims(tas_mean,axis=0),axis=0)
tas_global = np.expand_dims(np.expand_dims(tas_global_mean,axis=0),axis=0)
tas_mean   = np.expand_dims(tas_mean,axis=0)

# Get other metadata
methods = ['ERA-20C']
ens_spatial = np.array([1])
ens_global = ens_spatial

# If this data can't be reformatted to the standard format, add a note here 
notes = ['']


#%% FORMAT DATA

# Create new array
data_xarray_output = xr.Dataset(
    {
        'tas_global':(['method','ens_global','age'],             tas_global,{'units':'degrees Celsius'}),
        'tas_mean':  (['method','age','lat','lon'],              tas_mean,  {'units':'degrees Celsius'}),
        'tas_ens':   (['method','ens_spatial','age','lat','lon'],tas_ens,   {'units':'degrees Celsius'})
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
)


#%% SAVE DATA

# Save new array
output_dir = '/projects/pd_lab/data/paleoclimate_reconstructions/presto_format/'
output_name = 'era20c_v1_0_0_tas_'+quantity_txt+'.nc'
data_xarray_output.to_netcdf(output_dir+output_name)
