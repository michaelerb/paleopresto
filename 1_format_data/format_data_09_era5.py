#==============================================================================
# Make a standardized netcdf file for the ERA5 reanalysis.
#    author: Michael P. Erb
#    date  : 3/27/2024
#=============================================================================

import sys
import numpy as np
import xarray as xr
import utils
import calendar


#%% SETTINGS

#var_txt = 'tas'
#var_txt = 'precip'
#var_txt = 'slp'
#quantity_txt = 'annual'
#quantity_txt = 'jja'
#quantity_txt = 'djf'
var_txt = sys.argv[1]; quantity_txt = sys.argv[2]


#%% LOAD DATA
# Note: data can be downloaded at: https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels-monthly-means?tab=form

if   var_txt == 'tas':    file_var_txt = 't2m'
elif var_txt == 'precip': file_var_txt = 'tp'
elif var_txt == 'slp':    file_var_txt = 'msl'

print('=== Processing ERA5 ===')

# Load mean data
data_dir = '/projects/pd_lab/data/modern_datasets/ERA5/'
data_xarray = xr.open_dataset(data_dir+'era5_monthly_'+file_var_txt+'.nc')
var_spatial_monthly = data_xarray[file_var_txt].values[:,0,:,:]
lon = data_xarray['longitude'].values
lat = data_xarray['latitude'].values
month_lengths = data_xarray.time.dt.days_in_month.values

# Compute the global mean for ensemble members
lat_weights = np.cos(np.deg2rad(data_xarray.latitude))
var_global_monthly = data_xarray[file_var_txt].weighted(lat_weights).mean(('longitude','latitude')).values[:,0]
data_xarray.close()

years = np.arange(1940,2024,1)
age = 1950-years


#%% CALCULATIONS 1

# Remove months in incomplete years
nmonths = var_spatial_monthly.shape[0]
ind_end = int(np.floor(nmonths/12))*12
var_spatial_monthly = var_spatial_monthly[:ind_end,:,:]
month_lengths       = month_lengths[:ind_end]
var_global_monthly  = var_global_monthly[:ind_end]

# Change units and set unit text
if var_txt == 'tas':
    # Convert temperature from K to degC
    var_spatial_monthly = var_spatial_monthly-273.15
    var_global_monthly  = var_global_monthly-273.15
    unit_txt = 'degrees Celsius'
elif var_txt == 'precip':
    # Convert precipitation from m/day to mm/day
    var_spatial_monthly = var_spatial_monthly*1000
    var_global_monthly  = var_global_monthly*1000
    unit_txt = 'mm/day'
elif var_txt == 'slp':
    # Convert SLP from Pa to hPa
    var_spatial_monthly = var_spatial_monthly/100
    var_global_monthly  = var_global_monthly/100
    unit_txt = 'hPa'

# Get dimensions
ntime = var_spatial_monthly.shape[0]
nlat  = var_spatial_monthly.shape[1]
nlon  = var_spatial_monthly.shape[2]


#%% CALCULATIONS 2

# Reshape the data
nyears = int(ntime/12)
var_spatial_monthly_reshape = np.reshape(var_spatial_monthly,(nyears,12,nlat,nlon))
var_global_monthly_reshape  = np.reshape(var_global_monthly, (nyears,12))
month_lengths_reshape       = np.reshape(month_lengths,      (nyears,12))

# Compute the means weighted by the correct number of days in each month.
var_spatial_mean = np.zeros((nyears,nlat,nlon)); var_spatial_mean[:] = np.nan
var_global_mean  = np.zeros((nyears));           var_global_mean[:]  = np.nan
i=0;year=years[i]
for i,year in enumerate(years):
    days_in_months = month_lengths_reshape[i,:]
    if quantity_txt == 'annual':
        var_spatial_mean[i,:,:] = np.average(var_spatial_monthly_reshape[i,:,:,:],axis=0,weights=days_in_months)
        var_global_mean[i]      = np.average(var_global_monthly_reshape[i,:],     axis=0,weights=days_in_months)
    if quantity_txt == 'jja':
        var_spatial_mean[i,:,:] = np.average(var_spatial_monthly_reshape[i,[5,6,7],:,:],axis=0,weights=days_in_months[[5,6,7]])
        var_global_mean[i]      = np.average(var_global_monthly_reshape[i,[5,6,7]],     axis=0,weights=days_in_months[[5,6,7]])
    if quantity_txt == 'djf':
        # Get DJF weights
        try: weights_DJF = np.concatenate((np.expand_dims(days_in_months[11],axis=0),month_lengths_reshape[i+1,:2]),axis=0)
        except:
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

"""
#%% CALCULATIONS 3

# Find the years with data in all 12 months
var_spatial_mean_mean = np.nanmean(np.nanmean(var_spatial_mean,axis=2),axis=1)
years_with_data = np.isfinite(var_spatial_mean_mean)

# Select only the years with data in all 12 months
var_spatial_mean = var_spatial_mean[years_with_data,:,:]
var_global_mean  = var_global_mean[years_with_data]
age = age[years_with_data]
"""

#%% CALCULATIONS 4

# Calculate lat and lon bounds
lat_bounds,lon_bounds = utils.bounding_latlon(lat,lon)

# Format the variables
var_spatial_mean = np.expand_dims(var_spatial_mean,axis=0)
var_global_mean  = np.expand_dims(var_global_mean, axis=0)

# The ensemble members for ERA5 have a different grid with half resolution. Since this is difficult to account for,
# only the means is used. The mean and the ensemble member are set to be the same.
var_spatial_members = np.expand_dims(var_spatial_mean,axis=1)
var_global_members  = np.expand_dims(var_global_mean, axis=1)

# Get other metadata
methods = ['ERA5']
ens_spatial = np.array([1])
ens_global = ens_spatial

# If this data can't be reformatted to the standard format, add a note here 
notes = ['ERA5 has ensemble members on a different grid. they are not included here. The mean and the ensemble member are set to be the same.']

# Check the shape of the variables
print(var_spatial_members.shape)
print(var_spatial_mean.shape)
print(var_global_members.shape)
print(var_global_mean.shape)


#%% FORMAT DATA

# Create new array
data_xarray_output = xr.Dataset(
    {
        var_txt+'_global_mean':    (['method','age'],                          var_global_mean,    {'units':unit_txt}),
        var_txt+'_global_members': (['method','ens_global','age'],             var_global_members, {'units':unit_txt}),
        var_txt+'_spatial_mean':   (['method','age','lat','lon'],              var_spatial_mean,   {'units':unit_txt}),
        var_txt+'_spatial_members':(['method','ens_spatial','age','lat','lon'],var_spatial_members,{'units':unit_txt})
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
        'dataset_name':      'ERA5',
        'dataset_source_url':'https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels-monthly-means?tab=overview',
    },
)


#%% SAVE DATA

# Save new array
output_dir = '/projects/pd_lab/data/paleoclimate_reconstructions/presto_format/'
data_xarray_output.to_netcdf(output_dir+'era5_v1_0_0_'+var_txt+'_'+quantity_txt+'.nc')

