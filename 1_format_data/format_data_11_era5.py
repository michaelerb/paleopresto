#==============================================================================
# Make a standardized netcdf file for the ERA5 reanalysis.
#    author: Michael P. Erb
#    date  : 9/12/2023
#=============================================================================

import sys
import numpy as np
import xarray as xr
import utils
import calendar


#%% SETTINGS

var_txt = 'tas'
#var_txt = 'precip'
#var_txt = 'slp'
#var_txt = 'u10'
#var_txt = 'v10'
quantity_txt = 'annual'
#quantity_txt = 'jja'
#quantity_txt = 'djf'
#var_txt = sys.argv[1]; quantity_txt = sys.argv[2]


#%% LOAD DATA
# Note: data can be downloaded at: https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels-monthly-means?tab=overview

if   var_txt == 'tas':    file_var_txt = 't2m';         load_var_txt = 't2m'
elif var_txt == 'precip': file_var_txt = 'totalprecip'; load_var_txt = 'tp'
elif var_txt == 'slp':    file_var_txt = 'msl';         load_var_txt = 'msl'
elif var_txt == 'u10':    file_var_txt = 'u10';         load_var_txt = 'u10'
elif var_txt == 'v10':    file_var_txt = 'v10';         load_var_txt = 'v10'

print('=== Processing ERA5 ===')

# Load data
data_dir = '/projects/pd_lab/data/modern_datasets/ERA5/'
data_xarray = xr.open_dataset(data_dir+'era5_monthly_'+file_var_txt+'.nc')
var  = data_xarray[load_var_txt].values[:,0,:,:]
lon  = data_xarray['longitude'].values
lat  = data_xarray['latitude'].values
time = data_xarray['time'].values
data_xarray.close()

data_xarray = xr.open_dataset(data_dir+'era5_monthly_'+file_var_txt+'_ens.nc')
var_ens = data_xarray[load_var_txt].values[:,0,:,:,:]
lon_ens = data_xarray['longitude'].values
lat_ens = data_xarray['latitude'].values

years = np.arange(1959,2023,1)
age = 1950-years


#%% CALCULATIONS 1

# Change units
if   var_txt == 'tas':    pass # Units: K
elif var_txt == 'precip': var = var*1000; var_ens = var_ens*1000  # Convert precipitation from m/day to mm/day
elif var_txt == 'slp':    var = var/100;  var_ens = var_ens/100   # Convert SLP from Pa to hPa
elif var_txt == 'u10':    pass # Units: m/s
elif var_txt == 'v10':    pass # Units: m/s

# Compute the global mean
lat_weights = np.cos(np.deg2rad(data_xarray.latitude))
var_global = data_xarray[load_var_txt].weighted(lat_weights).mean(('longitude','latitude')).values[:,0,:]
data_xarray.close()

# Put the var_ens data on the same grid as var
ind_lat = np.where(np.in1d(lat,lat_ens))[0]
ind_lon = np.where(np.in1d(lon,lon_ens))[0]
ntime = var.shape[0]
nlat  = var.shape[1]
nlon  = var.shape[2]
n_ens = var_ens.shape[1]

var_ens_new = np.zeros((ntime,n_ens,nlat,nlon)); var_ens_new[:] = np.nan
for counter,i in enumerate(ind_lat):
    #print(lat[i],lat_ens[counter])
    var_ens_new[:,:,i,ind_lon] = var_ens[:,:,counter,:]


#%% CALCULATIONS 2

# Get number of days in each month
month_lengths = data_xarray.time.dt.days_in_month.values

# Reshape the data
nyears = int(ntime/12)
var_reshape        = np.reshape(var,(nyears,12,nlat,nlon)) 
var_ens_reshape    = np.reshape(var_ens_new,(nyears,12,n_ens,nlat,nlon)) 
var_global_reshape = np.reshape(var_global,(nyears,12,n_ens)) 
ndays_reshape      = np.reshape(month_lengths,(nyears,12))

# Compute the means weighted by the correct number of days in each month.
var_season        = np.zeros((nyears,nlat,nlon));       var_season[:]        = np.nan
var_ens_season    = np.zeros((nyears,n_ens,nlat,nlon)); var_ens_season[:]    = np.nan
var_global_season = np.zeros((nyears,n_ens));           var_global_season[:] = np.nan
i=0;year=years[i]
for i,year in enumerate(years):
    #
    days_in_months = ndays_reshape[i,:]
    #
    # Annual-mean
    if quantity_txt == 'annual':
        var_season[i,:,:]       = np.average(var_reshape[i,:,:,:],      axis=0,weights=days_in_months)
        var_ens_season[i,:,:,:] = np.average(var_ens_reshape[i,:,:,:,:],axis=0,weights=days_in_months)
        var_global_season[i,:]  = np.average(var_global_reshape[i,:,:], axis=0,weights=days_in_months)
    #
    # JJA-mean
    if quantity_txt == 'jja':
        var_season[i,:,:]       = np.average(var_reshape[i,[5,6,7],:,:],      axis=0,weights=days_in_months[[5,6,7]])
        var_ens_season[i,:,:,:] = np.average(var_ens_reshape[i,[5,6,7],:,:,:],axis=0,weights=days_in_months[[5,6,7]])
        var_global_season[i,:]  = np.average(var_global_reshape[i,[5,6,7],:], axis=0,weights=days_in_months[[5,6,7]])
    #
    # DJF-mean
    if quantity_txt == 'djf':
        # Get DJF weights
        try: weights_DJF = np.concatenate((np.expand_dims(days_in_months[11],axis=0),ndays_reshape[i+1,:2]),axis=0)
        except:
            if calendar.isleap(year+1): weights_DJF = np.array([31,31,29])
            else:                       weights_DJF = np.array([31,31,28])
        #
        data_D = np.expand_dims(var_reshape[i,11,:,:],axis=0)
        try:    data_JF = var_reshape[i+1,[0,1],:,:]
        except: data_JF = np.zeros((2,nlat,nlon)); data_JF[:] = np.nan
        data_DJF = np.concatenate((data_D,data_JF),axis=0)
        var_season[i,:,:] = np.average(data_DJF,axis=0,weights=weights_DJF)
        #
        data_ens_D = np.expand_dims(var_ens_reshape[i,11,:,:,:],axis=0)
        try:    data_ens_JF = var_ens_reshape[i+1,[0,1],:,:,:]
        except: data_ens_JF = np.zeros((2,n_ens,nlat,nlon)); data_ens_JF[:] = np.nan
        data_ens_DJF = np.concatenate((data_ens_D,data_ens_JF),axis=0)
        var_ens_season[i,:,:,:] = np.average(data_ens_DJF,axis=0,weights=weights_DJF)
        #
        data_global_D = np.expand_dims(var_global_reshape[i,11,:],axis=0)
        try:    data_global_JF = var_global_reshape[i+1,[0,1],:]
        except: data_global_JF = np.zeros((2,n_ens)); data_global_JF[:] = np.nan
        data_global_DJF = np.concatenate((data_global_D,data_global_JF),axis=0)
        var_global_season[i,:] = np.average(data_global_DJF,axis=0,weights=weights_DJF)


#%% CALCULATIONS 3

# Find the years with data in all 12 months
var_season_mean = np.nanmean(np.nanmean(var_season,axis=2),axis=1)
years_with_data = np.isfinite(var_season_mean)


#%% CALCULATIONS 4

# Calculate lat and lon bounds
lat_bounds,lon_bounds = utils.bounding_latlon(lat,lon)

# Format the variables
var_ens_season    = np.swapaxes(var_ens_season,   0,1)
var_global_season = np.swapaxes(var_global_season,0,1)
var_ens_season    = np.expand_dims(var_ens_season,   axis=0)
var_season        = np.expand_dims(var_season,       axis=0)
var_global_season = np.expand_dims(var_global_season,axis=0)

# Get other metadata
methods = ['ERA-20C']
ens_spatial = np.arange(n_ens)+1
ens_global = ens_spatial

# If this data can't be reformatted to the standard format, add a note here 
notes = ['']


#%% FORMAT DATA

# Create new array
data_xarray_output = xr.Dataset(
    {
        var_txt+'_global':(['method','ens_global','age'],             var_global_season[:,:,years_with_data]),
        var_txt+'_mean':  (['method','age','lat','lon'],              var_season[:,years_with_data,:,:]),
        var_txt+'_ens':   (['method','ens_spatial','age','lat','lon'],var_ens_season[:,:,years_with_data,:,:])
        #var_txt+'_global_mean':    (['method','age'],                          var_global_mean,    {'units':'degrees Celsius'}),
        #var_txt+'_global_members': (['method','ens_global','age'],             var_global_members, {'units':'degrees Celsius'}),
        #var_txt+'_spatial_mean':   (['method','age','lat','lon'],              var_spatial_mean,   {'units':'degrees Celsius'}),
        #var_txt+'_spatial_members':(['method','ens_spatial','age','lat','lon'],var_spatial_members,{'units':'degrees Celsius'})
    },
    coords={
        'method':     (['method'],methods),
        'notes':      (['notes'],notes),
        'ens_global': (['ens_global'],ens_global,{'description':'ensemble members'}),
        'ens_spatial':(['ens_spatial'],ens_spatial,{'description':'ensemble members'}),
        'age':        (['age'],age[years_with_data],{'units':'yr BP'}),
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

