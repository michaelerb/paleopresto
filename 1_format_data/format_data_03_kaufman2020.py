#==============================================================================
# Make a standardized netcdf file for the Kaufman et al., 2020 composites.
#    author: Michael P. Erb
#    date  : 9/12/2023
#==============================================================================

import numpy as np
import xarray as xr


#%% LOAD DATA
# Note: data can be downloaded at: https://www.ncei.noaa.gov/access/paleo-search/study/29712

print('=== Processing Kaufman et al., 2020 composites ===')

# Load data
data_dir = '/projects/pd_lab/data/paleoclimate_reconstructions/Kaufman2020/'
data_xarray = xr.open_dataset(data_dir+'temp12k_alldata.nc')

# Set lats and lons
print(data_xarray.latband_ranges.values)
lat = np.arange(75,-90,-30)
lon = np.array([0])
lat_bounds = np.arange(90,-91,-30)
lon_bounds = np.array([-180,180])

# Get other metadata
methods = ['SCC','DCC','CPS','PAI','GAM']
ens_spatial = data_xarray.ens.values
ens_global  = ens_spatial
age = data_xarray.age.values


#%% CALCULATIONS

# Get dimensions
n_methods = len(methods)
n_ens     = len(ens_spatial)
n_ages    = len(age)
n_lat     = len(lat)
n_lon     = len(lon)

# Create a single spatial variable with the chosen dimensions
var_spatial_members = np.zeros((n_methods,n_ens,n_ages,n_lat,n_lon)); var_spatial_members[:] = np.nan
var_spatial_members[0,:,:,:,0] = np.swapaxes(data_xarray['scc_latbands'].values,0,2)
var_spatial_members[1,:,:,:,0] = np.swapaxes(data_xarray['dcc_latbands'].values,0,2)
var_spatial_members[2,:,:,:,0] = np.swapaxes(data_xarray['cps_latbands'].values,0,2)
var_spatial_members[3,:,:,:,0] = np.swapaxes(data_xarray['pai_latbands'].values,0,2)
var_spatial_members[4,:,:,:,0] = np.swapaxes(data_xarray['gam_latbands'].values,0,2)

# Create a single spatial variable with the chosen dimensions
var_global_members = np.zeros((n_methods,n_ens,n_ages)); var_global_members[:] = np.nan
var_global_members[0,:,:] = np.swapaxes(data_xarray['scc_globalmean'].values,0,1)
var_global_members[1,:,:] = np.swapaxes(data_xarray['dcc_globalmean'].values,0,1)
var_global_members[2,:,:] = np.swapaxes(data_xarray['cps_globalmean'].values,0,1)
var_global_members[3,:,:] = np.swapaxes(data_xarray['pai_globalmean'].values,0,1)
var_global_members[4,:,:] = np.swapaxes(data_xarray['gam_globalmean'].values,0,1)

# Create an average across ensemble members
var_spatial_mean = np.mean(var_spatial_members,axis=1)
var_global_mean  = np.mean(var_global_members,axis=1)

# If this data can't be reformatted to the standard format, add a note here 
notes = ['']

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
        'dataset_name':      'Kaufman et al., 2020',
        'dataset_source_url':'https://www.ncei.noaa.gov/access/paleo-search/study/29712',
    },
)


#%% SAVE DATA

# Save new array
output_dir = '/projects/pd_lab/data/paleoclimate_reconstructions/presto_format/'
data_xarray_output.to_netcdf(output_dir+'kaufman2020_v1_0_0_tas_annual.nc')

