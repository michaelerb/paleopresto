#==============================================================================
# Make a standardized netcdf file for the Kaufman et al., 2020 composites.
#    author: Michael P. Erb
#    date  : 5/4/2023
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
tas_ens = np.zeros((n_methods,n_ens,n_ages,n_lat,n_lon)); tas_ens[:] = np.nan
tas_ens[0,:,:,:,0] = np.swapaxes(data_xarray['scc_latbands'].values,0,2)
tas_ens[1,:,:,:,0] = np.swapaxes(data_xarray['dcc_latbands'].values,0,2)
tas_ens[2,:,:,:,0] = np.swapaxes(data_xarray['cps_latbands'].values,0,2)
tas_ens[3,:,:,:,0] = np.swapaxes(data_xarray['pai_latbands'].values,0,2)
tas_ens[4,:,:,:,0] = np.swapaxes(data_xarray['gam_latbands'].values,0,2)

# Create a single spatial variable with the chosen dimensions
tas_global = np.zeros((n_methods,n_ens,n_ages)); tas_global[:] = np.nan
tas_global[0,:,:] = np.swapaxes(data_xarray['scc_globalmean'].values,0,1)
tas_global[1,:,:] = np.swapaxes(data_xarray['dcc_globalmean'].values,0,1)
tas_global[2,:,:] = np.swapaxes(data_xarray['cps_globalmean'].values,0,1)
tas_global[3,:,:] = np.swapaxes(data_xarray['pai_globalmean'].values,0,1)
tas_global[4,:,:] = np.swapaxes(data_xarray['gam_globalmean'].values,0,1)

# Create an average across ensemble members
tas_mean = np.mean(tas_ens,axis=1)

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
output_name = 'kaufman2020_v1_0_0_tas_annual.nc'
data_xarray_output.to_netcdf(output_dir+output_name)
