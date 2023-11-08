#==============================================================================
# Make a standardized netCDF file for PHYDA.
#    author: Michael P. Erb
#    date  : 9/12/2023
#==============================================================================

import sys
import numpy as np
import xarray as xr
import utils


#%% SETTINGS

#var_txt = 'tas'
#var_txt = 'pdsi'
#var_txt = 'spei'
#season_txt = 'annual'
#season_txt = 'jja'
#season_txt = 'djf'
var_txt = sys.argv[1]; season_txt = sys.argv[2]


#%% LOAD DATA
# Note: data can be downloaded at: https://www.ncei.noaa.gov/access/paleo-search/study/33112

print('=== Processing PHYDA ===')
data_dir = '/projects/pd_lab/data/paleoclimate_reconstructions/PHYDA/'
if   season_txt == 'annual': data_file = 'da_hydro_AprMar_r.1-2000_d.05-Jan-2018.nc'
elif season_txt == 'jja':    data_file = 'da_hydro_JunAug_r.1-2000_d.05-Jan-2018.nc'
elif season_txt == 'djf':    data_file = 'da_hydro_DecFeb_r.1-2000_d.05-Jan-2018.nc'

# Load data
data_xarray = xr.open_dataset(data_dir+data_file)
years = data_xarray['time'].values
lat   = data_xarray['lat'].values
lon   = data_xarray['lon'].values
lat_weights = np.cos(np.deg2rad(data_xarray.lat))

if var_txt == 'tas':
    var_spatial_mean    = data_xarray['tas_mn'].values
    var_spatial_members = data_xarray['tas_pc'].values  # The 5th, 50th, and 95th percentiles
    var_global_mean     = data_xarray['gmt_mn'].values
    var_global_members  = data_xarray['gmt_pc'].values
    var_units           = data_xarray['tas_mn'].units
elif var_txt == 'pdsi':
    var_spatial_mean    = data_xarray['pdsi_mn'].values
    var_spatial_members = data_xarray['pdsi_pc'].values
    var_global_mean     = data_xarray.pdsi_mn.weighted(lat_weights).mean(('lon','lat'))
    var_global_members  = data_xarray.pdsi_pc.weighted(lat_weights).mean(('lon','lat'))  #TODO: Check this.
    var_units           = data_xarray['pdsi_mn'].units
elif var_txt == 'spei':
    var_spatial_mean    = data_xarray['spei_mn'].values
    var_spatial_members = data_xarray['spei_pc'].values
    var_global_mean     = data_xarray.spei_mn.weighted(lat_weights).mean(('lon','lat'))
    var_global_members  = data_xarray.spei_pc.weighted(lat_weights).mean(('lon','lat'))  #TODO: Check this.
    var_units           = data_xarray['spei_mn'].units

data_xarray.close()

age = 1950 - years


"""
#%% FIGURE
import matplotlib.pyplot as plt

# Global mean temperature
plt.figure(figsize=(20,14))
ax1 = plt.subplot2grid((1,1),(0,0))
#line1, = ax1.plot(years,var_global,color='k',linewidth=3,label='Mean')
#range1 = ax1.fill_between(years,var_global_members[0,:],tas_global_pc[2,:],facecolor='k',alpha=0.15,label='95% range')
line1, = ax1.plot(years,var_global_mean,color='k',linewidth=3,label='Mean')
ax1.axhline(y=0,color='k',linewidth=1,linestyle='dashed',alpha=0.5)
ax1.set_ylabel('Global mean $\Delta$T ($^\circ$C)',fontsize=20)
ax1.tick_params(axis='both',which='major',labelsize=20)
plt.show()
"""

#%% CALCULATIONS

# Calculate lat and lon bounds
lat_bounds,lon_bounds = utils.bounding_latlon(lat,lon)

# Format the variables
var_spatial_mean    = np.expand_dims(var_spatial_mean,   axis=0)
var_spatial_members = np.expand_dims(var_spatial_members,axis=0)
var_global_mean     = np.expand_dims(var_global_mean,    axis=0)
var_global_members  = np.expand_dims(var_global_members, axis=0)

# Set other metadata
methods = ['PHYDA']
ens_global  = ['p5','p50','p95']
ens_spatial = ['p5','p50','p95']

# If this data can't be reformatted to the standard format, add a note here 
notes = ['Ensemble members are actually the 5th, 50th, and 95th percentiles']

# Check the shape of the variables
print(var_spatial_members.shape)
print(var_spatial_mean.shape)
print(var_global_members.shape)
print(var_global_mean.shape)


#%% FORMAT DATA

# Create new arrays
data_xarray_output = xr.Dataset(
    {
        var_txt+'_global_mean':    (['method','age'],                          var_global_mean,    {'units':var_units}),
        var_txt+'_global_members': (['method','ens_global','age'],             var_global_members, {'units':var_units}),
        var_txt+'_spatial_mean':   (['method','age','lat','lon'],              var_spatial_mean,   {'units':var_units}),
        var_txt+'_spatial_members':(['method','ens_spatial','age','lat','lon'],var_spatial_members,{'units':var_units})
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
        'dataset_name':      'PHYDA',
        'dataset_source_url':'https://www.ncei.noaa.gov/access/paleo-search/study/33112',
    },
)


#%% SAVE DATA

# Save new array
output_dir = '/projects/pd_lab/data/paleoclimate_reconstructions/presto_format/'
data_xarray_output.to_netcdf(output_dir+'phyda_v1_0_0_'+var_txt+'_'+season_txt+'.nc')

