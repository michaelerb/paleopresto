#==============================================================================
# Make a standardized netcdf file for the Holocene Hydroclimate composites.
#    author: Michael P. Erb
#    date  : 5/4/2023
#=============================================================================

import numpy as np
import xarray as xr
import glob
import regionmask


#%% SETTINGS

var_txt = 'tas'
#var_txt = 'precip'


#%% WGI REGIONS

# Get all WGI regions
ar6_all = regionmask.defined_regions.ar6.all
ar6_abbreviations = ar6_all.abbrevs
n_regions = len(ar6_abbreviations)


#%% COMPUTE AREAS

# Make a lat/lon grid
lat = np.arange(-89.5,90,1)
lon = np.arange(-179.5,180,1)

ar6_mask = regionmask.defined_regions.ar6.all.mask(lon,lat).values
lon_2d,lat_2d = np.meshgrid(lon,lat)
lon_flatten = lon_2d.flatten()
lat_flatten = lat_2d.flatten()
ar6_mask_flatten = ar6_mask.flatten()

ar6_area = np.zeros(n_regions); ar6_area[:] = np.nan
for i in range(n_regions):
    ind_selected = ar6_mask_flatten == i
    ar6_area[i] = np.sum(np.cos(np.radians(lat_flatten[ind_selected])))

"""
# Get coordinates for each region
i=0
for i in range(n_regions):
    region_shape  = ar6_all[i]._polygon
    if region_shape.geom_type == 'Polygon':
        region_lons,region_lats = region_shape.exterior.coords.xy
    elif region_shape.geom_type == 'MultiPolygon':
        for shape in region_shape:
            region_lons,region_lats = shape.exterior.coords.xy
"""


#%% SET UP

# Set ages and ens
age = np.arange(0,12001,100)
ens = np.arange(1,501)

# Get dimensions
n_methods = 1
n_ens     = len(ens)
n_ages    = len(age)


#%% LOAD DATA
# Note: data can be downloaded at: https://github.com/clhancock/HoloceneHydroclimate  #TODO: Update this with a more permanant link once the results are archived.

data_dir = '/projects/pd_lab/mpe32/HoloceneHydroclimate/Data/RegionComposites/'

# Get the names of the available files
if   var_txt == 'tas':    data_subdir = 'T/'
elif var_txt == 'precip': data_subdir = 'HC/'
filenames_all = glob.glob(data_dir+data_subdir+'*.csv')
filenames_regions = [filename.split('/')[-1].split('.')[0] for filename in filenames_all if filename[-12:] != 'byRegion.csv']

# Load regional composites
var_ens = np.zeros((n_methods,n_ens,n_ages,n_regions,1)); var_ens[:] = np.nan
for i,region_name in enumerate(ar6_abbreviations):
    if region_name not in filenames_regions: print(i,region_name,'NO FILE FOR '+var_txt+' DATA')
    else:
        var_for_region = np.loadtxt(data_dir+data_subdir+region_name+'.csv',dtype=str,delimiter=',',skiprows=1)
        var_for_region[var_for_region == 'NA'] = np.nan
        var_for_region = var_for_region.astype(float)
        var_ens[0,:,:,i,0] = np.swapaxes(var_for_region,0,1)


#%% CALCULATIONS

# Compute mean of ensemble members
var_mean = np.mean(var_ens,axis=1)

# Compute global means
var_global = np.zeros((n_methods,n_ens,n_ages)); var_global[:] = np.nan
for i in range(n_ages):
    ind_valid_hc = np.isfinite(var_mean[0,i,:,0])
    var_global[:,:,i] = np.average(var_ens[:,:,i,ind_valid_hc,0],axis=2,weights=ar6_area[ind_valid_hc])

# Get other metadata
methods = ['Holocene Hydroclimate']
ens_spatial = ens
ens_global  = ens

# If this data can't be reformatted to the standard format, add a note here 
notes = ['lat_is_ar6_region_names']


#%% FORMAT DATA

# Create new array
message_txt = np.array(['Use IPCC AR4 regions'])
data_xarray_output = xr.Dataset(
    {
        var_txt+'_global':(['method','ens_global','age'],             var_global,{'units':'Z score'}),
        var_txt+'_mean':  (['method','age','lat','lon'],              var_mean,  {'units':'Z score'}),
        var_txt+'_ens':   (['method','ens_spatial','age','lat','lon'],var_ens,   {'units':'Z score'})
    },
    coords={
        'method':     (['method'],methods),
        'notes':      (['notes'],notes),
        'ens_global': (['ens_global'],ens_global,{'description':'ensemble members'}),
        'ens_spatial':(['ens_spatial'],ens_spatial,{'description':'selected ensemble members'}),
        'age':        (['age'],age,{'units':'yr BP'}),
        'lat':        (['lat'],ar6_abbreviations,{'units':'none'}),
        'lon':        (['lon'],message_txt,{'units':'none'}),
        'lat_bounds': (['lat_bounds'],message_txt,{'units':'none'}),
        'lon_bounds': (['lon_bounds'],message_txt,{'units':'none'}),
    },
)


#%% SAVE DATA

# Save new array
output_dir = '/projects/pd_lab/data/paleoclimate_reconstructions/presto_format/'
output_name = 'holocenehydroclimate_v1_0_0_'+var_txt+'_annual.nc'
data_xarray_output.to_netcdf(output_dir+output_name)
