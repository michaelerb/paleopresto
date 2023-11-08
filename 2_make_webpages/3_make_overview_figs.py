#==============================================================================
# Make some overview figures
#    author: Michael P. Erb
#    date  : 11/7/2023
#==============================================================================

import sys
import numpy as np
#import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
import xarray as xr
import copy

save_instead_of_plot = True

# Choose dataset and version
#dataset_txt = 'daholocene';  version_txt = '1_0_0'  # Global: ensemble members; Spatial: ensemble members
#dataset_txt = 'lgmr';        version_txt = '1_0_0'  # Global: ensemble members; Spatial: ensemble members
#dataset_txt = 'kaufman2020'; version_txt = '1_0_0'  # Global: ensemble members; Spatial: ensemble members
#dataset_txt = 'lmr';         version_txt = '2_0_0'  # Global: ensemble members; Spatial: mean
#dataset_txt = 'lmr';         version_txt = '2_1_0'  # Global: ensemble members; Spatial: mean
#dataset_txt = 'phyda';       version_txt = '1_0_0'  # Global: p5 p50 p95; Spatial: p5 p50 p95
#dataset_txt = 'neukom2019';  version_txt = '1_0_0'  # Global: ensemble members; Spatial: ensemble members
#dataset_txt = 'era20c';      version_txt = '1_0_0'  # Global: ensemble members; Spatial: ensemble members
#dataset_txt = 'era5';        version_txt = '1_0_0'  # Not created yet
dataset_txt = sys.argv[1];   version_txt = sys.argv[2]


#%% SET SPECIFIC PARAMETERS

# Set some parameters by dataset
if   dataset_txt == 'daholocene':  x_range = [12000,0];   y_range = [-3,1.25]; maxval = 5; ref_txt = '0-1 ka';       line_weight = 3; fig_type = 'contourf'
elif dataset_txt == 'lgmr':        x_range = [12000,0];   y_range = [-3,1.25]; maxval = 5; ref_txt = '0-1 ka';       line_weight = 3; fig_type = 'contourf'
elif dataset_txt == 'kaufman2020': x_range = [12000,0];   y_range = [-3,1.25]; maxval = 5; ref_txt = '0-1 ka';       line_weight = 3; fig_type = 'pcolormesh'
elif dataset_txt == 'lmr':         x_range = [0,2000];    y_range = [-.5,.85]; maxval = 2; ref_txt = '0-1 ka';       line_weight = 1; fig_type = 'contourf'
elif dataset_txt == 'phyda':       x_range = [0,2000];    y_range = [-.5,.85]; maxval = 2; ref_txt = '0-1 ka';       line_weight = 1; fig_type = 'contourf'
elif dataset_txt == 'neukom2019':  x_range = [0,2000];    y_range = [-.5,.85]; maxval = 2; ref_txt = '0-1 ka';       line_weight = 1; fig_type = 'contourf'
elif dataset_txt == 'era20c':      x_range = [1900,2021]; y_range = [-.5,1];   maxval = 2; ref_txt = '1951-1980 CE'; line_weight = 1; fig_type = 'contourf'
elif dataset_txt == 'era5':        x_range = [1900,2021]; y_range = [-.5,1];   maxval = 2; ref_txt = '1951-1980 CE'; line_weight = 1; fig_type = 'contourf'


#%% LOAD DATA

data_dir   = '/projects/pd_lab/data/paleoclimate_reconstructions/presto_format/'
filename_txt = dataset_txt+'_v'+version_txt+'_tas_annual'
print('===== LOADING '+filename_txt+' =====')

data_xarray = xr.open_dataset(data_dir+filename_txt+'.nc')
tas_spatial_mean   = data_xarray['tas_spatial_mean'].values
tas_global_mean    = data_xarray['tas_global_mean'].values
tas_global_members = data_xarray['tas_global_members'].values
ens_global         = data_xarray['ens_global'].values
method             = data_xarray['method'].values
age                = data_xarray['age'].values
lat                = data_xarray['lat'].values
lon                = data_xarray['lon'].values
lat_bounds         = data_xarray['lat_bounds'].values
lon_bounds         = data_xarray['lon_bounds'].values
dataset_name       = data_xarray.attrs['dataset_name']
data_xarray.close()
year = 1950-age


#%% Set more parameters

# Set time values
if min(year) > -2000: time_txt = 'Year (CE)';   time_var = year; time_start = min(time_var); time_end = max(time_var)
else:                 time_txt = 'Age (yr BP)'; time_var = age;  time_start = max(time_var); time_end = min(time_var)

# Set other properties
levels = np.arange(-maxval,maxval+.1,maxval/10)
colormap = plt.colormaps['bwr']
norm = BoundaryNorm(levels,ncolors=colormap.N,clip=True)


#%% Process data

# Remove the chosen reference period from the reconstruction
if   ref_txt == '0-1 ka':       ind_ref = np.where((age  >= 0)    & (age  < 1000))[0]
elif ref_txt == '1951-1980 CE': ind_ref = np.where((year >= 1951) & (year <= 1980))[0]

# Process spatial mean
tas_spatial_mean = tas_spatial_mean - np.nanmean(tas_spatial_mean[:,ind_ref,:,:],axis=1)[:,None,:,:]
tas_spatial_mean_allmethods = np.nanmean(tas_spatial_mean,axis=0)

# Get dimensions
n_methods    = tas_global_members.shape[0]
n_ens_global = tas_global_members.shape[1]
n_time       = tas_global_members.shape[2]

# Remove a reference period
tas_ref_period_mean = np.nanmean(tas_global_mean[:,ind_ref],axis=1)
tas_global_members = tas_global_members - tas_ref_period_mean[:,None,None]
tas_global_mean    = tas_global_mean    - tas_ref_period_mean[:,None]

# Calculate global mean
gmt_mean = np.nanmean(tas_global_mean,axis=0)

# Process global mean bounds (PHYDA has percentiles, the others have ensemble members)
if dataset_txt != 'phyda':
    tas_global_members_allmethods = np.reshape(tas_global_members,(n_methods*n_ens_global,n_time))
    gmt_lowerbound = np.percentile(tas_global_members_allmethods,2.5, axis=0)
    gmt_upperbound = np.percentile(tas_global_members_allmethods,97.5,axis=0)
    global_uncertainty_txt = '2.5 - 97.5th percentile range'
    #
elif dataset_txt == 'phyda':
    #print(ens_global)
    gmt_lowerbound = tas_global_members[0,0,:]
    gmt_upperbound = tas_global_members[0,2,:]
    global_uncertainty_txt = '5 - 95th percentile range'

# Compute bounds of the time var
time_bounds = copy.deepcopy(time_var)
time_bounds = (time_var[:-1] + time_var[1:])/2
time_diff = time_bounds[1] - time_bounds[0]
time_bounds = np.insert(time_bounds,0,time_bounds[0]-time_diff)
time_bounds = np.append(time_bounds,time_bounds[-1]+time_diff)


#%% FIGURES

# Make an overview figure
plt.figure(figsize=(20,14))
ax1 = plt.subplot2grid((1,1),(0,0))

# Zonal mean temperature
if fig_type == 'contourf':
    hov1 = ax1.contourf(time_var,lat,np.transpose(np.mean(tas_spatial_mean_allmethods,axis=2)),cmap=colormap,levels=levels,extend='both')
    colorbar1 = plt.colorbar(hov1,orientation='horizontal',ax=ax1,fraction=0.08,pad=0.09)
elif fig_type == 'pcolormesh':
    hov1 = ax1.pcolormesh(time_bounds,lat_bounds,np.transpose(np.mean(tas_spatial_mean_allmethods,axis=2)),cmap=colormap,norm=norm)
    colorbar1 = plt.colorbar(hov1,orientation='horizontal',ax=ax1,fraction=0.08,pad=0.09)
ax1.set_xlim(x_range)
ax1.set_ylim(-90,90)
ax1.set_xlabel(time_txt,fontsize=20)
ax1.set_ylabel('Latitude ($^\circ$)',fontsize=20)
ax1.tick_params(axis='both',which='major',labelsize=20)
colorbar1.ax.tick_params(labelsize=20)
colorbar1.set_label('Zonal mean $\Delta$T ($^\circ$C)',fontsize=20)
ax1.set_title('$\Delta$T in '+dataset_name+' v'+str(version_txt).replace('_','.'),fontsize=32)

# Global mean temperature
ax_twin1 = ax1.twinx()
line1, = ax_twin1.plot(time_var,gmt_mean,color='k',linewidth=line_weight,label='Mean')
range1 = ax_twin1.fill_between(time_var,gmt_lowerbound,gmt_upperbound,facecolor='k',alpha=0.15,label='95% range')
ax_twin1.axhline(y=0,color='k',linewidth=1,linestyle='dashed',alpha=0.5)
ax_twin1.set_ylim(y_range)
ax_twin1.set_ylabel('Global mean $\Delta$T ($^\circ$C) and '+global_uncertainty_txt,fontsize=20)
ax_twin1.tick_params(axis='both',which='major',labelsize=20)

ax1.yaxis.tick_right(); ax1.yaxis.set_label_position('right')
ax_twin1.yaxis.tick_left(); ax_twin1.yaxis.set_label_position('left')

if save_instead_of_plot:
    plt.savefig('figures/overview_'+filename_txt+'.png',dpi=150,format='png',bbox_inches='tight',pad_inches=0.0)
    plt.close()
else:
    plt.show()

