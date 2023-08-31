#==============================================================================
# Make some overview figures
#    author: Michael P. Erb
#    date  : 5/5/2023
#==============================================================================

import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
import xarray as xr
import copy

save_instead_of_plot = True

# Choose dataset and version
#dataset_txt = 'daholocene';  version_txt = '1_0_0'
#dataset_txt = 'lgmr';        version_txt = '1_0_0'
#dataset_txt = 'kaufman2020'; version_txt = '1_0_0'
#dataset_txt = 'lmr';         version_txt = '2_0_0'
#dataset_txt = 'lmr';         version_txt = '2_1_0'
#dataset_txt = 'neukom2019';  version_txt = '1_0_0'
#dataset_txt = 'era20c';      version_txt = '1_0_0'
#dataset_txt = 'era5';        version_txt = '1_0_0'
#dataset_txt = 'phyda';       version_txt = '1_0_0'
dataset_txt = sys.argv[1];   version_txt = sys.argv[2]


#%% LOAD DATA

data_dir   = '/projects/pd_lab/data/paleoclimate_reconstructions/presto_format/'
filename_txt = dataset_txt+'_v'+version_txt+'_tas_annual'
print('===== LOADING '+filename_txt+' =====')

handle = xr.open_dataset(data_dir+filename_txt+'.nc')
tas_global = handle['tas_global'].values
tas_mean   = handle['tas_mean'].values
method     = handle['method'].values
age        = handle['age'].values
lat        = handle['lat'].values
lon        = handle['lon'].values
lat_bounds = handle['lat_bounds'].values
lon_bounds = handle['lon_bounds'].values
method     = handle['method'].values
handle.close()
year = 1950-age


#%% Set parameters

# Set some parameters by dataset
if   dataset_txt == 'daholocene':  data_txt = 'Holocene Reconstruction'; x_range = [12000,0];   y_range = [-3,1.25];  maxval = 5; time_var = age;  time_txt = 'Age (yr BP)'; ref_txt = '0-1 ka';       line_weight = 3; fig_type = 'contourf'
elif dataset_txt == 'lgmr':        data_txt = 'LGMR';                    x_range = [12000,0];   y_range = [-3,1.25];  maxval = 5; time_var = age;  time_txt = 'Age (yr BP)'; ref_txt = '0-1 ka';       line_weight = 3; fig_type = 'contourf'
elif dataset_txt == 'kaufman2020': data_txt = 'Kaufman et al., 2020';    x_range = [12000,0];   y_range = [-3,1.25];  maxval = 5; time_var = age;  time_txt = 'Age (yr BP)'; ref_txt = '0-1 ka';       line_weight = 3; fig_type = 'pcolormesh'
elif dataset_txt == 'lmr':         data_txt = 'LMR';                     x_range = [0,2000];    y_range = [-.5,.85];  maxval = 2; time_var = year; time_txt = 'Year (CE)';   ref_txt = '0-1 ka';       line_weight = 1; fig_type = 'contourf'
elif dataset_txt == 'phyda':       data_txt = 'PHYDA';                   x_range = [0,2000];    y_range = [-.5,.85];  maxval = 2; time_var = year; time_txt = 'Year (CE)';   ref_txt = '0-1 ka';       line_weight = 1; fig_type = 'contourf'
elif dataset_txt == 'neukom2019':  data_txt = 'Neukom et al., 2019';     x_range = [0,2000];    y_range = [-.5,.85];  maxval = 2; time_var = year; time_txt = 'Year (CE)';   ref_txt = '0-1 ka';       line_weight = 1; fig_type = 'contourf'
elif dataset_txt == 'era20c':      data_txt = 'ERA-20C';                 x_range = [1900,2021]; y_range = [-.75,.75]; maxval = 2; time_var = year; time_txt = 'Year (CE)';   ref_txt = '1951-1980 CE'; line_weight = 1; fig_type = 'contourf'
elif dataset_txt == 'era5':        data_txt = 'ERA5';                    x_range = [1900,2021]; y_range = [-.75,.75]; maxval = 2; time_var = year; time_txt = 'Year (CE)';   ref_txt = '1951-1980 CE'; line_weight = 1; fig_type = 'contourf'

levels = np.arange(-maxval,maxval+.1,maxval/10)

colormap = plt.colormaps['bwr']
norm = BoundaryNorm(levels,ncolors=colormap.N,clip=True)


#%% Process data

# Remove the chosen reference period from the reconstruction
if   ref_txt == '0-1 ka':       ind_ref = np.where((age >= 0)     & (age < 1000))[0]
elif ref_txt == '1951-1980 CE': ind_ref = np.where((year >= 1951) & (age <= 1980))[0]
tas_mean   = tas_mean   - np.nanmean(tas_mean[:,ind_ref,:,:],axis=1)[:,None,:,:]
tas_global = tas_global - np.nanmean(np.nanmean(tas_global[:,:,ind_ref],axis=2),axis=1)[:,None,None]

# Compute the mean of all methods
tas_mean_allmethods = np.nanmean(tas_mean,axis=0)
tas_global_allmethods = np.reshape(tas_global,(tas_global.shape[0]*tas_global.shape[1],tas_global.shape[2]))

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
    hov1 = ax1.contourf(time_var,lat,np.transpose(np.mean(tas_mean_allmethods,axis=2)),cmap=colormap,levels=levels,extend='both')
    colorbar1 = plt.colorbar(hov1,orientation='horizontal',ax=ax1,fraction=0.08,pad=0.09)
elif fig_type == 'pcolormesh':
    hov1 = ax1.pcolormesh(time_bounds,lat_bounds,np.transpose(np.mean(tas_mean_allmethods,axis=2)),cmap=colormap,norm=norm)
    colorbar1 = plt.colorbar(hov1,orientation='horizontal',ax=ax1,fraction=0.08,pad=0.09)
ax1.set_xlim(x_range)
ax1.set_ylim(-90,90)
ax1.set_xlabel(time_txt,fontsize=20)
ax1.set_ylabel('Latitude ($^\circ$)',fontsize=20)
ax1.tick_params(axis='both',which='major',labelsize=20)
colorbar1.ax.tick_params(labelsize=20)
colorbar1.set_label('Zonal mean $\Delta$T ($^\circ$C)',fontsize=20)
ax1.set_title('$\Delta$T in '+data_txt+' v'+str(version_txt).replace('_','.'),fontsize=32)

# Global mean temperature
ax_twin1 = ax1.twinx()
line1, = ax_twin1.plot(time_var,np.nanmean(tas_global_allmethods,axis=0),color='k',linewidth=line_weight,label='Mean')
range1 = ax_twin1.fill_between(time_var,np.nanpercentile(tas_global_allmethods,2.5,axis=0),np.nanpercentile(tas_global_allmethods,97.5,axis=0),facecolor='k',alpha=0.15,label='95% range')
ax_twin1.axhline(y=0,color='k',linewidth=1,linestyle='dashed',alpha=0.5)
ax_twin1.set_ylim(y_range)
ax_twin1.set_ylabel('Global mean $\Delta$T ($^\circ$C)',fontsize=20)
ax_twin1.tick_params(axis='both',which='major',labelsize=20)

ax1.yaxis.tick_right(); ax1.yaxis.set_label_position('right')
ax_twin1.yaxis.tick_left(); ax_twin1.yaxis.set_label_position('left')

if save_instead_of_plot:
    plt.savefig('overview_'+filename_txt+'.png',dpi=150,format='png',bbox_inches='tight',pad_inches=0.0)
    plt.close()
else:
    plt.show()

