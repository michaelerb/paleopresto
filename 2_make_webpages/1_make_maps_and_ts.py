#==============================================================================
# Make a set of maps and time series for the visualizer webpage.
#    author: Michael P. Erb
#    date  : 5/8/2023
#==============================================================================

import sys
sys.path.append('/home/mpe32/analysis/15_Holocene_Reconstruction/data_assimilation')
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
import cartopy.util as cutil
from cartopy.feature import ShapelyFeature
import xarray as xr
import copy
import regionmask
from bokeh.plotting import figure
from bokeh.resources import CDN
from bokeh.embed import file_html
from bokeh.models import HoverTool
from bokeh.models import Range1d
from bokeh.models import Span

save_instead_of_plot = True
plt.style.use('ggplot')

# Choose dataset
#dataset_txt = 'daholocene'
#dataset_txt = 'holocenehydroclimate'
#dataset_txt = 'lgmr'
#dataset_txt = 'kaufman2020'
#dataset_txt = 'lmr'
#dataset_txt = 'neukom2019'
#dataset_txt = 'era20c'
#dataset_txt = 'era5'
dataset_txt = sys.argv[1]

# Choose version
#version_txt = '1_0_0'
#version_txt = '2_0_0'
#version_txt = '2_1_0'
version_txt = sys.argv[2]

# Choose variable
#var_txt = 'tas'
#var_txt = 'precip'
#var_txt = 'd18Op'
#var_txt = 'slp'
#var_txt = 'winds'
var_txt = sys.argv[3]

# Choose quantity
#quantity_txt = 'Annual'
#quantity_txt = 'JJA'
#quantity_txt = 'DJF'
quantity_txt = sys.argv[4]


#%% LOAD DATA

data_dir   = '/projects/pd_lab/data/paleoclimate_reconstructions/presto_format/'
output_dir = '/projects/pd_lab/mpe32/figures_presto/'
filename_txt = dataset_txt+'_v'+version_txt+'_'+var_txt+'_'+quantity_txt.lower()
print('===== LOADING '+filename_txt+' =====')

handle = xr.open_dataset(data_dir+filename_txt+'.nc')
var_mean   = handle[var_txt+'_mean']
var_ens    = handle[var_txt+'_ens']
var_global = handle[var_txt+'_global'].values
method     = handle['method'].values
age        = handle['age'].values
lat        = handle['lat'].values
lon        = handle['lon'].values
lat_bounds = handle['lat_bounds'].values
lon_bounds = handle['lon_bounds'].values
handle.close()
year = 1950-age

#TODO: Improve this later.
if dataset_txt == 'lmr':
    #
    if   version_txt == '2_0_0': version_txt_lmr = '2.0'
    elif version_txt == '2_1_0': version_txt_lmr = '2.1'
    if   var_txt == 'tas':    var_txt_lmr = 'air'
    elif var_txt == 'precip': var_txt_lmr = 'prate'
    #
    data_dir_lmr = '/projects/pd_lab/data/paleoclimate_reconstructions/LMR_final/'
    handle_stdev = xr.open_dataset(data_dir_lmr+var_txt_lmr+'_MCruns_ensemble_spread_LMRv'+version_txt_lmr+'.nc')
    var_ens_stdev_array = handle_stdev[var_txt_lmr].values
    var_ens_stdev_array = np.swapaxes(var_ens_stdev_array,0,1)
    var_ens_stdev_array = np.expand_dims(var_ens_stdev_array,axis=0)
    handle_stdev.close()
    #
    # Put in an xarray
    data_xarray_stdev = xr.Dataset(
        {
            'var_ens_stdev': (['method','ens','age','lat','lon'],var_ens_stdev_array)
        },
        coords={
            'method':    (['method'],['LMR']),
            'ens':       (['ens'],np.arange(20)+1,{'description':'standard deviations of the 20 Monte Carlo interations'}),
            'age':       (['age'],age,{'units':'yr BP'}),
            'lat':       (['lat'],lat,{'units':'degrees_north'}),
            'lon':       (['lon'],lon,{'units':'degrees_east'}),
            'lat_bounds':(['lat_bounds'],lat_bounds,{'units':'degrees_north'}),
            'lon_bounds':(['lon_bounds'],lon_bounds,{'units':'degrees_east'}),
        },
    )
    #
    var_ens_stdev = data_xarray_stdev['var_ens_stdev']
    if var_txt == 'precip': var_ens_stdev = var_ens_stdev*60*60*24 # Convert precipitation units from kg/m2/s to mm/day


#%% Set parameters
# Set some parameters by variable type
if var_txt == 'tas':
    unit_txt        = '$^\circ$C'
    html_unit_txt   = '&#176C'
    colorbar_txt    = quantity_txt+' 2m air temperature ($^\circ$C)'
    title_txt_bokeh = quantity_txt+' 2m air temperature (\u00B0C)'
    cmap            = 'bwr'
    levels          = np.array([-10,-5,-2,-1,-.5,-.2,-.1,0,.1,.2,.5,1,2,5,10])
elif var_txt == 'precip':
    unit_txt        = 'mm/day'
    html_unit_txt   = 'mm/day'
    colorbar_txt    = quantity_txt+' precipitation (mm day$^{-1}$)'
    title_txt_bokeh = quantity_txt+' precipitation (mm day\u207B\u00B9)'
    cmap            = 'BrBG'
    levels          = np.linspace(-1,1,21)
elif var_txt == 'd18Op':
    unit_txt        = 'permil'
    html_unit_txt   = 'permil'
    colorbar_txt    = quantity_txt+' d18Op (permil)'
    title_txt_bokeh = quantity_txt+' d18Op (permil)'
    cmap            = 'PiYG'
    levels          = np.array([-10,-5,-2,-1,-.5,-.2,-.1,0,.1,.2,.5,1,2,5,10])
elif var_txt == 'slp':
    unit_txt        = 'hPa'
    html_unit_txt   = 'hPa'
    colorbar_txt    = quantity_txt+' sea level pressure (hPa)'
    title_txt_bokeh = quantity_txt+' sea level pressure (hPa)'
    cmap            = 'PuOr_r'
    levels          = np.arange(-5,5.1,.5)

# Set some parameters by dataset
if dataset_txt == 'daholocene':
    dataset_name = 'Holocene Reconstruction'
    x_range      = [12000,0]
    ts_height    = 300
    ts_yrange    = [-5,1]
    #
elif dataset_txt == 'holocenehydroclimate':
    dataset_name = 'Holocene Hydroclimate'
    x_range      = [12000,0]
    ts_height    = 300
    ts_yrange    = [-5,1]
    #
elif dataset_txt == 'lgmr':
    dataset_name = 'LGMR'
    x_range      = [24000,0]
    ts_height    = 300
    if   var_txt == 'tas':   ts_yrange = [-30,5]
    elif var_txt == 'd18Op': ts_yrange = [-15,5]
    #
elif dataset_txt == 'kaufman2020':
    dataset_name = 'Kaufman et al., 2020'
    ts_height    = 500
    x_range      = [12000,0]
    ts_yrange    = [-3,3]
    #
elif dataset_txt == 'lmr':
    dataset_name = 'LMR'
    x_range      = [0,2000]
    ts_height    = 300
    if   var_txt == 'tas':    ts_yrange = [-2,3]
    elif var_txt == 'precip': ts_yrange = [-2,2]
    #
elif dataset_txt == 'neukom2019':
    dataset_name = 'Neukom et al., 2019'
    x_range      = [0,2000]
    ts_height    = 500
    ts_yrange    = [-2,3]
    #
elif dataset_txt == 'era20c':
    dataset_name = 'ERA-20C'
    x_range      = [1900,2010]
    ts_height    = 300
    ts_yrange    = [-2,3]
    #
elif dataset_txt == 'era5':
    dataset_name = 'ERA5'
    x_range      = [1959,2022]
    ts_height    = 300
    ts_yrange    = [-2,3]

# Set some variables for the holocenehydroclimate dataset
if dataset_txt == 'holocenehydroclimate':
    if var_txt == 'tas':
        unit_txt        = 'Z-score'
        html_unit_txt   = 'Z-score'
        colorbar_txt    = quantity_txt+' temperature (Z-score)'
        title_txt_bokeh = quantity_txt+' temperature (Z-score)'
        ts_yrange       = [-4,4]
    elif var_txt == 'precip':
        unit_txt        = 'Z-score'
        html_unit_txt   = 'Z-score'
        colorbar_txt    = quantity_txt+' hydroclimate (Z-score)'
        title_txt_bokeh = quantity_txt+' hydroclimate (Z-score)'
        ts_yrange       = [-4,4]

# Set some miscellaneous parameters
if   dataset_txt in ['era20c','era5']:    time_var = year; time_name_txt = 'Year'; time_unit_txt = 'CE';    ref_period_txt = '1951-1980 CE'
elif dataset_txt in ['lmr','neukom2019']: time_var = year; time_name_txt = 'Year'; time_unit_txt = 'CE';    ref_period_txt = '0-1 ka'
else:                                     time_var = age;  time_name_txt = 'Age';  time_unit_txt = 'yr BP'; ref_period_txt = '0-1 ka'
if   dataset_txt == 'kaufman2020':          map_type = 'pcolormesh';       make_gridded_ts = True;  make_regional_ts = False
elif dataset_txt == 'holocenehydroclimate': map_type = 'regions_ipcc_ar6'; make_gridded_ts = False; make_regional_ts = True
else:                                       map_type = 'contourf';         make_gridded_ts = True;  make_regional_ts = True


#%% Process data
# Remove the chosen reference period from the reconstruction
if   ref_period_txt == '0-1 ka':       ind_ref = np.where((age >= 0)     & (age < 1000))[0]
elif ref_period_txt == '1951-1980 CE': ind_ref = np.where((year >= 1951) & (age <= 1980))[0]
var_ens    = var_ens    - np.nanmean(var_mean[:,ind_ref,:,:],axis=1)[:,None,None,:,:]
var_mean   = var_mean   - np.nanmean(var_mean[:,ind_ref,:,:],axis=1)[:,None,:,:]
var_global = var_global - np.nanmean(np.nanmean(var_global[:,:,ind_ref],axis=2),axis=1)[:,None,None]

if dataset_txt == 'lmr': #TODO: Improve this later
    var_upper_2std = var_ens + (2*var_ens_stdev)
    var_lower_2std = var_ens - (2*var_ens_stdev)

# Compute the mean of all methods
var_mean_allmethods = np.nanmean(var_mean,axis=0)

# Put ensemble members an methods on the same axis
var_ens_allmethods = np.reshape(var_ens.values,(var_ens.shape[0]*var_ens.shape[1],var_ens.shape[2],var_ens.shape[3],var_ens.shape[4]))

# Compute global means
globalmean_all_ens = np.mean(var_global,axis=0)
globalmean_all     = np.mean(globalmean_all_ens,axis=0)


#%% Compute grid aspects
n_total = len(lat)*len(lon)

reduce_number_of_ts = True
if reduce_number_of_ts & (dataset_txt not in ['holocenehydroclimate','kaufman2020']):
    #
    # Print stats for original grid
    lat_grid_mean = np.abs(np.mean(lat[1:]-lat[:-1]))
    lon_grid_mean = np.abs(np.mean(lon[1:]-lon[:-1]))
    print('=== ORIGINAL GRID ===')
    print('Original grid size (degrees):')
    print('Lat:',lat_grid_mean)
    print('Lon:',lon_grid_mean)
    print('Total points:',n_total)
    #
    # Figure out skip values
    lat_grid_desired = 5
    lon_grid_desired = 5
    lat_skip_factor = max(round(lat_grid_desired/lat_grid_mean),1)
    lon_skip_factor = max(round(lon_grid_desired/lon_grid_mean),1)
    print('/n--- Skip factors ---')
    print('Lat:',lat_skip_factor)
    print('Lon:',lon_skip_factor)
    #
    lat_start_value = int(np.floor(np.remainder(len(lat)-1,lat_skip_factor)/2))
    lon_start_value = int(np.floor(np.remainder(len(lon)-1,lon_skip_factor)/2))
    if dataset_txt == 'era5': lon_start_value += 1
    j_for_ts = np.arange(lat_start_value,len(lat),lat_skip_factor)
    i_for_ts = np.arange(lon_start_value,len(lon),lon_skip_factor)
    print('/n=== GRID WITH SKIPPED POINTS ===')
    print('Lat:',np.abs(np.mean(lat[j_for_ts][1:]-lat[j_for_ts][:-1])))
    print('Lon:',np.abs(np.mean(lon[i_for_ts][1:]-lon[i_for_ts][:-1])))
    print('Total points:',len(lat[j_for_ts])*len(lon[i_for_ts]))
    #
else:
    #
    j_for_ts = np.arange(len(lat))
    i_for_ts = np.arange(len(lon))


#%% Print grid lats and lons

if dataset_txt != 'holocenehydroclimate':
    #
    # Make a version of the longitude that goes from -180 to 180
    lon_neg = copy.deepcopy(lon)
    lon_neg[lon_neg > 180] = lon_neg[lon_neg > 180] - 360
    #
    # Make versions for printing
    lat_string = ','.join([str('{:.1f}'.format(value)) for value in lat[j_for_ts]])
    lon_string = ','.join([str('{:.1f}'.format(value)) for value in lon_neg[i_for_ts]])
    #
    lat_string = "          var lat_all = ["+lat_string+"];"
    lon_string = "          var lon_all = ["+lon_string+"];"
    #
    # Save the latitudes and longitudes to a file
    with open('latlon_'+filename_txt+'.txt','w') as f:
        f.write(lat_string+'\n')
        f.write(lon_string)


#%% FIGURES

# Get the colors from the colorbar. This is important for non-linear colorbars
colors_from_cmap = matplotlib.colormaps[cmap]
n_colors = len(levels)+1
colors_selected = colors_from_cmap(np.linspace(0,1,n_colors))

# Make 2D lon bound variables
lon_bounds_2d,lat_bounds_2d = np.meshgrid(lon_bounds,lat_bounds)


#%% If the reconstruction uses the IPCC AR6 regions, do some processing
if map_type == 'regions_ipcc_ar6':
    #
    # Get all WGI regions
    ar6_all = regionmask.defined_regions.ar6.all
    ar6_abbreviations = ar6_all.abbrevs
    regions_all = lat
    #
    # Create some variables to use later
    n_methods = var_ens.shape[0]
    n_ens     = var_ens.shape[1]
    n_time    = var_ens.shape[2]
    n_regions = len(ar6_abbreviations)
    var_regional_ens = np.zeros((n_methods,n_ens,n_time,n_regions)); var_regional_ens[:] = np.nan
    var_regional     = np.zeros((n_methods,n_time,n_regions));       var_regional[:]     = np.nan
    for i,region_txt in enumerate(ar6_abbreviations):
        ind_region = np.where(region_txt==regions_all)[0][0]
        var_regional_ens[:,:,:,i] = var_ens[:,:,:,ind_region,0]
        var_regional[:,:,i]       = var_mean[:,:,ind_region,0]


#%% MAPS
# Make a map of values
print('Step 1: Making maps: '+str(len(time_var)))
i=0;time=time_var[i]
for i,time in enumerate(time_var):
    #
    # Make a text box to show on the website
    plt.figure(figsize=(4,2))
    ax1 = plt.subplot2grid((1,1),(0,0))
    ax1.axis('off')
    ax1.fill_between(time_var,np.percentile(globalmean_all_ens,2.5,axis=0),np.percentile(globalmean_all_ens,97.5,axis=0),color='gray',alpha=0.25)
    ax1.plot(time_var,globalmean_all)
    ax1.axvline(x=time,color='gray',alpha=1,linestyle='--',linewidth=1)
    ax1.axhline(y=0,color='gray',alpha=0.5,linestyle='--',linewidth=1)
    ax1.set_xlim(x_range[0],x_range[1]+(x_range[1]-x_range[0])/100)
    if dataset_txt == 'holocenehydroclimate': chosen_txt = 'Mean'
    else: chosen_txt = 'Global'
    ax1.set_title(chosen_txt+': '+str('{:.2f}'.format(globalmean_all[i]))+' '+unit_txt,fontsize=18)
    if save_instead_of_plot:
        plt.savefig(output_dir+'info_'+filename_txt+'_'+str(int(np.ceil(time))).zfill(5)+'.png',dpi=50,format='png',bbox_inches='tight')
        plt.close()
    else:
        plt.show()
    #
    #
    #%%
    # Make the primary map to show
    plt.figure(figsize=(10,10))
    ax1 = plt.subplot2grid((1,1),(0,0),projection=ccrs.Mercator(central_longitude=0,min_latitude=-85,max_latitude=85))
    #ax1.set_extent([-180,180,-85,85],crs=ccrs.PlateCarree())
    #ax1.set_extent([-179.9,179.9,-85,85],crs=ccrs.PlateCarree())
    if map_type == 'contourf':
        var_cyclic,lon_cyclic = cutil.add_cyclic_point(var_mean_allmethods[i,:,:],coord=lon)
        map1 = ax1.contourf(lon_cyclic,lat,var_cyclic,colors=colors_selected,levels=levels,extend='both',transform=ccrs.PlateCarree())
        colorbar = plt.colorbar(map1,ticks=levels,orientation='horizontal',ax=ax1,fraction=0.01,pad=-0.07)
    elif map_type == 'pcolormesh':
        map1 = ax1.pcolormesh(lon_bounds_2d,lat_bounds_2d,var_mean_allmethods[i,:,:],cmap=cmap,vmin=-2,vmax=2,transform=ccrs.PlateCarree())
        colorbar = plt.colorbar(map1,orientation='horizontal',ax=ax1,fraction=0.01,pad=-0.07)
    elif map_type == 'regions_ipcc_ar6':
        norm = matplotlib.colors.Normalize(vmin=-2,vmax=2,clip=True)
        for j in range(len(ar6_all)):
            region_txt = ar6_abbreviations[j]
            if region_txt in regions_all:
                ind_region = np.where(region_txt==regions_all)[0][0]
                value_for_region = var_mean_allmethods[i,ind_region,0]
            else: value_for_region = np.nan
            if np.isnan(value_for_region): facecolor = 'lightgray'
            else: facecolor = colors_from_cmap(norm(value_for_region))
            region = ShapelyFeature([ar6_all[j]._polygon],facecolor=facecolor,crs=ccrs.PlateCarree())
            ax1.add_feature(region)
        colorbar = plt.colorbar(matplotlib.cm.ScalarMappable(norm=norm,cmap=cmap),orientation='horizontal',ax=ax1,fraction=0.01,pad=-0.07)
    if np.isnan(var_mean_allmethods[i,:,:]).all(): ax1.text(0.5,0.5,'Please select another year.',fontsize=12,ha='center',va='center',transform=ax1.transAxes)
    ax1.coastlines()
    gl = ax1.gridlines(color='gray',linestyle=':',draw_labels=False)
    gl.ylocator = mticker.FixedLocator(np.arange(-90,91,30))
    gl.xlocator = mticker.FixedLocator(np.arange(-180,181,60))
    colorbar.set_label(colorbar_txt+', rel. '+ref_period_txt,fontsize=6)
    colorbar.ax.tick_params(labelsize=3)
    plt.text(0,-82,dataset_name+', v.'+version_txt.replace('_','.')+', '+str(time_var[i])+' '+time_unit_txt,fontsize=7,horizontalalignment='center',transform=ccrs.PlateCarree())
    #
    if save_instead_of_plot:
        plt.savefig(output_dir+'map_'+filename_txt+'_'+str(int(np.ceil(time))).zfill(5)+'.png',dpi=150,format='png',bbox_inches='tight',pad_inches=0.0)
        plt.close()
    else:
        plt.show()


#%% COLORBAR
i=0;time=time_var[i]
plt.figure(figsize=(10,10))
ax1 = plt.subplot2grid((1,1),(0,0),projection=ccrs.Mercator(central_longitude=0,min_latitude=-85,max_latitude=85))
#ax1.set_extent([-180,180,-85,85],crs=ccrs.PlateCarree())
if map_type == 'contourf':
    var_cyclic,lon_cyclic = cutil.add_cyclic_point(var_mean_allmethods[i,:,:],coord=lon)
    map1 = ax1.contourf(lon_cyclic,lat,var_cyclic,colors=colors_selected,levels=levels,extend='both',transform=ccrs.PlateCarree())
    colorbar = plt.colorbar(map1,ticks=levels,orientation='horizontal',ax=ax1,fraction=0.08,pad=0.02)
elif map_type == 'pcolormesh':
    map1 = ax1.pcolormesh(lon_bounds_2d,lat_bounds_2d,var_mean_allmethods[i,:,:],cmap=cmap,vmin=-2,vmax=2,transform=ccrs.PlateCarree())
    colorbar = plt.colorbar(map1,orientation='horizontal',ax=ax1,fraction=0.08,pad=0.02)
elif map_type == 'regions_ipcc_ar6':
    norm = matplotlib.colors.Normalize(vmin=-2,vmax=2,clip=True)
    colorbar = plt.colorbar(matplotlib.cm.ScalarMappable(norm=norm,cmap=cmap),orientation='horizontal',ax=ax1,fraction=0.08,pad=0.02)
plt.gca().set_visible(False)
colorbar.set_label(colorbar_txt+', rel. '+ref_period_txt,fontsize=16)
colorbar.ax.tick_params(labelsize=12)

if save_instead_of_plot:
    plt.savefig(output_dir+'colorbar_'+filename_txt+'.png',dpi=150,format='png',bbox_inches='tight')
    plt.close()
else:
    plt.show()


#%% TIME SERIES

color_list = ['black','royalblue','salmon','olive','orange','darkseagreen']

# Make a timeseries at every location
if make_gridded_ts:
    j,i = 0,0
    n_total = len(j_for_ts)*len(i_for_ts)
    print('Step 2: Making time series at points: '+str(n_total))
    for j in j_for_ts:
        for i in i_for_ts:
            #
            # Make an interactive time series with bokeh
            lat_txt = str('{:.1f}'.format(lat[j]))
            lon_txt = str('{:.1f}'.format(lon_neg[i]))
            p1 = figure(width=1200,
                        height=ts_height,
                        title=title_txt_bokeh+' for '+dataset_name+', v.'+version_txt.replace('_','.')+' near '+lat_txt+'\u00B0N, '+lon_txt+'\u00B0E',
                        tools='pan,box_zoom,hover,save,reset',
                        active_drag='box_zoom',active_inspect='hover',
                        x_range=Range1d(bounds=(min(x_range),max(x_range))))
            #
            p1.title.text_color = 'black'
            p1.xaxis.axis_label = time_name_txt+' ('+time_unit_txt+')'
            p1.yaxis.axis_label = title_txt_bokeh
            p1.x_range.start = x_range[0]
            p1.x_range.end   = x_range[1]
            p1.y_range.start = ts_yrange[0]
            p1.y_range.end   = ts_yrange[1]
            #
            for k,method_chosen in enumerate(method):
                if dataset_txt == 'era20c':
                    pass  # ERA20C does not have an ensemble
                elif dataset_txt == 'lmr': #TODO: Improve this later
                    p1.varea(time_var,np.mean(var_lower_2std[k,:,:,j,i].values,axis=0),np.mean(var_upper_2std[k,:,:,j,i].values,axis=0),color=color_list[k],alpha=0.1,legend_label=method_chosen)
                else:
                    p1.varea(time_var,np.percentile(var_ens[k,:,:,j,i].values,2.5,axis=0),np.percentile(var_ens[k,:,:,j,i].values,97.5,axis=0),color=color_list[k],alpha=0.1,legend_label=method_chosen)
                p1.line(time_var,var_mean.values[k,:,j,i],color=color_list[k],line_width=1,legend_label=method_chosen)
            line0 = Span(location=0,dimension='width',line_color='gray',line_width=1)
            p1.renderers.extend([line0])
            p1.background_fill_color           = 'white'
            p1.grid.grid_line_color            = '#e0e0e0'
            p1.axis.axis_label_text_font_style = 'normal'
            p1.axis.axis_label_text_font_size  = '16px'
            p1.title.text_font_size            = '16px'
            p1.title.align                     = 'center'
            p1.legend.location     = 'bottom_right'
            p1.legend.click_policy = 'hide'
            #
            hover = p1.select_one(HoverTool)
            hover.tooltips = [
                    (time_name_txt,'@x{int} '+time_unit_txt),
                    ('Temp','@y \u00B0C'),
                    ]
            #hover.mode='vline'
            #
            # Save as html
            outputfile_txt = 'ts_'+filename_txt+'_lat_'+lat_txt+'_lon_'+lon_txt+'.html'
            html = file_html(p1,CDN,outputfile_txt)
            output_file = open(output_dir+outputfile_txt,'w')
            output_file.write(html)
            output_file.close()


#%% Make regional time series plots, if requested
if make_regional_ts:
    #
    if dataset_txt != 'holocenehydroclimate':
        #
        ### Compute regional means
        # Get all WGI regions
        ar6_all = regionmask.defined_regions.ar6.all
        #
        # Make a mask for the different regions
        mask_3D = ar6_all.mask_3D(lon,lat)
        #
        # Calculate weights for every gridcell
        lon_2d,lat_2d = np.meshgrid(lon,lat)
        lat_weights = np.cos(np.deg2rad(lat_2d))
        #
        # Compute regional means
        var_regional     = var_mean.weighted(mask_3D * lat_weights).mean(dim=('lat','lon')).values
        var_regional_ens = var_ens.weighted(mask_3D * lat_weights).mean(dim=('lat','lon')).values
        if dataset_txt == 'lmr': #TODO: Improve this later
            var_regional_ens_stdev = var_ens_stdev.weighted(mask_3D * lat_weights).mean(dim=('lat','lon')).values
            var_upper_2std = var_regional_ens + (2*var_regional_ens_stdev)
            var_lower_2std = var_regional_ens - (2*var_regional_ens_stdev)
    #
    #
    ### Make regional plots
    n_regions = len(ar6_all.abbrevs)
    print('Step 3: Making time series for regions: '+str(n_regions))
    j = 0
    for j in range(n_regions):
        #
        region_number = str(j+1)
        region_abbrev = ar6_all.abbrevs[j]
        region_name   = ar6_all.names[j]
        ts_to_plot    = var_regional[:,:,j]
        #
        # Make an interactive time series with bokeh
        p1 = figure(width=1200,
                    height=ts_height,
                    title=title_txt_bokeh+' for '+dataset_name+', v.'+version_txt.replace('_','.')+' for region '+region_abbrev+' ('+region_name+')',
                    tools='pan,box_zoom,hover,save,reset',
                    active_drag='box_zoom',active_inspect='hover',
                    x_range=Range1d(bounds=(min(x_range),max(x_range))))
        #
        p1.title.text_color = 'green'
        p1.xaxis.axis_label = time_name_txt+' ('+time_unit_txt+')'
        p1.yaxis.axis_label = title_txt_bokeh
        p1.x_range.start = x_range[0]
        p1.x_range.end   = x_range[1]
        p1.y_range.start = ts_yrange[0]
        p1.y_range.end   = ts_yrange[1]
        #
        for k,method_chosen in enumerate(method):
            if dataset_txt == 'era20c':
                pass  # ERA20C does not have an ensemble
            elif dataset_txt == 'lmr': #TODO: Improve this later
                p1.varea(time_var,np.mean(var_lower_2std[k,:,:,j].values,axis=0),np.mean(var_upper_2std[k,:,:,j].values,axis=0),color=color_list[k],alpha=0.1,legend_label=method_chosen)
            else:
                p1.varea(time_var,np.percentile(var_regional_ens[k,:,:,j],2.5,axis=0),np.percentile(var_regional_ens[k,:,:,j],97.5,axis=0),color=color_list[k],alpha=0.1,legend_label=method_chosen)
            p1.line(time_var,ts_to_plot[k,:],color=color_list[k],line_width=1,legend_label=method_chosen)
        line0 = Span(location=0,dimension='width',line_color='gray',line_width=1)
        p1.renderers.extend([line0])
        p1.background_fill_color           = 'white'
        p1.grid.grid_line_color            = '#e0e0e0'
        p1.axis.axis_label_text_font_style = 'normal'
        p1.axis.axis_label_text_font_size  = '16px'
        p1.title.text_font_size            = '16px'
        p1.title.align                     = 'center'
        p1.legend.location     = 'bottom_right'
        p1.legend.click_policy = 'hide'
        #
        hover = p1.select_one(HoverTool)
        hover.tooltips = [
                (time_name_txt,'@x{int} '+time_unit_txt),
                ('Temp','@y \u00B0C'),
                ]
        #hover.mode='vline'
        #
        # Save as html
        outputfile_txt = 'ts_'+filename_txt+'_region_'+region_abbrev+'.html'
        html = file_html(p1,CDN,outputfile_txt)
        output_file = open(output_dir+outputfile_txt,'w')
        output_file.write(html)
        output_file.close()


#%%
print('=== FINISHED ===')

