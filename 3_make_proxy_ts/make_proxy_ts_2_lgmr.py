#==============================================================================
# Make a set of dynamic html files of proxies.
#    author: Michael P. Erb
#    date  : 5/5/2023
#==============================================================================

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from bokeh.plotting import figure
from bokeh.resources import CDN
from bokeh.embed import file_html
from bokeh.models import HoverTool,Div
from bokeh.models.widgets import Panel,Tabs

plt.style.use('ggplot')

# Choose dataset
version_txt = '1_0_0'


#%% LOAD DATA

data_dir   = '/projects/pd_lab/data/paleoclimate_reconstructions/Holocene_reconstructions/Osman_etal_2021/'
output_dir = '/projects/pd_lab/mpe32/figures_presto/'

# Load the proxy data
proxy_database = Dataset(data_dir+'proxyDatabase.nc')


#%% Figure out what's in the proxy database
proxy_names = list(proxy_database.groups.keys())
n_proxies = len(proxy_names)
attributes_all      = []
variables_all       = []
data_attributes_all = []
for i in range(n_proxies):
    attributes = list(proxy_database[proxy_names[i]].ncattrs())
    attributes_all = attributes_all + attributes
    variables = proxy_database[proxy_names[i]].groups['data'].variables.keys()
    variables_all = variables_all + list(variables)
    for var in variables:
        if (var[:3] == 'age') or (var == 'depth'): continue
        data_attributes = list(proxy_database[proxy_names[i]].groups['data'].variables[var].ncattrs())
        data_attributes_all = data_attributes_all + data_attributes

# Print some values
print('N_proxies:',n_proxies)
print(' === PROXY ATTRIBUTES ===')
print(np.unique(attributes_all,return_counts=True))
print(' === PROXY VARIABLES ===')
print(np.unique(variables_all,return_counts=True))
print(' === DATA ATTRIBUTES ===')
print(np.unique(data_attributes_all,return_counts=True))


#%%
lats_all = np.zeros(n_proxies); lats_all[:] = np.nan
lons_all = np.zeros(n_proxies); lons_all[:] = np.nan
latlon_all = []
counters = {}
for i in range(n_proxies):
    name = proxy_names[i]
    lats_all[i] = proxy_database[name].latitude
    lons_all[i] = proxy_database[name].longitude
    latlon = str(lats_all[i])+'_'+str(lons_all[i])
    if latlon not in latlon_all:
        counters[latlon] = 1
    else:
        # If the lat/lon already appears in the list, shift the new point east by 0.1 degree
        lons_all[i] = lons_all[i] + counters[latlon]*0.1
        counters[latlon] += 1
    #
    latlon_all.append(latlon)


#%% TIME SERIES
# Make a timeseries at every location
i=0
for i in range(n_proxies):
    #
    # Get metadata
    name = proxy_names[i]
    proxy_name    = proxy_database[name].site_name
    proxy_lat     = proxy_database[name].latitude
    proxy_lon     = proxy_database[name].longitude
    proxy_elev    = proxy_database[name].elevation
    proxy_comment = proxy_database[name].comment
    proxy_year    = proxy_database[name].collection_year
    proxy_ref     = proxy_database[name].reference
    proxy_ref_short = proxy_ref[:30]+'...'
    #
    # Get proxy data
    data = proxy_database[name].groups['data'].variables
    proxy_ages = data['age_median'][:]
    age_units  = data['age_median'].units
    age_name   = data['age_median'].long_name
    data_variables = data.keys()
    var = 'd18o_pachyderma'
    tab = {}
    tab_list = []
    for var in data_variables:
        if (var[:3] == 'age') or (var == 'depth'): continue
        #
        proxy_data = data[var][:]
        data_units = data[var].units
        try:    data_species = data[var].species
        except: data_species = ''
        try:    data_cleaning = data[var].mgca_cleaning_protocol
        except: data_cleaning = ''
        #
        title_txt = 'Proxy data: '+name+'  |  '+var+'  |  '+data_species+'  |  '+proxy_ref_short
        #
        # Make an interactive time series with bokeh
        p1 = figure(plot_width=1200,
                    plot_height=300,
                    title=title_txt,
                    tools='pan,box_zoom,hover,save,reset',
                    active_drag='box_zoom',active_inspect='hover')
        #
        p1.title.text_color = 'purple'
        p1.xaxis.axis_label = 'Age ('+age_units+')'
        p1.yaxis.axis_label = data_units
        p1.x_range.start = 24000
        p1.x_range.end   = 0
        #p1.y_range.start = ts_yrange[0]
        #p1.y_range.end   = ts_yrange[1]
        #
        p1.scatter(proxy_ages,proxy_data,color='purple',size=10,marker='dot')
        p1.line(proxy_ages,proxy_data,color='purple',line_width=1)
        p1.background_fill_color           = 'white'
        p1.grid.grid_line_color            = '#e0e0e0'
        p1.axis.axis_label_text_font_style = 'normal'
        p1.axis.axis_label_text_font_size  = '16px'
        p1.title.text_font_size            = '16px'
        p1.title.align                     = 'center'
        #
        hover = p1.select_one(HoverTool)
        hover.tooltips = [
                #('Age','@x{int} '+age_units),
                ('Age','@x '+age_units),
                ('Data','@y '+data_units),
                ]
        hover.mode='vline'
        #
        tab[var] = Panel(child=p1,title=var)
        tab_list.append(tab[var])
    #
    # Put the tabs together
    full_plot = Tabs(tabs=tab_list)
    #
    # Save as html
    outputfile_txt = 'proxy_lgmr_v'+version_txt+'_'+str(i).zfill(5)+'.html'
    html = file_html(full_plot,CDN,outputfile_txt)
    output_file = open(output_dir+outputfile_txt,'w')
    output_file.write(html)
    output_file.close()


#%% PRINT LAT/LONS FOR WEBSITE

# Make versions for printing
lat_string = ','.join([str('{:.2f}'.format(value)) for value in lats_all])
lon_string = ','.join([str('{:.2f}'.format(value)) for value in lons_all])

lat_string = "              lat_proxies = ["+lat_string+"];"
lon_string = "              lon_proxies = ["+lon_string+"];"

# Save the latitudes and longitudes to a file
with open('latlon_lgmr_proxies.txt','w') as f:
    f.write(lat_string+'\n')
    f.write(lon_string)

print('=== FINISHED ===')

