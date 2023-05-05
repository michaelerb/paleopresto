#==============================================================================
# Make a set of dynamic html files of proxies.
#    author: Michael P. Erb
#    date  : 5/5/2023
#==============================================================================

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from bokeh.plotting import figure
from bokeh.resources import CDN
from bokeh.embed import file_html
from bokeh.models import HoverTool,Div
#from bokeh.layouts import column
from bokeh.models.widgets import Panel,Tabs

plt.style.use('ggplot')

# Choose dataset
version_txt = '1_0_0'


#%% LOAD DATA

data_dir   = '/projects/pd_lab/data/paleoclimate_reconstructions/LMR_final/input_data/proxies/'
output_dir = '/projects/pd_lab/mpe32/figures_presto/'

# Load the proxy data
proxy_data_all     = pd.read_pickle(data_dir+'LMRdb_v1.0.0_Proxies.df.pckl')
proxy_metadata_all = pd.read_pickle(data_dir+'LMRdb_v1.0.0_Metadata.df.pckl')


#%% TIME SERIES

proxy_names = list(proxy_data_all.keys())
n_proxies = len(proxy_names)

# Set up arrays for lat/lons
lats_all = np.zeros(n_proxies); lats_all[:] = np.nan
lons_all = np.zeros(n_proxies); lons_all[:] = np.nan
latlon_all = []
counters = {}
#for key in proxy_metadata.keys(): print(key,proxy_metadata[key].values)

# Put all proxy records with the same dataset name together and make a figure at every location
i=0
for i in range(n_proxies):
    #
    # Get proxy data and metadata
    proxy_name = proxy_names[i]
    proxy_data = proxy_data_all[proxy_name]
    proxy_metadata = proxy_metadata_all.loc[proxy_metadata_all['Proxy ID'] == proxy_name]
    #
    # Get lat and lons
    lats_all[i] = proxy_metadata['Lat (N)']
    lons_all[i] = proxy_metadata['Lon (E)']
    latlon = str(lats_all[i])+'_'+str(lons_all[i])
    if latlon not in latlon_all:
        counters[latlon] = 1
    else:
        # If the lat/lon already appears in the list, shift the new point east by 0.1 degree
        lons_all[i] = lons_all[i] + counters[latlon]*0.1
        counters[latlon] += 1
    #
    latlon_all.append(latlon)
    #
    #print(j,indices_selected)
    tab = {}
    tab_list = []
    #
    # Get data and metadata
    proxy_years  = np.array(proxy_data.index).astype(float)
    proxy_values = np.array(proxy_data.values).astype(float)
    #
    # Get proxy metadata
    metadata = {}
    metadata['Dataset name']   = str(proxy_metadata['Site name'].values[0])
    metadata['Archive type']   = str(proxy_metadata['Archive type'].values[0])
    metadata['Proxy type']     = str(proxy_metadata['Proxy measurement'].values[0])
    metadata['TSid']           = str(proxy_metadata['Proxy ID'].values[0])
    metadata['paper_title']    = str(proxy_metadata['Study name'].values[0])
    metadata['seasonality']    = str(proxy_metadata['Seasonality'].values[0])
    data_units = str(proxy_metadata['Proxy measurement'].values[0])
    age_units  = 'C.E.'
    if len(metadata['TSid']) <= 40: tsid_str = metadata['TSid']
    else:                           tsid_str = metadata['TSid'][0:40]+'...'
    #
    title_txt = 'Proxy data: '+metadata['Archive type']+'  |  '+metadata['Proxy type']+'  |  '+metadata['seasonality']+'  |  '+metadata['Dataset name']+'  |  '+tsid_str
    #
    # Make an interactive time series with bokeh
    p1 = figure(plot_width=1200,
                plot_height=300,
                title=title_txt,
                tools='pan,box_zoom,hover,save,reset',
                active_drag='box_zoom',active_inspect='hover')
    #
    p1.title.text_color = 'purple'
    p1.xaxis.axis_label = 'Year ('+age_units+')'
    p1.yaxis.axis_label = data_units
    p1.x_range.start = 0
    p1.x_range.end   = 2000
    #p1.y_range.start = ts_yrange[0]
    #p1.y_range.end   = ts_yrange[1]
    #
    p1.scatter(proxy_years,proxy_values,color='purple',size=10,marker='dot')
    p1.line(proxy_years,proxy_values,color='purple',line_width=1)
    p1.background_fill_color           = 'white'
    p1.grid.grid_line_color            = '#e0e0e0'
    p1.axis.axis_label_text_font_style = 'normal'
    p1.axis.axis_label_text_font_size  = '16px'
    p1.title.text_font_size            = '16px'
    p1.title.align                     = 'center'
    #
    hover = p1.select_one(HoverTool)
    hover.tooltips = [
            ('Year','@x '+age_units),
            ('Data','@y '+data_units),
            ]
    hover.mode='vline'
    #
    season = metadata['seasonality']
    tab[season] = Panel(child=p1,title=metadata['Proxy type']+', '+season)
    tab_list.append(tab[season])
    #
    # Put the tabs together
    all_tabs = Tabs(tabs=tab_list)
    #
    # Save as html
    full_plot = all_tabs
    outputfile_txt = 'proxy_lmrdb_v'+version_txt+'_'+str(i).zfill(5)+'.html'
    #html = file_html(p1,CDN,outputfile_txt)
    html = file_html(full_plot,CDN,outputfile_txt)
    output_file = open(output_dir+outputfile_txt,'w')
    output_file.write(html)
    output_file.close()


#%% PRINT LAT/LONS FOR WEBSITE

# Make sure that the longitude values go from -180 to 180
lons_all[lons_all > 180] = lons_all[lons_all > 180] - 360

# Make versions for printing
lat_string = ','.join([str('{:.2f}'.format(value)) for value in lats_all])
lon_string = ','.join([str('{:.2f}'.format(value)) for value in lons_all])

lat_string = "              lat_proxies = ["+lat_string+"];"
lon_string = "              lon_proxies = ["+lon_string+"];"

# Save the latitudes and longitudes to a file
with open(output_dir+'latlon_lmrdb_proxies.txt','w') as f:
    f.write(lat_string+'\n')
    f.write(lon_string)

print('=== FINISHED ===')

