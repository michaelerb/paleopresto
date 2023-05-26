# PaleoPresto

This directory contains tutorials and other code for the PReSto project. When launched, the PReSto website will be located at [paleopresto.com](https://paleopresto.com/).

## Tutorials

The "tutorials" directory contains examples and instructions for analyzing paleoclimate data

## Paleoclimate data processing

The other directories contain scripts to create images and html code for the paleopresto.com website.

### Data

To download the paleoclimate data that these scripts use, see comments within each "format_data" script.

### Python environment

This code can be run using Python3. To create a new environment for this code, install anaconda and then run the following lines on the command line:

```
conda create -n python3_presto numpy xarray matplotlib cartopy bokeh netCDF4 pandas spyder geopandas conda-forge::regionmask
conda activate python3_presto
pip install LiPD
```
