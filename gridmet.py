# Libraries
import pandas as pd
import numpy as np
import zipfile
import pygeos
import folium
import pickle
import time
import hvplot.pandas
import hvplot.xarray
import geopandas as gpd

import grd2shp_xagg
import xagg as xa
import xarray as xr

### Variables
start_date = '1960-01-01'
end_date = '1980-01-07'
lon_min = -79
lon_max = -67
lat_min = 36
lat_max = 42

tmax_url = 'http://thredds.northwestknowledge.net:8080/thredds/dodsC/agg_met_tmmx_1979_CurrentYear_CONUS.nc'
tmin_url = 'http://thredds.northwestknowledge.net:8080/thredds/dodsC/agg_met_tmmn_1979_CurrentYear_CONUS.nc'
prcp_url = 'http://thredds.northwestknowledge.net:8080/thredds/dodsC/agg_met_pr_1979_CurrentYear_CONUS.nc'

ds_prcp = xr.open_dataset(prcp_url+'#fillmismatch')
ds_tmin = xr.open_dataset(tmin_url+'#fillmismatch')
ds_tmax = xr.open_dataset(tmax_url+'#fillmismatch')

### subset gridMET data
ds_prcp_1 = ds_prcp.sel(day=slice(start_date, end_date),
                        lon=slice(lon_min, lon_max),
                        lat=slice(lat_max, lat_min))
ds_tmin_1 = ds_tmin.sel(day=slice(start_date, end_date),
                        lon=slice(lon_min, lon_max),
                        lat=slice(lat_max, lat_min))
ds_tmax_1 = ds_tmax.sel(day=slice(start_date, end_date),
                        lon=slice(lon_min, lon_max),
                        lat=slice(lat_max, lat_min))

ds_tmax_1.load()
ds_tmin_1.load()
ds_prcp_1.load()

gdf_nhru02 = gpd.read_file('./data/nhru_02/nhru_02.shp')
# gdf_gdb = gpd.read_file('./data/GeospatialFabricFeatures_02.gdb', layer="nhru")

### Produce weightmap
# uncomment below if first time through notebook to generate weights
start = time.perf_counter()
weightmap = xa.pixel_overlaps(ds_tmax_1, gdf_nhru02)
end = time.perf_counter()
print(f'finished agg in {round(end-start, 2)} second(s)')

# save weights for future use
with open('./data/nhru_02_weights.pickle', 'wb') as file:
    pickle.dump(weightmap, file)

### initialize Grd2ShpXagg()
g2s = grd2shp_xagg.Grd2ShpXagg()
g2s.initialize(
    grd=[ds_tmax_1, ds_tmin_1, ds_prcp_1],
    shp=gdf_nhru02,
    wght_file='./data/nhru_02_weights.pickle',
    time_var='day',
    lat_var='lat',
    lon_var='lon',
    var=['daily_maximum_temperature', 'daily_minimum_temperature', 'precipitation_amount'],
    var_output=['tmax', 'tmin', 'prcp'],
    ctype=0,
)

### Regridding

start = time.perf_counter()
g2s.run_weights()
end = time.perf_counter()
print(f'finished agg in {round(end-start, 2)} second(s)')

print(type(g2s))
g2s.write_gm_file(opath="./data", prefix="t_")

xr_mapped = xr.open_dataset('./data/t_climate_2022_01_13.nc', decode_times=False)

gdf_nhru02["tmax"] = xr_mapped["tmax"][:,0]
ds_tmax_1["daily_maximum_temperature"] = ds_tmax_1["daily_maximum_temperature"] * 9/5 - 459.67

### plotting