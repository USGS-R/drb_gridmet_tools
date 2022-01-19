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

# Variables
gdf_nhru02 = gpd.read_file('./data/nhru_02/nhru_02.shp')
start_date = '1979-01-01'
end_date = '1989-01-01'
lon_min = -79
lon_max = -67
lat_min = 36
lat_max = 42

data_vars = ['tmax', 'tmin', 'prcp', 'swave_rad', 'wdspeed','rmax','rmin','spec_hum']

## API urls
tmax_url = 'http://thredds.northwestknowledge.net:8080/thredds/dodsC/agg_met_tmmx_1979_CurrentYear_CONUS.nc'
tmin_url = 'http://thredds.northwestknowledge.net:8080/thredds/dodsC/agg_met_tmmn_1979_CurrentYear_CONUS.nc'
prcp_url = 'http://thredds.northwestknowledge.net:8080/thredds/dodsC/agg_met_pr_1979_CurrentYear_CONUS.nc'
swave_rad_url = 'http://thredds.northwestknowledge.net:8080/thredds/dodsC/agg_met_srad_1979_CurrentYear_CONUS.nc'
wind_vel_url = 'http://thredds.northwestknowledge.net:8080/thredds/dodsC/agg_met_vs_1979_CurrentYear_CONUS.nc'
rmax_url = 'http://thredds.northwestknowledge.net:8080/thredds/dodsC/agg_met_rmax_1979_CurrentYear_CONUS.nc'
rmin_url = 'http://thredds.northwestknowledge.net:8080/thredds/dodsC/agg_met_rmin_1979_CurrentYear_CONUS.nc'
sph_url = 'http://thredds.northwestknowledge.net:8080/thredds/dodsC/agg_met_sph_1979_CurrentYear_CONUS.nc'

data_url_list = [tmax_url, tmin_url, prcp_url, swave_rad_url, wind_vel_url, rmax_url, rmin_url, sph_url]

## Creating ds
for i in data_url_list:
    print(i)
    df = xr.open_dataset(i+'#fillmismatch')
    df_subset = df.sel(day=slice(start_date, end_date),
                         lon=slice(lon_min, lon_max),
                         lat=slice(lat_max, lat_min))

    data_xr_ds_list.append(df_subset)

data_dict = dict(zip(data_vars, data_xr_ds_list))

print(data_dict)

data_xr_weightmap_list = []

## Did this for all the files bit not necessary
for j in data_dict:
    # Produce weightmap
    # uncomment below if first time through notebook to generate weights
    start = time.perf_counter()
    print(data_dict[j])
    print(j)
    weightmap = xa.pixel_overlaps(data_dict[j], gdf_nhru02)
    end = time.perf_counter()
    print(f'finished agg in {round(end-start, 2)} second(s)')

    # data_xr_weightmap_list.append(weightmap)

    name = str(j)
    print(name)
    print('{}'.format(name))

    with open('./data/nhru_02_weights_{}.pickle'.format(name), 'wb') as file:
        pickle.dump(weightmap, file)

# save weights for future use


# initialize Grd2ShpXagg()
g2s = grd2shp_xagg.Grd2ShpXagg()
g2s.initialize(
    grd= data_dict.values,
    shp=gdf_nhru02,
    wght_file='./data/nhru_02_weights_tmax.pickle',
    time_var='day',
    lat_var='lat',
    lon_var='lon',
    var=['daily_maximum_temperature', 'daily_minimum_temperature',
         'precipitation_amount', 'daily_mean_shortwave_radiation_at_surface',
         'daily_mean_wind_speed', 'daily_maximum_relative_humidity',
         'daily_minimum_relative_humidity','daily_mean_specific_humidity'],

    var_output= data_vars,
    ctype=0,
)

# Regridding
start = time.perf_counter()
g2s.run_weights()
end = time.perf_counter()
print(f'finished agg in {round(end-start, 2)} second(s)')

print(type(g2s))
g2s.write_gm_file(opath="./data", prefix="t_")

xr_mapped = xr.open_dataset('./data/t_climate_2022_01_13.nc', decode_times=False)

gdf_nhru02["tmax"] = xr_mapped["tmax"][:,0]
ds_tmax_1["daily_maximum_temperature"] = ds_tmax_1["daily_maximum_temperature"] - 459.67

# plotting