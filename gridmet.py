# Libraries
import pandas as pd
import pickle
import time
import hvplot.pandas
import hvplot.xarray
import geopandas as gpd
import grd2shp_xagg
import xagg as xa
import xarray as xr

# Variables
gdf_nhru02_path = './data/nhru_02/nhru_02.shp'
start_date = '1979-01-01'
end_date = '1989-01-01'
lon_min = -79
lon_max = -67
lat_min = 36
lat_max = 42
output_path = './data/'
## official list of variables needed for drb-inland-salinity model
data_vars_shrt_all = ['tmmx', 'tmmn', 'pr', 'srad', 'vs','rmax','rmin','sph']
data_vars_long_all = ['daily_maximum_temperature', 'daily_minimum_temperature',
         'precipitation_amount', 'daily_mean_shortwave_radiation_at_surface',
         'daily_mean_wind_speed', 'daily_maximum_relative_humidity',
         'daily_minimum_relative_humidity','daily_mean_specific_humidity']

## Shorter list of data variables for testing:
data_vars_grd2shp_shrt = ['tmmn', 'vs']

## Full function
def g2shp_regridding(polygon_file_path,
                     output_data_folder,
                     var_short,
                     start_date, end_date,
                     lon_min, lon_max, lat_min, lat_max,
                     g2s_file_prefix = 'tmp'):
    """
    :param str polygon_file_path:
    :param str output_data_folder:
    :param str list var_short: data variable short name or list of data variable short names. Must be one of the following: ['tmmx', 'tmmn', 'pr', 'srad', 'vs', 'rmax', 'rmin', 'sph']
    :param str start_date:
    :param str end_date:
    :param str lon_min:
    :param str lon_max:
    :param str lat_min:
    :param str lat_max:
    :return: g2s array
    """
    ## Create dictionary of short variable names with official long variable name
    data_vars_shrt_all = ['tmmx', 'tmmn', 'pr', 'srad', 'vs', 'rmax', 'rmin', 'sph']
    data_vars_long_all = ['daily_maximum_temperature', 'daily_minimum_temperature',
                          'precipitation_amount', 'daily_mean_shortwave_radiation_at_surface',
                          'daily_mean_wind_speed', 'daily_maximum_relative_humidity',
                          'daily_minimum_relative_humidity', 'daily_mean_specific_humidity']

    all_data_vars_dict = dict(zip(data_vars_shrt_all, data_vars_long_all))
    subset = {key: all_data_vars_dict[key] for key in var_short}
    print(subset)
    del all_data_vars_dict

    ## geospatial vector file
    gdf = gpd.read_file(polygon_file_path)

    ## check if only 1 var (short and long) is provided, change to list
    if type(var_short) == str:
        var_short = [var_short]
        print(var_short)
    ## check that var_short is a list
    elif type(var_short) != list:
        var_short = [var_short]
        print(var_short)

    ## empty list for the xarray.datasets
    data_xr_ds_list = []

    ## Creating ds
    for i in var_short:
        # Get url:
        print(i)
        url = 'http://thredds.northwestknowledge.net:8080/thredds/dodsC/agg_met_{}_1979_CurrentYear_CONUS.nc'.format(i)
        print(url)

        # call data from url
        start = time.perf_counter()
        df = xr.open_dataset(url + '#fillmismatch')

        # Subset to timeframe and bbox of interest:
        df_subset = df.sel(day=slice(start_date, end_date),
                           lon=slice(lon_min, lon_max),
                           lat=slice(lat_max, lat_min))
        end = time.perf_counter()
        print('finish in {} seconds'.format(round(start - end), 2))

        # Append to list of xr.datasets
        data_xr_ds_list.append(df_subset)

    ## Creating dictionary for next process
    data_dict = dict(zip(var_short, data_xr_ds_list))
    ## taking first dict item for creation of weight's file.
    wghtmap_var = list(data_dict.keys())[0]

    ## Produce weightmap
    start = time.perf_counter()
    weightmap = xa.pixel_overlaps(data_dict[wghtmap_var], gdf)
    end = time.perf_counter()
    print(f'finished agg in {round(end - start, 2)} second(s)')

    with open(output_data_folder + 'grd2shp_weights_{}.pickle'.format(wghtmap_var), 'wb') as file:
        pickle.dump(weightmap, file)

    ## initialize Grd2ShpXagg()
    g2s = grd2shp_xagg.Grd2ShpXagg()
    g2s.initialize(
        grd=list(data_dict.values()),
        shp=gdf,
        wght_file= output_data_folder + 'grd2shp_weights_{}.pickle'.format(wghtmap_var),
        time_var='day',
        lat_var='lat',
        lon_var='lon',
        var= list(subset.values()),
        var_output=list(data_dict.keys()),
        ctype=0,
    )

    # Regridding
    start = time.perf_counter()
    g2s.run_weights()
    end = time.perf_counter()
    print(f'finished agg in {round(end - start, 2)} second(s)')

    print(type(g2s))
    g2s.write_gm_file(opath=output_data_folder, prefix= g2s_file_prefix)

    return g2s


g2shp_regridding(polygon_file_path = gdf_nhru02_path,
                var_short = data_vars_grd2shp_shrt,
                 output_data_folder=output_path,
                 start_date=start_date, end_date=end_date,
                 lat_max=lat_max, lat_min=lat_min,
                 lon_max=lon_max, lon_min=lon_min,
                 g2s_file_prefix = 'test3_')

# xr_mapped = xr.open_dataset('t_srad_prclimate_2022_01_19.nc', decode_times=False)
#
# gdf_nhru02 = gpd.read_file(polygon_file_path)
# gdf_nhru02["pr"] = xr_mapped["pr"][:,0]


######## OLD CODE ################################

# ## API urls
# # tmax_url = 'http://thredds.northwestknowledge.net:8080/thredds/dodsC/agg_met_tmmx_1979_CurrentYear_CONUS.nc'
# # tmin_url = 'http://thredds.northwestknowledge.net:8080/thredds/dodsC/agg_met_tmmn_1979_CurrentYear_CONUS.nc'
# # prcp_url = 'http://thredds.northwestknowledge.net:8080/thredds/dodsC/agg_met_pr_1979_CurrentYear_CONUS.nc'
# # swave_rad_url = 'http://thredds.northwestknowledge.net:8080/thredds/dodsC/agg_met_srad_1979_CurrentYear_CONUS.nc'
# # wind_vel_url = 'http://thredds.northwestknowledge.net:8080/thredds/dodsC/agg_met_vs_1979_CurrentYear_CONUS.nc'
# # rmax_url = 'http://thredds.northwestknowledge.net:8080/thredds/dodsC/agg_met_rmax_1979_CurrentYear_CONUS.nc'
# # rmin_url = 'http://thredds.northwestknowledge.net:8080/thredds/dodsC/agg_met_rmin_1979_CurrentYear_CONUS.nc'
# # sph_url = 'http://thredds.northwestknowledge.net:8080/thredds/dodsC/agg_met_sph_1979_CurrentYear_CONUS.nc'
#
# if type(data_vars)== str:
#     data_vars = [data_vars]
#     print(data_vars)
#
# data_xr_ds_list = []
#
# ## Creating ds
# for i in data_vars:
#     print(i)
#     url = 'http://thredds.northwestknowledge.net:8080/thredds/dodsC/agg_met_{}_1979_CurrentYear_CONUS.nc'.format(i)
#     print(url)
#     start = time.perf_counter()
#     df = xr.open_dataset(url +'#fillmismatch')
#     df_subset = df.sel(day=slice(start_date, end_date),
#                          lon=slice(lon_min, lon_max),
#                          lat=slice(lat_max, lat_min))
#     end = time.perf_counter()
#     data_xr_ds_list.append(df_subset)
#
#     print('finish in {} seconds'.format(round(start - end), 2))
#
# # Creating dictionary for next process
# data_dict = dict(zip(data_vars, data_xr_ds_list))
#
# print(data_dict)
#
# data_xr_weightmap_list = []
#
# ## Did this for all the files bit not necessary
# for j in data_dict:
#     # Produce weightmap
#     # uncomment below if first time through notebook to generate weights
#     #start = time.perf_counter()
#     print(data_dict[j])
#     print(j)
#     weightmap = xa.pixel_overlaps(data_dict[j], gdf_nhru02)
#     end = time.perf_counter()
#     print(f'finished agg in {round(end-start, 2)} second(s)')
#
#     name = str(j)
#     print(name)
#     print('{}'.format(name))
#
#     # save weights for future use
#     with open('./data/nhru_02_weights_{}.pickle'.format(name), 'wb') as file:
#         pickle.dump(weightmap, file)
#
# # initialize Grd2ShpXagg()
# g2s = grd2shp_xagg.Grd2ShpXagg()
# g2s.initialize(
#     grd= list(data_dict.values()),
#     shp=gdf_nhru02,
#     wght_file='./data/nhru_02_weights_srad.pickle',
#     time_var='day',
#     lat_var='lat',
#     lon_var='lon',
#     var=['daily_mean_shortwave_radiation_at_surface', 'precipitation_amount'],
#
#     var_output= list(data_dict.keys()),
#     ctype=0,
# )
#
# ## RUNS SUCCESSFULLY TILL HER ##
# # Regridding
# start = time.perf_counter()
# g2s.run_weights()
# end = time.perf_counter()
# print(f'finished agg in {round(end-start, 2)} second(s)')
#
# print(type(g2s))
# g2s.write_gm_file(opath=".", prefix="t_sradpr_")
#
# xr_mapped = xr.open_dataset('t_srad_prclimate_2022_01_19.nc', decode_times=False)
#
# gdf_nhru02["tmax"] = xr_mapped["tmax"][:,0]
# ds_tmax_1["daily_maximum_temperature"] = ds_tmax_1["daily_maximum_temperature"] - 459.67
#
