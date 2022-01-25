# Libraries
import pickle
import time
import geopandas as gpd
import grd2shp_xagg
import xagg as xa
import xarray as xr
import os

# Variables
gdf_nhru02_path = './data/nhru_02/nhru_02.shp'
gdf = gpd.read_file(gdf_nhru02_path)
start_date = '1979-01-01'
end_date = '1980-01-01'
lon_min, lat_min, lon_max, lat_max = gdf.total_bounds
output_path = './data/'

## official list of variables needed for drb-inland-salinity model
data_vars_shrt_all = ['tmmx', 'tmmn', 'pr', 'srad', 'vs','rmax','rmin','sph']

def get_gridmet_datasets(variables, start_date, end_date, polygon_bbox = None, lon_min = None, lat_min = None, lon_max = None, lat_max = None):

    ## check bounds for data slicing
    if polygon_bbox is not None:
        lon_min, lat_min, lon_max,lat_max = polygon_bbox
        print(lon_min, lat_min, lon_max,lat_max)
    else:
        print(lon_min, lat_min, lon_max,lat_max)

    ## check if only 1 var (short) is provided, change to list
    if not isinstance(variables, list):
        variables = [variables]
        print(variables)

    ## empty list for the xarray.datasets
    data_xr_ds_list = []

    for var in variables:
        ## Pulling data
        # source: http://thredds.northwestknowledge.net:8080/thredds/reacch_climate_MET_aggregated_catalog.html
        url = f'http://thredds.northwestknowledge.net:8080/thredds/dodsC/agg_met_{var}_1979_CurrentYear_CONUS.nc'
        # call data from url
        start = time.perf_counter()
        ds = xr.open_dataset(url + '#fillmismatch')

        ## Subset to timeframe and bbox of interest:
        ds_subset = ds.sel(day = slice(start_date, end_date),
                           lon = slice(lon_min, lon_max),
                           lat = slice(lat_max, lat_min))
        end = time.perf_counter()
        print(f'finish in {round(start - end, 2)} seconds')

        # Append to list of xr.datasets
        data_xr_ds_list.append(ds_subset)

    ## Creating output dictionary
    xarray_dict = dict(zip(variables, data_xr_ds_list))

    return xarray_dict

xarray_dict = get_gridmet_datasets(variables = data_vars_shrt_all,
                                   start_date = start_date, end_date = end_date,
                                   polygon_bbox = gdf.total_bounds)
#                     lon_min = -80, lat_min = 36, lon_max = -71, lat_max = 45

def create_weightmap(polygon, xarray_dict, output_data_folder, weightmap_var = None):

    if weightmap_var is None:
    ## taking first dict item for creation of weight's file.
        weightmap_file = os.path.join(output_data_folder, f'grd2shp_weights.pickle')
        weightmap_var = list(xarray_dict.keys())[0]
        print(weightmap_file, weightmap_var)

    else:
        weightmap_file = os.path.join(output_data_folder, f'grd2shp_weights_{weightmap_var}.pickle')
        print(weightmap_file)

    ## Produce weightmap
    start = time.perf_counter()
#    ds = xarray_dict[weightmap_var].load()

    weightmap = xa.pixel_overlaps(xarray_dict[weightmap_var], polygon)
    end = time.perf_counter()
    print('finished agg in {} second(s)'.format(round(end - start, 2)))

    with open(weightmap_file, 'wb') as file:
        pickle.dump(weightmap, file)

    return weightmap_file

create_weightmap(polygon = gdf, xarray_dict = xarray_dict, output_data_folder = output_path, weightmap_var = 'tmmx')

def g2shp_regridding(xarray_dict, weightmap_file, output_data_folder, g2s_file_prefix):

    ## params for g2s.initialise()
    # grab long name of data vars from dict and place in list
    vars_long_list = []
    for key in xarray_dict:
        vars_long_list = vars_long_list + list(xarray_dict[key])
    # grab short name of data vars from dict keys
    vars_short_list = list(xarray_dict.keys())
    vars_grd_list = list(xarray_dict.values())

    ## initialize Grd2ShpXagg()
    g2s = grd2shp_xagg.Grd2ShpXagg()

    g2s.initialize(
        grd = vars_grd_list,
        shp = gdf,
        wght_file = weightmap_file,
        time_var = 'day',
        lat_var = 'lat',
        lon_var = 'lon',
        var = vars_long_list,
        var_output = vars_short_list,
        ctype=0,
    )

    # Regridding
    start = time.perf_counter()
    g2s.run_weights()
    end = time.perf_counter()
    print('finished agg in {} second(s)'.format(round(end - start, 2)))

    print(type(g2s))

    g2s.write_gm_file(opath=output_data_folder, prefix=g2s_file_prefix)

    return g2s

g2shp_regridding(xarray_dict = xarray_dict, weightmap_file= './data/grd2shp_weights_tmmx.pickle', g2s_file_prefix='gsgrb', output_data_folder=output_path)