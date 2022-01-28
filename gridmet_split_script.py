# Libraries
import pickle
import time
import geopandas as gpd
import geopandas.geodataframe
import grd2shp_xagg
import xagg as xa
import xarray as xr
import os
import datetime

# Function definitions
def get_gridmet_datasets(variable, start_date, end_date, polygon_for_bbox = None, lon_min = None, lat_min = None, lon_max = None, lat_max = None):

    """
    :param str/list variable: data variable short name or list of data variables in short name. Must be one of the following: ['tmmx', 'tmmn', 'pr', 'srad', 'vs', 'rmax', 'rmin', 'sph']
    :param str start_date: Start date of data collection yyyy-mm-dd
    :param str end_date: End date of data collection yyyy-mm-dd
    :param gdf.total_bounds polygon_for_bbox: total bounds of a geodataframe. Ex: geodataframe.total_bounds. If not known, None
    :param str lon_min: bbox of aoi longitude min. None if polygon_for_bbox is given (i.e. polygon_bbox != None)
    :param str lat_min: bbox of aoi latitude min. None if polygon_for_bbox is given (i.e. polygon_bbox != None)
    :param str lon_max: bbox of aoi longitude max. None if polygon_for_bbox is given (i.e. polygon_bbox != None)
    :param str lat_max: bbox of aoi latitude max. None if polygon_for_bbox is given (i.e. polygon_bbox != None)
    :return: dictionary of xarray stored by variable name
    """

    ## check/define bounds for data slicing
    if polygon_for_bbox is not None:
        if isinstance(polygon_for_bbox, geopandas.geodataframe.GeoDataFrame):
            print('polygon is geodataframe')
            print(polygon_for_bbox.total_bounds)
            pass
        elif os.path.isfile(polygon_for_bbox):
            print('polygon is path to shapefile')
            polygon_for_bbox = gpd.read_file(polygon_for_bbox)
            print(polygon_for_bbox.total_bounds)
        else:
             raise TypeError('polygon_for_bbox should str path to shapefile or geodataframe')

        if polygon_for_bbox.crs != 'EPSG:4326':
            polygon_for_bbox = polygon_for_bbox.to_crs('EPSG:4326')
            print('geodataframe crs changed to EPSG:4326')

        lon_min, lat_min, lon_max, lat_max = polygon_for_bbox.total_bounds

    ## if only 1 var (short) is provided, change to list
    if not isinstance(variable, list):
        variable = [variable]

    ## Initiate list of empty datasets
    xarray_dict = dict()

    ## Loop through variables and pull data + place in dict
    for var in variable:
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
        print(f'finish in {round(end - start, 2)} seconds')

        # Append to dict of xr.datasets
        xarray_dict[var] = ds_subset

    return xarray_dict

def create_weightmap(xarray_dict, polygon, output_data_folder, weightmap_var = None):

    """
    :param geodataframe polygon: geodataframe polygon or multipolygon
    :param dict xarray_dict: dictionary of gridmet data. the output of get_gridmet_datasets()
    :param str output_data_folder: path to folder for function output
    :param str weightmap_var: variable used to make weightfile + naming of output e.i. weightmap_var = 'tmin'. If none, first variable of dict keys will be used.
    :return: output path to weightmap pickle file
    """
    if isinstance(polygon, geopandas.geodataframe.GeoDataFrame):
        print('polygon is geodataframe')
        pass
    elif os.path.isfile(polygon):
        print('polygon is path to shapefile')
        polygon = gpd.read_file(polygon)
    else:
        raise TypeError('polygon should str path to shapefile or geodataframe')

    ## Name output weightmap_file
    if weightmap_var is None:
        # If non var given
        weightmap_file = os.path.join(output_data_folder, 'grd2shp_weights.pickle')
        weightmap_var = list(xarray_dict.keys())[0]
     else:
        weightmap_file = os.path.join(output_data_folder, f'grd2shp_weights_{weightmap_var}.pickle')


    ## Produce weightmap
    start = time.perf_counter()
    weightmap = xa.pixel_overlaps(xarray_dict[weightmap_var], polygon)
    end = time.perf_counter()
    print(f'finished agg in {round(end - start, 2)} second(s)')

    with open(weightmap_file, 'wb') as file:
        pickle.dump(weightmap, file)

    return weightmap_file

def g2shp_regridding(xarray_dict, polygon, weightmap_file, output_data_folder, g2s_file_prefix,
                     g2s_time_var = 'day', g2s_lat_var = 'lat', g2s_lon_var = 'lon'):

    """
    :param dict xarray_dict: dictionary of gridmet data. the output of get_gridmet_datasets()
    :paran str weightmap_file: path to weight file for redridding
    :param str output_data_folder: path to folder for function output
    :param str weightmap_var: if len(xarray_dict) selected
    :param str g2s_time_var: time variable for g2shp initialization (i.e.g2s_time_var = 'day')
    :param str g2s_lat_var: latitude variable for g2shp initialization (i.e.g2s_time_var = 'lat')
    :param str g2s_lon_var: time variable for g2shp initialization (i.e.g2s_time_var = 'lon')
    :return:
    """

    ## params for g2s.initialise()
    # grab long name of data vars from dict and place in list
    vars_long_list = [list(xarray_dict[key])[0] for key in xarray_dict]
    # List of short names of data vars from dict keys
    vars_short_list = list(xarray_dict.keys())
    # Lsit of xarray.datasets
    vars_grd_list = list(xarray_dict.values())

    ## initialize Grd2ShpXagg()
    g2s = grd2shp_xagg.Grd2ShpXagg()
    g2s.initialize(
        grd = vars_grd_list,
        shp= polygon,
        wght_file= weightmap_file,
        time_var = g2s_time_var,
        lat_var = g2s_lat_var,
        lon_var = g2s_lon_var,
        var = vars_long_list,
        var_output = vars_short_list,
        ctype=0,
    )

    # Regridding
    start = time.perf_counter()
    g2s.run_weights()
    end = time.perf_counter()
    print('finished agg in {} second(s)'.format(round(end - start, 2)))
    print(g2s)

    ## output
    g2s.write_gm_file(opath=output_data_folder, prefix=g2s_file_prefix)

    ## Print location of output
    todays_date = str(datetime.datetime.now().strftime("%Y_%m_%d"))
    file_name = g2s_file_prefix + "climate_" + todays_date + ".nc"
    print('netcdf4 file path: ' + os.path.join(output_data_folder, file_name))

    return g2s

# Variables and run functions
if __name__ =='__main__':

    ## Variable definitions
    ### official list of variables needed for drb-inland-salinity model
    data_vars_shrt_all = ['tmmx', 'tmmn', 'pr', 'srad', 'vs','rmax','rmin','sph']
    ### drb catchment polygons
    gdf_nhru02_path = './data/nhru_01/nhru_01.shp'
#    gdf_nhru02_path = './data/nhru_02/nhru_02.shp'
    gdf = gpd.read_file(gdf_nhru02_path)
    ### date range
    start_date = '1979-01-01'
    end_date = '2021-01-01'
    ### lat lon
    lon_min, lat_min, lon_max, lat_max = gdf.total_bounds
    ### output folder
    output_path = './data/'

    xarray_dict = get_gridmet_datasets(variable = data_vars_shrt_all,
                                       start_date= start_date, end_date = end_date,
                                       polygon_for_bbox = gdf)
                                      # lat_min=round(lat_min,2), lon_min=round(lon_min,2),
                                      # lat_max=round(lat_max,2), lon_max=round(lon_max,2))
    
    ### Subset for streamlined testing 
    subset = {key: xarray_dict[key] for key in ['tmmx']}

    weight_map_path = create_weightmap(xarray_dict = xarray_dict,
                     polygon=gdf,
                     output_data_folder = output_path,
                     weightmap_var = 'tmmx')


    g2shp_regridding(xarray_dict= subset,
                     polygon=gdf,
                     weightmap_file= weight_map_path,
                     g2s_file_prefix='t_',
                     output_data_folder= output_path,
                     g2s_time_var = 'day',
                     g2s_lat_var = 'lat',
                     g2s_lon_var = 'lon')


######### Plot ##########

## Viewing output
#xr_mapped = xr.open_dataset('./data/all_climate_2022_01_25.nc', decode_times=False)
#gdf["tmmn"] = xr_mapped["tmmn"][:,0]
# gdf.plot(column = 'tmmn', legend = True)
# xr_mapped.isel(geomid=3).tmmn.plot()

