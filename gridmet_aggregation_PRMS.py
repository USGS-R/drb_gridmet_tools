# libraries
import os
import geopandas as gpd
import xarray as xr
import time

# ncdf_to_gdf() function converts ncdf to dataset and merged with shpfile information (geometry + area)
def ncdf_to_gdf(ncdf_path, shp, left_on = 'geomid', right_on_index = True, gpkg_layer = None):

    """
    :param str ncdf_path: path to regridded ncdf file (output of g2shp_regridding())
    :param gpd.GeoDataFrame or str shp: shpfile or geopackage with polygons used for regridding
    :param str left_on: Selected ncdf timeseries data col to enable left merge with shp geodataframe. Default 'geomid'
    :param str right_on_index: Boolean. Default True. If True, index from shp var is used as right side col for left merge
    :param str gpkg_layer: Provide gpkg_layer name if shp is path to gpkg
    :return: Geodataframe of the original ncdf file with geospatial information (e.g. geom)
    """

    ## Read in
    xr_mapped_df = xr.open_dataset(ncdf_path, decode_times=True).to_dataframe().reset_index()
    if isinstance(shp, gpd.geodataframe.GeoDataFrame):
        gdf = shp
    elif os.path.isfile(shp):
        gdf_prms_path_edited = shp
        gdf = gpd.read_file(gdf_prms_path_edited, layer=gpkg_layer)
    else:
        print('shp must be path to geospatial file or a geodataframe')

    ## Merge ncdf w/ shapefile (the shpfile has area info) & convert to GeoDataFrame
    gridmet_drb_df = xr_mapped_df.merge(gdf, how ='left', left_on = left_on, right_index = right_on_index)
    gridmet_drb_gdf = gpd.GeoDataFrame(gridmet_drb_df)

    return gridmet_drb_gdf

# gridmet_prms_area_avg_agg() function aggregates to groupby_cols scale and outputs the aggregated dataset
def gridmet_prms_area_avg_agg(df, groupby_cols, val_colnames, wgt_col, output_path = None):
    """
    :param pd.DataFrame or gpd.GeoDataFrame df:
    :param  list groupby_cols: group by cols of df
    :param list val_colnames: list of columns to perform area-weighted average calculation on
    :param str wgt_col: weight column
    :param str output_path: if not None, path to save output df dataset in local dir.
    :return: df with gridmet metrics per PRMS_segid
    """

    ## Weighted average calc = sum((val[1] * weight[1]), ... , (val[n] * weight[n]) ) / sum(weight[1],...,weight[n])

    start = time.perf_counter()

    df = df.copy()
    ## 1. Calc value x weight col
    df.loc[:, val_colnames] = df[val_colnames].values * df[[wgt_col]].values

    ## 2. Create group by, by summing new value cols weighted by area
    df_grouped = df.groupby(groupby_cols).sum()

    ## 3. Get weighted average
    df_grouped.loc[:, val_colnames] = df_grouped[val_colnames].values / df_grouped[[wgt_col]].values

    df_final = df_grouped[val_colnames]

    if output_path:
        df_final.to_csv(output_path, sep = ',')
        print('Output saved in: ' + output_path)

    end = time.perf_counter()

    print('Time (sec) elapsed:', end - start)

    return df_final

# Define variables and run
if __name__ =='__main__':

    ## Variable definitions
    gdf_prms_path_edited = 'https://github.com/USGS-R/drb-network-prep/blob/940073e8d77c911b6fb9dc4e3657aeab1162a158/2_process/out/GFv1_catchments_edited.gpkg?raw=true'
    gdf = gpd.read_file(gdf_prms_path_edited, layer='GFv1_catchments_edited')
    gridmet_ncdf = './data/t_climate_2022_03_31.nc'
    data_vars_shrt_all = ['tmmx', 'tmmn', 'pr', 'srad', 'vs', 'rmax', 'rmin', 'sph']

    ## Create dataframe and merge with shapefile information
    gridmet_drb_gdf = ncdf_to_gdf(ncdf_path=gridmet_ncdf,
                                 shp = gdf,
                                 left_on = 'geomid',
                                 right_on_index = True)

    ## run aggregation on PRMS_segid and time
    df_agg = gridmet_prms_area_avg_agg(gridmet_drb_gdf,
                                         groupby_cols = ['PRMS_segid',"time"],
                                         val_colnames = data_vars_shrt_all,
                                         wgt_col='hru_area_m2',
                                         output_path= None)

    ## Uncomment to run
    # df_agg.reset_index().to_csv('../drb-inland-salinity-ml/1_fetch/in/grdmet_drb_agg_032321.csv', index = False)