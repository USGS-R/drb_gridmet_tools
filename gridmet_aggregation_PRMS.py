import geopandas as gpd
import geopandas.geodataframe
import xarray as xr

def ncdf_to_df(ncdf_path, shp, left_on = 'geomid', right_on_index = True, gpkg_layer = None):

    """
    :param str ncdf_path: path to regridded ncdf file (output of g2shp_regridding())
    :param gpd.GeoDataFrame or str shp: shpfile or geopackage with polygons used for regridding
    :param str left_on: Selected ncdf timeseries data col to enable left merge with shp geodataframe. Default 'geomid'
    :param str right_on_index: Boolean. Default True. If True, index from shp var is used as right side col for left merge
    :param str gpkg_layer: Provide gpkg_layer name if shp is path to gpkg
    :return: Geodataframe of the original ncdf file with geospatial information (e.g. geom)
    """

    # Read in
    xr_mapped_df = xr.open_dataset(ncdf_path, decode_times=True).to_dataframe().reset_index()

    if isinstance(shp, geopandas.geodataframe.GeoDataFrame):
        gdf = shp
    elif os.path.isfile(shp):
        gdf_prms_path_edited = shp
        gdf = gpd.read_file(gdf_prms_path_edited, layer=gpkg_layer)

    else:
        print('shp must be path to geospatial file or a geodataframe')

    ## Merge & convert to GeoDataFrame
    gridmet_drb_df = xr_mapped_df.merge(gdf, how ='left', left_on = left_on, right_index = right_on_index)

    gridmet_drb_gdf = gpd.GeoDataFrame(gridmet_drb_df)

    return gridmet_drb_gdf

## function #2
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
    # Calc value x weight col
    for col in val_colnames:
           df[f'{col}_x_area'] = df[col]*df[wgt_col]

    # Create group by by suming new value cols weighted by are
    df_grouped = df.groupby(groupby_cols).agg(area_m2_sum = (wgt_col, 'sum'),
                                              tmmx_wgtd_sum = ('tmmx_x_area', 'sum'),
                                              tmmn_wgtd_sum = ('tmmn_x_area', 'sum'),
                                              pr_wgtd_sum = ('pr_x_area', 'sum'),
                                              srad_wgtd_sum = ('srad_x_area', 'sum'),
                                              vs_wgtd_sum = ('vs_x_area', 'sum'),
                                              rmax_wgtd_sum = ('rmax_x_area', 'sum'),
                                              rmin_wgtd_sum = ('rmin_x_area', 'sum'),
                                              sph_wgtd_sum = ('sph_x_area', 'sum')
                                              )
    # Get weighted average
    for col in val_colnames:
        df_grouped[f'{col}_wgtd_avg'] = df_grouped[f'{col}_wgtd_sum'] / df_grouped['area_m2_sum']

    # remove interim col
    df_final = df_grouped.loc[:, ~(df_grouped.columns.str.contains('wgtd_sum'))]

    if output_path is not None:
        df_final.to_csv(output_path, sep = ',')
        print('Output saved in: ' + output_path)

    end = time.perf_counter()

    print('Time (sec) ellapsed:' + end - start)

    return df_final

# Define variables and run functions
if __name__ =='__main__':

    ## Variable definitions
    gdf_prms_path_edited = './data/GFv1_catchments_edited.gpkg'
    gdf = gpd.read_file(gdf_prms_path_edited, layer = 'GFv1_catchments_edited')
    data_vars_shrt_all = ['tmmx', 'tmmn', 'pr', 'srad', 'vs', 'rmax', 'rmin', 'sph']

    ## Create dataframe and merge with shapefile information
    gridmet_drb_gdf = ncdf_to_df(ncdf_path='./data/t_climate_2022_03_09.nc',
                                 shp = gdf,
                                 left_on = 'geomid',
                                 right_on_index = True)

    ## run aggregation on PRMS_segid and time
    df_agg = gridmet_prms_area_avg_agg(gridmet_drb_gdf,
                                         groupby_cols = ['PRMS_segid',"time"],
                                         val_colnames = data_vars_shrt_all,
                                         wgt_col='hru_area_m2',
                                         output_path= 'data/gridmet_drb_agg.csv')

