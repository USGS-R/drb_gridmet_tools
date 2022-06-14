import geopandas as gpd
import xarray as xr
from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
from gridmet_split_script import get_gridmet_datasets, create_weightmap, g2shp_regridding
from gridmet_aggregation_PRMS import ncdf_to_gdf, gridmet_prms_area_avg_agg
import requests
from datetime import datetime

todays_date = datetime.today().strftime('%Y_%m_%d')

final_files = []

nc_file_path = "drb-gridmet/{fabric_id}/{todays_date}/{run_prefix}_climate_{todays_date}.nc"
final_nc_file_path = nc_file_path.format(fabric_id = config['fabric_id'],
                                         todays_date=todays_date,
                                         run_prefix=config['run_prefix'])

seg_file_path = "drb-gridmet/{fabric_id}/{todays_date}/{run_prefix}_climate_{todays_date}_segments.csv"
final_seg_file_path = seg_file_path.format(fabric_id = config['fabric_id'],
                                           todays_date=todays_date,
                                           run_prefix=config['run_prefix'])

final_files.append(final_nc_file_path)

if config['fabric_id'] == 'nhm':
    final_files.append(final_seg_file_path)

if config['use_S3']:
    S3 = S3RemoteProvider(keep_local=True)
    final_files = [S3.remote(f) for f in final_files]
    nc_file_path = S3.remote(nc_file_path)
    seg_file_path = S3.remote(seg_file_path)


rule all:
    input:
        final_files


rule make_weight_map:
    output:
        "drb-gridmet/{fabric_id}/grd2shp_weights.pickle"
    run:
        gdf = gpd.read_file(config['catchment_file_path'])
        # getting just one date and one variable to make the weight map.
        # the same weight map applies to all dates and all variables
        data_dict = get_gridmet_datasets(variable="tmmn",
                                         start_date="2001-01-01",
                                         end_date="2001-01-02",
                                         polygon_for_bbox=gdf)
        create_weightmap(xarray_dict=data_dict,
                         polygon=gdf,
                         output_data_folder = os.path.split(output[0])[0])


rule aggregate_gridmet_to_polygons_one_var:
    input:
        "drb-gridmet/{fabric_id}/grd2shp_weights.pickle"
    output:
        "drb-gridmet/{fabric_id}/{todays_date}/{run_prefix}_var_{variable}_climate_{todays_date}.nc"
    run:
        gdf = gpd.read_file(config['catchment_file_path'])
        data_dict = get_gridmet_datasets(variable=wildcards.variable,
                                         start_date=config.get('start_date', "1979-01-01"),
                                         end_date=config.get('end_date', todays_date.replace("_", "-")),
                                         polygon_for_bbox=gdf)
        g2shp_regridding(xarray_dict=data_dict,
                         polygon=gdf,
                         weightmap_file= input[0],
                         g2s_file_prefix=f'{wildcards.run_prefix}_var_{wildcards.variable}_',
                         output_data_folder= os.path.split(output[0])[0],
                         g2s_time_var = 'day',
                         g2s_lat_var = 'lat',
                         g2s_lon_var = 'lon')


rule gather_gridmets:
    input:
        expand("drb-gridmet/{fabric_id}/{todays_date}/{run_prefix}_var_{variable}_climate_{todays_date}.nc",
               fabric_id=config['fabric_id'],
               todays_date=todays_date,
               run_prefix=config['run_prefix'],
               variable=config['data_vars'])
    output:
        nc_file_path
    run:
        ds_list = [xr.open_dataset(nc_file) for nc_file in input]
        xr.merge(ds_list).to_netcdf(output[0])
        
        
rule aggregate_gridmet_polygons_to_flowlines:
    input:
        "drb-gridmet/{fabric_id}/{todays_date}/{run_prefix}_climate_{todays_date}.nc"
    output:
        seg_file_path
    run:
        gdf = gpd.read_file(config['catchment_file_path'])
        gridmet_drb_gdf = ncdf_to_gdf(ncdf_path=input[0],
                                      shp = gdf,
                                      left_on = 'geomid',
                                      right_on_index = True)
        df_agg = gridmet_prms_area_avg_agg(gridmet_drb_gdf,
                                           groupby_cols = ['PRMS_segid',"time"],
                                           val_colnames = config['data_vars'],
                                           wgt_col='hru_area_m2',
                                           output_path= output[0])

