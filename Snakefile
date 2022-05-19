import pickle
import geopandas as gpd
from gridmet_split_script import get_gridmet_datasets, create_weightmap
from gridmet_aggregation_PRMS import ncdf_to_gdf, gridmet_prms_area_avg_agg

run_date = "2022_04_07"
out_dir = "data/out"
data_dir = "data/"
drb_hru_file = "data/GFv1_catchments_edited.gpkg"

rule all:
    input:
        f"{out_dir}/drb_climate_{run_date}.nc",
        f"{out_dir}/drb_climate_{run_date}_segments.csv"


rule fetch_drb_catchments:
    output:
        drb_hru_file
    shell:
        "wget https://github.com/USGS-R/drb-network-prep/blob/940073e8d77c911b6fb9dc4e3657aeab1162a158/2_process/out/GFv1_catchments_edited.gpkg?raw=true -O {output}"


rule make_dataset_dict:
    input:
        drb_hru_file,
    params:
        data_vars = ['tmmx', 'tmmn', 'pr', 'srad', 'vs','rmax','rmin','sph'],
        start_date = "1979-01-01-",
        end_date = run_date.replace("_", "-")
    output:
        "{outdir}/dataset_dict_{run_date}.pickle"
    run:
        gdf = gpd.read_file(input[0], layer="GFv1_catchments_edited")
        data_dict = get_gridmet_datasets(variable = params.data_vars,
                                         start_date = params.start_date,
                                         end_date = params.end_date,
                                         polygon_for_bbox = gdf)
        with open(output[0], "wb") as f:
            pickle.dump(data_dict, f)


rule make_weight_map:
    input:
        drb_hru_file,
        f"{{outdir}}/dataset_dict_{run_date}.pickle"
    output:
        "{outdir}/grd2shp_weights.pickle"
    run:
        gdf = gpd.read_file(input[0], layer="GFv1_catchments_edited")
        with open(input[1], "rb") as f:
            xarray_dict = pickle.load(f)
        create_weightmap(xarray_dict = xarray_dict,
                         polygon=gdf,
                         output_data_folder = os.path.split(output[0])[0],
                         weightmap_var = 'tmmn')


rule aggregate_gridmet_to_polygons:
    input:
        drb_hru_file,
        "{outdir}/dataset_dict_{run_date}.pickle",
        "{outdir}/grd2shp_weights.pickle"
    output:
        "{outdir}/drb_climate_{run_date}.nc"
    run:
        gdf = gpd.read_file(input[0], layer="GFv1_catchments_edited")
        with open(input[1], "rb") as f:
            xarray_dict = pickle.load(f)
        g2shp_regridding(xarray_dict= xarray_dict,
                         polygon=gdf,
                         weightmap_file= input[2],
                         g2s_file_prefix='drb_',
                         output_data_folder= os.path.split(output[0])[0],
                         g2s_time_var = 'day',
                         g2s_lat_var = 'lat',
                         g2s_lon_var = 'lon')
        
        
rule make_nc_gdf:
    input:
        drb_hru_file,
        "{outdir}/drb_climate_{run_date}.nc"
    output:
        "{outdir}/drb_climate_{run_date}_gdf.pickle"
    run:
        gdf = gpd.read_file(input[0], layer="GFv1_catchments_edited")
        gridmet_drb_gdf = ncdf_to_gdf(ncdf_path=input[1],
                                      shp = gdf,
                                      left_on = 'geomid',
                                      right_on_index = True)
        with open(output[0], "wb") as f:
            pickle.dump(gridmet_drb_gdf, f)


rule aggregate_gridmet_polygons_to_flowlines:
    input:
        "{outdir}/drb_climate_{run_date}_gdf.pickle",
        "{outdir}/drb_climate_{run_date}.nc"
    params:
        data_vars = ['tmmx', 'tmmn', 'pr', 'srad', 'vs','rmax','rmin','sph'],
    output:
        "{outdir}/drb_climate_{run_date}_segments.csv"
    run:
        with open(input[0], "rb") as f:
            gridmet_drb_gdf = pickle.load(f)
        df_agg = gridmet_prms_area_avg_agg(gridmet_drb_gdf,
                                           groupby_cols = ['PRMS_segid',"time"],
                                           val_colnames = params.data_vars,
                                           wgt_col='hru_area_m2',
                                           output_path= output[0])

