# drb_gridmet_tools

	Repository with functions to aggregate [gridmet climate raster data](https://www.climatologylab.org/gridmet.html) from pixel grid to vector aoi (hru polygon grid, polylines). This repository heavily relies on the [grd2shp_xagg](https://github.com/rmcd-mscb/grd2shp_xagg) library, by [rmcd-mscb](https://github.com/rmcd-mscb), which in turn relies oni [xagg](https://github.com/ks905383/xagg).

# Accessing the re-gridded files
## On Caldera
On Caldera, re-gridded file paths are structured like this:`/caldera/projects/usgs/water/impd/pump/gridmet/drb_gridmet_tools/drb-gridmet/{fabric}/{run_date}/drb_climate_{run_date}.nc`. For example: ``/caldera/projects/usgs/water/impd/pump/gridmet/drb_gridmet_tools/drb-gridmet/nhm/2022_06_14/drb_climate_2022_06_14.nc`.

## On S3
As of 6/2022, the results of the workflow are, by default, also stored S3 in the `drb-gridmet` bucket. For example, `s3://drb-gridmet/nhm/2022_06_14/drb_climate_2022_06_14.nc`.:w

# Running re-gridding for the Delaware River Basin via Snakemake 
The re-gridding process for the Delaware River Basin has been run on USGS's Tallgrass via Singularity. It has been run for the National Hydrologic Model (NHM) fabric and the National Hydrographic Database (NHD) fabric. It should be able to be run via Docker as well as Singularity, but it has not been tried yet. Instructions for running and modifying the pipeline are below.

_Parallelization in Snakemake_
The Snakemake workflow parallelizes the re-gridding of the 8 Gridmet variables. As long as you provide at least 8 cores, these tasks will run in parallel.

_S3_
By default, the results of the workflow are stored locally and on S3 in the `drb-gridmet` bucket. If you would like to only store the resulting files locally, change the `use_S3` option in the config file to `False`.

NOTE: Although this workflow was first run on DRB catchments, it should be able to be run on any set of polygons in the conterminous US. This however, has not be tested.

## Run via Singularity
To run the workflow with Singularity, you can either use the image already in Caldera via Tallgrass, or pull the image into a new directory.

### Option A. Use the existing image and cloned repo
1. move to the correct directory

```
cd /caldera/projects/usgs/water/impd/pump/gridmet/drb_gridmet_tools/
```

2. Decide on which fabric you want to use. You will use a different `config` file depending on which fabric you use: `config_nhm.yml` for the NHM fabric and `config_nhd.yml` for the `nhd` fabric.
3. [Optional] edit the options in the config file 
4. Run the workflow
You can run the workflow on Tallgrass either in batch mode or interactively. Subtite your HPC account (may not be 'iidd') and desired config file (may not be `config_nhd.yml`)

You can run the workflow via `sbatch` 
```
sbatch -A iidd slurm/launch_snakemake.slurm config_nhd.yml
```
It may be helpful instead to run it interactively

```
salloc -N 1 -n 8 -t 10:00:00 -p cpu -A iidd
module load singularity
singularity exec-agg_v0.3.sif /opt/conda/bin/snakemake -j --configfile config_nhd.yml 
```

### Option B. Executing in a different directory 
1. Clone the repo and move to the `drb_gridmet_tools` directory
```
git clone git@github.com:USGS-R/drb_gridmet_tools.git
cd drb_gridmet_tools
```
2. Pull down Docker image into a Singularity file
```
singularity pull docker://jsadler2/gridmet-agg:v0.3
```
3. Do Steps 2-4 from Option A.


# Selected gridMET variables: 

tmmx:
* Description: Daily Maximum Temperature (2m)\
* Units: degrees Fahrenheit

tmmn: 
* Description: Daily Minimum Temperature (2m)\
* Units: degrees Fahrenheit
       
pr: 
* Description: Daily Accumulated Precipitation\
* Units: inch

srad: 
* Description: Daily Mean downward shortwave radiation at surface\
* Units: W m-2
      
vs:
* Description: Daily Mean Wind Speed (10m) 
* Units: m/s

rmax:
* Description: Daily Maximum Relative Humidity (2m) 
* Units: %
      
rmin:
* Description: Daily Minimum Relative Humidity (2m) 
* Units: %

sph:
* Description: Daily mean specific humidity (2m)
* Units: kg/kg



