# drb_gridmet_tools

Repository with functions to aggregate [gridmet climate raster data](https://www.climatologylab.org/gridmet.html) from pixel grid to vector aoi (hru polygon grid, polylines). This repository heavily relies on [grd2shp_xagg](https://github.com/rmcd-mscb/grd2shp_xagg) library, by [rmcd-mscb](https://github.com/rmcd-mscb)


## Selected gridMET variables: 

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


## running the code
`gridmet_split_script.py` processes the gridmet raster dataset values to the scale of the input multi-polygon shapefile.

`gridmet_aggregation_PRMS.py` processes the output of `gridmet_splot_script.py` and aggregates the PRMS_segid scale calculating an area weighted average. 

The snakemake pipeline specified in `Snakefile` will call functions from the .py scripts to aggregate gridmet data to the polygons in the geo file specified. Currently, it's set up to pull the gridmet data from the start of the data archive (1979-01-01) to the present day.

### running with Docker
```
docker run jsadler2/gridmet-agg:v0.1 snakemake
```

### running with Singularity
We have to specify the full path to `snakemake` because it is not part of the system PATH in the Singularity container
```
singularity pull jsadler2/gridmet-agg:v0.1
singularity exec gridmet-agg_v0_1.sif /opt/conda/snakemake
```
>>>>>>> [#7] update readme with container instructions
