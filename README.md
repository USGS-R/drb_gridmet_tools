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


## Running re-gridding for the Delaware River Basin

`gridmet_split_script.py` processes the gridmet raster dataset values to the scale of the input multi-polygon shapefile.
