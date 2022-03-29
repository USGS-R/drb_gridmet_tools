# drb_gridmet_tools

Repository with functions to aggregate raster data from pixel grid to vector aoi (hru polygon grid, polylines). This repository heavily relies on [grd2shp_xagg](https://github.com/rmcd-mscb/grd2shp_xagg) library, by [rmcd-mscb](https://github.com/rmcd-mscb)


Run gridmet_split_script.py to get output ncdf at polygon grid.
Run gridmet_aggregation.py to aggregate to coarser polygon scale.For example, hru_catchment polygon to prms_segment. 

