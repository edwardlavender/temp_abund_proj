##############################
##############################
#### process_temp.R
#### Edward Lavender (el72@st-andrews.ac.uk)

#### This code:
# 1) Processes historical and future temperature projections 
# ... Gets CMIP5 projections from /data-raw
# ... Rotates onto 'standard' grid and forces extent from c(-180, 180, -90, 90)
# ... Masks predictions on land, using coastline data from NE (as for SDMs)
# ... Saves files 

#### Steps preceding this code:
# 1) Acquisition of raw temperature projections from NOAA's climate web portal for CMIP5 
# ... (https://psl.noaa.gov/ipcc/ocn/ccwp.html)
# 2) Acquisition of coastline data (from Natural Earth) and processing of these data (process_coastline.R)


##############################
##############################
#### Set up 

#### Wipe workspace
rm(list = ls())

#### Load coastline 
coastline <- readRDS("./data/spatial/coastline/coastline.rds")

#### Define root directories to raw and processed files 
root_raw_sst <- "./data-raw/temperature/sst/"
root_raw_sbt <- "./data-raw/temperature/sbt/"
root_pro_sst <- "./data/temperature/sst/"
root_pro_sbt <- "./data/temperature/sbt/"
# Check files 
list.files(root_raw_sst, recursive = TRUE)
list.files(root_raw_sbt, recursive = TRUE)

#### Define helper function for RCP projection processing 
# This function loads the projections for a specific variable 
# ... i.e., anomaly (the difference between the historical baseline and the scenario)
# ... rotates these, 
# ... forces the full extent (to match other spatial datasets), 
# ... adds the differences to the baseline, 
# ... masks these by the coastline
# ... and saves the processed raster in the appropriate directory. 
process_rcp <- function(baseline, 
                        load_file, 
                        varname = "anomaly", 
                        mask = coastline, 
                        save_file, 
                        overwrite = TRUE){
  # Load file 
  delta <- raster::raster(load_file, varname = varname)
  # Rotate
  delta <- raster::rotate(delta)
  # Force full extent 
  delta <- raster::extend(delta, raster::extent(baseline), value = NA)
  # Define predictions 
  pred <- baseline + delta
  # Mask predictions
  pred <- raster::mask(pred, mask)
  # Check predictions are sensible 
  raster::plot(pred)
  print(raster::cellStats(pred, range))
  # Save file 
  raster::writeRaster(pred, save_file, overwrite = overwrite)
  return(invisible(pred))
}


##############################
##############################
#### SST

##############################
#### Historical 

#### HADISST 
hadisst <- raster::raster(paste0(root_raw_sst, "historical/hadisst/SSTMean.asc"))
hadisst <- hadisst/100
hadisst <- raster::mask(hadisst, coastline)
raster::plot(hadisst)
raster::cellStats(hadisst, range)
raster::writeRaster(hadisst, paste0(root_pro_sst, "historical/hadisst/hadisst.asc"), overwrite = TRUE)

#### Historical climatology 
historical <- raster::raster(paste0(root_raw_sst, "historical/historical.nc"), 
                             varname = "histclim")
historical <- raster::rotate(historical)
historical <- raster::extend(historical, raster::extent(-180, 180, -90, 90), value = NA)
historical <- raster::mask(historical, coastline)
raster::plot(historical)
raster::cellStats(historical, range)
raster::writeRaster(historical, paste0(root_pro_sst, "historical/historical.asc"), overwrite = TRUE)


##############################
#### Mid-century

#### RCP 4.5
sst_mid_rcp45 <- process_rcp(baseline = historical, 
                             load_file = paste0(root_raw_sst, "mid_century/rcp45/rcp45.nc"), 
                             save_file = paste0(root_pro_sst, "mid_century/rcp45/rcp45.asc"))

#### RCP 8.5
sst_mid_rcp85 <- process_rcp(baseline = historical, 
                             load_file = paste0(root_raw_sst, "mid_century/rcp85/rcp85.nc"), 
                             save_file = paste0(root_pro_sst, "mid_century/rcp85/rcp85.asc"))


##############################
#### Late-century 

#### RCP 4.5
sst_late_rcp45 <- process_rcp(baseline = historical, 
                              load_file = paste0(root_raw_sst, "late_century/rcp45/rcp45.nc"), 
                              save_file = paste0(root_pro_sst, "late_century/rcp45/rcp45.asc"))

#### RCP 8.5
sst_late_rcp85 <- process_rcp(baseline = historical, 
                              load_file = paste0(root_raw_sst, "late_century/rcp85/rcp85.nc"), 
                              save_file = paste0(root_pro_sst, "late_century/rcp85/rcp85.asc"))


##############################
##############################
#### SBT

##############################
#### Historical 

#### Historical climatology 
historical <- raster::raster(paste0(root_raw_sbt, "historical/historical.nc"), 
                             varname = "histclim")
historical <- raster::rotate(historical)
historical <- raster::extend(historical, raster::extent(-180, 180, -90, 90), value = NA)
historical <- raster::mask(historical, coastline)
raster::plot(historical)
raster::cellStats(historical, range)
raster::writeRaster(historical, paste0(root_pro_sbt, "historical/historical.asc"), overwrite = TRUE)


##############################
#### Mid-century

#### RCP 4.5
sbt_mid_rcp45 <- process_rcp(baseline = historical, 
                             load_file = paste0(root_raw_sbt, "mid_century/rcp45/rcp45.nc"), 
                             save_file = paste0(root_pro_sbt, "mid_century/rcp45/rcp45.asc"))

#### RCP 8.5
sbt_mid_rcp85 <- process_rcp(baseline = historical, 
                             load_file = paste0(root_raw_sbt, "mid_century/rcp85/rcp85.nc"), 
                             save_file = paste0(root_pro_sbt, "mid_century/rcp85/rcp85.asc"))


##############################
#### Late-century 

#### RCP 4.5
sbt_late_rcp45 <- process_rcp(baseline = historical, 
                              load_file = paste0(root_raw_sbt, "late_century/rcp45/rcp45.nc"), 
                              save_file = paste0(root_pro_sbt, "late_century/rcp45/rcp45.asc"))

#### RCP 8.5
sbt_late_rcp85 <- process_rcp(baseline = historical, 
                              load_file = paste0(root_raw_sbt, "late_century/rcp85/rcp85.nc"), 
                              save_file = paste0(root_pro_sbt, "late_century/rcp85/rcp85.asc"))


#### End of code. 
##############################
##############################