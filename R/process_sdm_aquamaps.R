##############################
##############################
#### process_sdm_aquamaps.R
#### Edward Lavender (el72@st-andrews.ac.uk)

#### This code:
# 1) Processes aquamaps SDMs for selected species

#### Steps preceding this code: 
# 1) process_spptraits.R defines a list of modelled species
# ... This script is designed to be implemented part-way through process_spptraits.R
# ... ... using a reduce list of species. 
# 2) manual acquisition of aquamaps files (in /data-raw/) from MTB
# 3) Manual acquisition and processing of coastline data used to mask predictions (process_coastline.R)


##############################
##############################
#### Set up 

#### Wipe workspace
rm(list = ls())

#### Load intermediate spptraits file 
# This contains a reduced list of species for which to create processed predictions 
spptraits <- readRDS("./data-raw/spptraits_for_process_sdm_aquamaps.rds")
head(spptraits)

#### Load coastline 
coastline <- readRDS("./data/spatial/coastline/coastline.rds")


##############################
##############################
#### Processing

#### Method
# Open the map for each species
# If necessary, resample to match the resolution of the temperature files
# ... This is not necessary for 1 dg files. 
# Re-set 0 'occupancy' to NA
# ... This simplifies acquisition of thermal niche parameters
# ... and temperature projections, by ensuring that we only ever use
# ... areas where a species is predicted to be present, and other areas are ignored as NAs. 
# Remove any cells that overlap with the land via raster::mask()
# ... With this method, some 'cells' will overlap with the coastline and land
# ... but these cells this is because of the low resolution of climate models 
# Save processed SDM map. 

#### Loop over each species and save processed predictions 
# [This code takes ~ 3 minutes with 11 cores]
spptraits$index <- 1:nrow(spptraits)
cl <- parallel::makeCluster(11L)
parallel::clusterExport(cl = cl, varlist = "coastline")
pbapply::pblapply(split(spptraits, spptraits$index), cl = cl, function(d){
  # d <- spptraits[1, ]
  r <- raster::raster(paste0("./data-raw/sdm_aquamaps/", d$number.asc.file))
  # raster::plot(r)
  r[r == 0] <- NA 
  r <- raster::mask(r, coastline, updatevalue = NA)
  # raster::plot(r); raster::lines(coastline);
  raster::writeRaster(r, paste0("data/sdm_aquamaps/", d$number.asc.file), overwrite = TRUE)
  return(invisible())
})
parallel::stopCluster(cl)

#### Manual checks
# E.g. 1
raster::plot(raster::raster(paste0("data-raw/sdm_aquamaps/", spptraits$number.asc.file[1])))
raster::plot(raster::raster(paste0("data/sdm_aquamaps/", spptraits$number.asc.file[1])))
# E.g. 2
raster::plot(raster::raster(paste0("data-raw/sdm_aquamaps/", spptraits$number.asc.file[10])))
raster::plot(raster::raster(paste0("data/sdm_aquamaps/", spptraits$number.asc.file[10])))
# E.g. 3
raster::plot(raster::raster(paste0("data-raw/sdm_aquamaps/", spptraits$number.asc.file[1000])))
raster::plot(raster::raster(paste0("data/sdm_aquamaps/", spptraits$number.asc.file[1000])))


#### End of code.
##############################
##############################