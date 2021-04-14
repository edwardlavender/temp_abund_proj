##############################
##############################
#### process_eez.R
#### Edward Lavender (el72@st-andrews.ac.uk)

#### This code:
# 1) Gets EEZ areas and latitudinal mid-points 

#### Steps preceding this code: 
# 1) Definition of global extent/modelling resolution 
# ... e.g., extracted from processed temperature data from process_temp.R


##############################
##############################
#### Set up 

#### Wipe workspace
rm(list = ls())

#### Load data
# Load SST data as a 'blank' raster that gives dimensions of area
sst_historical <- raster::raster("./data/temperature/sst/historical/historical.asc")
# Load EEZ data 
eez <- sf::read_sf("./data/spatial/eez", "eez_v11") 
eez <- as(eez, "Spatial")


##############################
##############################
#### Define a mask over EEZs 

#### Make a mask for temperature projections
# ... We will focus on areas within EEZs
# ... We will make a raster which defines, for each cell, whether or not it is inside an EEZ
# ... and then use this to mask temperature/abundance projections to focus on coastal areas. 
# ... This is much quicker than trying to mask files
# ... directly using the EEZ data. 

cover <- sst_historical
cover <- data.frame(raster::rasterToPoints(cover))
colnames(cover) <- c("x", "y", "z")
cover$z <- 0
sp::coordinates(cover) <- c("x", "y")
raster::crs(cover) <- raster::crs(eez)
cover$eez <- sp::over(cover, eez)$SOVEREIGN1
head(cover)
cover$z <- ifelse(!is.na(cover$eez), 0, NA)
cover$eez <- NA
cover <- raster::rasterFromXYZ(cover, res = 1)
cover <- raster::extend(cover, raster::extent(sst_historical))
raster::plot(cover)
raster::writeRaster(cover, "./data/spatial/eez/eez_mask.asc", overwrite = TRUE)


##############################
##############################
#### Obtain EEZ in each raster cell 

#### Convert raster to SpatialPointsDataFrame 
globe <- sst_historical
globe <- raster::setValues(globe, 1)
globe <- data.frame(raster::rasterToPoints(globe))
sp::coordinates(globe) <- c("x", "y")
raster::crs(globe) <- raster::crs(eez)

#### Get EEZ for each raster cell 
# It is necessary to loop over each EEZ, 
# ... rather than use sp::over() on all EEZs, 
# ... because some 1 dg cells contain multiple EEZs. 
# This code takes a long time to run:
# ... To just calculate the total number of cells, this takes 3.5 minutes
# ... But to get the mid-latitudinal coordinates of each state takes 3.5 hours 
eez_vec <- unique(eez$SOVEREIGN1)
eez_vec <- eez_vec[!is.na(eez_vec)]
eez_n_cell_by_eez <- 
  pbapply::pblapply(eez_vec, function(state){
    # state <- "Guinea-Bissau"
    eez_for_state <- subset(eez, SOVEREIGN1 == state)
    state_in_cell <- sp::over(globe, eez_for_state)$SOVEREIGN1 
    n_cell_tot <- sum(state_in_cell == state, na.rm = TRUE)
    lat_mid <- rgeos::gCentroid(eez_for_state, byid = FALSE)
    lat_mid <- sp::coordinates(lat_mid)[2]
    d <- data.frame(eez = state, 
                    lat_mid = lat_mid,
                    n_cell_tot = n_cell_tot)
    return(d)
  })
eez_n_cell <- do.call(rbind, eez_n_cell_by_eez)
head(eez_n_cell)

#### Save dataframe 
save <- FALSE
if(save){
  saveRDS(eez_n_cell, "./data/spatial/eez/eez_n_cell_tot.rds")
}


#### End of code 
##############################
##############################