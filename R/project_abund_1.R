##############################
##############################
#### project_abund_1.R 
#### Edward Lavender (el72@st-andrews.ac.uk)

#### This code:
# 1) Projects changes in abundance under different climate scenarios
# ... Saving outputs on a by-species basis 
# ... project_abund_2.R then synthesises these across species. 

#### Steps preceding this code:
# 1) Define species list (process_spptraits.R)
# 1) Temperature processing (process_temp.R)
# 2) SDMs processing (process_sdm_aquamaps.R)


##############################
##############################
#### Set up 

#### Wipe workspace
rm(list = ls())

#### Essential packages
library(magrittr)
options(dplyr.summarise.inform = FALSE)

#### Global param
root_sst <- "./data/temperature/sst"
root_sbt <- "./data/temperature/sbt"

#### Load data
## spptraits
spptraits <- readRDS("./data/spptraits.rds")
## SST baselines and projections
sst_historical <- raster::raster("./data/temperature/sst/historical/historical.asc")
sst_mid_rcp45  <- raster::raster("./data/temperature/sst/mid_century/rcp45/rcp45.asc")
sst_mid_rcp85  <- raster::raster("./data/temperature/sst/mid_century/rcp85/rcp85.asc")
sst_late_rcp45 <- raster::raster("./data/temperature/sst/late_century/rcp45/rcp45.asc")
sst_late_rcp85 <- raster::raster("./data/temperature/sst/late_century/rcp85/rcp85.asc")
## SBT baselines and projections 
sbt_historical <- raster::raster("./data/temperature/sbt/historical/historical.asc")
sbt_mid_rcp45  <- raster::raster("./data/temperature/sbt/mid_century/rcp45/rcp45.asc")
sbt_mid_rcp85  <- raster::raster("./data/temperature/sbt/mid_century/rcp85/rcp85.asc")
sbt_late_rcp45 <- raster::raster("./data/temperature/sbt/late_century/rcp45/rcp45.asc")
sbt_late_rcp85 <- raster::raster("./data/temperature/sbt/late_century/rcp85/rcp85.asc")
## EEZs
eez <- sf::read_sf("./data/spatial/eez", "eez_v11") 
eez <- as(eez, "Spatial")


##############################
##############################
#### Define helper functions 

#### Define helper function to predict species' abundance
# Predicts relative abundance onto a blank raster, x,
# ... across the range of the species' distribution
# ... defined by a species distribution map (sdm), 
# ... using a raster of temperatures (temperature), 
# ... and species-specific thermal parameters (sti, str_sd). 
# ... This assumes that x, sdm and temperature are defined across
# ... the same extent with the same resolution. 
predict_abund <- function(x, sdm, temperature, sti, str_sd, scale = dnorm(sti, mean = sti, sd = str_sd)){
  x[sdm == 1] <- dnorm(temperature[sdm == 1], mean = sti, sd = str_sd)/scale
  return(x)
}

##### Define helper function to get sovereign states
# Since species distributions remain the same in all scenarios, 
# ... This only needs to be implemented once for each species for all SST scenarios
# ... and once again for all SBT scenarios (since available projections differ slightly).
assign_eez <- function(abund, eez){
  # Convert predicted abundance to xyz file 
  xyz <- data.frame(raster::rasterToPoints(abund))
  # Define SpatialPointsDataFrame with identical CRS as for EEZs
  sp::coordinates(xyz) <- c("x", "y")
  raster::crs(xyz) <- raster::crs(eez)
  # Assign EEZs to coordinates
  xyz$eez <- sp::over(xyz, eez)$SOVEREIGN1 
  xyz$eez <- factor(xyz$eez)
  # Return vector of EEZs
  return(xyz$eez)
}

#### Define helper function to predict total abundance in EEZ
# This function requires the magrittr to be loaded for the pipe operator.  
# Used to work out:
# ... starting abundance in each EEZ [for each scenario...]
# ... ending abundance in each EEZ   [for each scenario...]
sum_abund_in_eez <- function(abund, sovereign){
  
  # Convert predicted abundance to xyz file 
  xyz <- data.frame(raster::rasterToPoints(abund))
  names(xyz) <- c("long","lat","abund")
  # Add sovereign states
  xyz$eez <- sovereign
  
  # Calculate the overall abundance in each eez 
  abund_by_eez <- 
    xyz %>% 
    dplyr::group_by(eez) %>% 
    dplyr::filter(!is.na(eez)) %>% 
    dplyr::summarise(abund = sum(abund)) %>% 
    data.frame()
  
  # Return dataframe
  return(abund_by_eez)
}


##############################
##############################
#### Projections for each species

#### Duration
# This code takes ~ 3 hours with 8 cores. 

#### Set up cluster 
cl <- parallel::makeCluster(8L)
parallel::clusterEvalQ(cl, {library(magrittr)}) # (for %>% in sum_abund_in_eez())
vl <- c("predict_abund", "assign_eez", "sum_abund_in_eez", 
        "eez",
        "sst_historical",
        "sst_mid_rcp45",
        "sst_mid_rcp85",
        "sst_late_rcp45",
        "sst_late_rcp85",
        "sbt_historical",
        "sbt_mid_rcp45",
        "sbt_mid_rcp85", 
        "sbt_late_rcp45", 
        "sbt_late_rcp85")
parallel::clusterExport(cl = cl, varlist = vl)
spptraits$index <- 1:nrow(spptraits)
# spptraits <- spptraits[1:10, ]

### Loop over each species and make projections 
ab_in_eez_by_spp <- pbapply::pblapply(split(spptraits, spptraits$index), cl = cl, function(d){
  
  #### Define param for projections 
  
  ## Load raster for species 
  # d <- spptraits[1000, ]
  sdm <- raster::raster(paste0("./data/sdm_aquamaps/",  d$number.asc.file))
  # Define blank raster
  # ... cells where the species occurs will be replaced with abundance predictions
  # ... cells were the species does not occur will remain as NA 
  blank <- sdm
  blank <- raster::setValues(blank, NA)
  
  ## Define species thermal niche parameters 
  sti    <- d$sst_t50
  str_sd <- (d$sst_t90 - d$sst_t10)/(2*1.281560031) 
  denom  <- dnorm(sti, mean = sti, sd = str_sd)
  
  #### SST projections 
  
  ## Predict abundance from baseline temps
  ab_sst_historical <- predict_abund(x = blank,  sdm = sdm, temperature = sst_historical, 
                                     sti = sti, str_sd = str_sd, scale = denom)
  # Get vector of EEZs
  sst_sovereigns <- assign_eez(ab_sst_historical, eez)
  # Abundance in EEZ
  ab_sst_historical_by_eez <- sum_abund_in_eez(ab_sst_historical, sst_sovereigns)
  
  ## Predict abundance from mid-century scenarios
  # Predictions (maps)
  ab_sst_mid_rcp45 <- predict_abund(x = blank,  sdm = sdm, temperature = sst_mid_rcp45, 
                                    sti = sti, str_sd = str_sd, scale = denom)
  ab_sst_mid_rcp85 <- predict_abund(x = blank,  sdm = sdm, temperature = sst_mid_rcp85, 
                                    sti = sti, str_sd = str_sd, scale = denom)
  # Predictions (eezs)
  ab_sst_mid_rcp45_by_eez <- sum_abund_in_eez(ab_sst_mid_rcp45, sst_sovereigns)
  ab_sst_mid_rcp85_by_eez <- sum_abund_in_eez(ab_sst_mid_rcp85, sst_sovereigns)
  # Mid-century change
  ab_delta_sst_mid_rcp45 <- ab_sst_mid_rcp45 - ab_sst_historical
  ab_delta_sst_mid_rcp85 <- ab_sst_mid_rcp85 - ab_sst_historical 
  
  ## Predict abundance from late-century scenarios
  # Predictions 
  ab_sst_late_rcp45 <- predict_abund(x = blank,  sdm = sdm, temperature = sst_late_rcp45, 
                                     sti = sti, str_sd = str_sd, scale = denom)
  ab_sst_late_rcp85 <- predict_abund(x = blank,  sdm = sdm, temperature = sst_late_rcp85, 
                                     sti = sti, str_sd = str_sd, scale = denom)
  # Predictions (eezs)
  ab_sst_late_rcp45_by_eez <- sum_abund_in_eez(ab_sst_late_rcp45, sst_sovereigns)
  ab_sst_late_rcp85_by_eez <- sum_abund_in_eez(ab_sst_late_rcp85, sst_sovereigns)
  # Late-century change 
  ab_delta_sst_late_rcp45 <- ab_sst_late_rcp45 - ab_sst_historical
  ab_delta_sst_late_rcp85 <- ab_sst_late_rcp85 - ab_sst_historical 
  
  ## Group EEZ predictions 
  ab_sst_by_eez <- ab_sst_historical_by_eez
  colnames(ab_sst_by_eez)       <- c("eez", "ab_historical")
  ab_sst_by_eez$spp             <- d$spp
  ab_sst_by_eez$number.asc.file <- d$number.asc.file
  ab_sst_by_eez                 <- ab_sst_by_eez[, c("spp", "number.asc.file", "eez", "ab_historical")]
  ab_sst_by_eez$ab_mid_rcp45    <- ab_sst_mid_rcp45_by_eez$abund
  ab_sst_by_eez$ab_mid_rcp85    <- ab_sst_mid_rcp85_by_eez$abund
  ab_sst_by_eez$ab_late_rcp45   <- ab_sst_late_rcp45_by_eez$abund
  ab_sst_by_eez$ab_late_rcp85   <- ab_sst_late_rcp85_by_eez$abund
  
  #### SBT projections 
  
  ## Predict abundance from baseline temps
  ab_sbt_historical <- predict_abund(x = blank,  sdm = sdm, temperature = sbt_historical, 
                                     sti = sti, str_sd = str_sd, scale = denom)
  # Get vector of EEZs (SBT data are slightly different from SST data)
  sbt_sovereigns <- assign_eez(ab_sbt_historical, eez)
  # Abundance in each EEZ
  ab_sbt_historical_by_eez <- sum_abund_in_eez(ab_sbt_historical, sbt_sovereigns)
  
  ## Predict abundance from mid-century scenarios
  # Predictions 
  ab_sbt_mid_rcp45 <- predict_abund(x = blank,  sdm = sdm, temperature = sbt_mid_rcp45, 
                                    sti = sti, str_sd = str_sd, scale = denom)
  ab_sbt_mid_rcp85 <- predict_abund(x = blank,  sdm = sdm, temperature = sbt_mid_rcp85, 
                                    sti = sti, str_sd = str_sd, scale = denom)
  # Predictions (eezs)
  ab_sbt_mid_rcp45_by_eez <- sum_abund_in_eez(ab_sbt_mid_rcp45, sbt_sovereigns)
  ab_sbt_mid_rcp85_by_eez <- sum_abund_in_eez(ab_sbt_mid_rcp85, sbt_sovereigns)
  # Mid-century change
  ab_delta_sbt_mid_rcp45 <- ab_sbt_mid_rcp45 - ab_sbt_historical
  ab_delta_sbt_mid_rcp85 <- ab_sbt_mid_rcp85 - ab_sbt_historical 
  
  ## Predict abundance from late-century scenarios
  # Predictions (maps)
  ab_sbt_late_rcp45 <- predict_abund(x = blank,  sdm = sdm, temperature = sbt_late_rcp45, 
                                     sti = sti, str_sd = str_sd, scale = denom)
  ab_sbt_late_rcp85 <- predict_abund(x = blank,  sdm = sdm, temperature = sbt_late_rcp85, 
                                     sti = sti, str_sd = str_sd, scale = denom)
  # Predictions (eezs)
  ab_sbt_late_rcp45_by_eez <- sum_abund_in_eez(ab_sbt_late_rcp45, sbt_sovereigns)
  ab_sbt_late_rcp85_by_eez <- sum_abund_in_eez(ab_sbt_late_rcp85, sbt_sovereigns)
  # Late-century change 
  ab_delta_sbt_late_rcp45 <- ab_sbt_late_rcp45 - ab_sbt_historical
  ab_delta_sbt_late_rcp85 <- ab_sbt_late_rcp85 - ab_sbt_historical 
  
  ## Group EEZ predictions 
  ab_sbt_by_eez <- ab_sbt_historical_by_eez
  colnames(ab_sbt_by_eez)       <- c("eez", "ab_historical")
  ab_sbt_by_eez$spp             <- d$spp
  ab_sbt_by_eez$number.asc.file <- d$number.asc.file
  ab_sbt_by_eez                 <- ab_sbt_by_eez[, c("spp", "number.asc.file", "eez", "ab_historical")]
  ab_sbt_by_eez$ab_historical   <- ab_sbt_historical_by_eez$abund
  ab_sbt_by_eez$ab_mid_rcp45    <- ab_sbt_mid_rcp45_by_eez$abund
  ab_sbt_by_eez$ab_mid_rcp85    <- ab_sbt_mid_rcp85_by_eez$abund
  ab_sbt_by_eez$ab_late_rcp45   <- ab_sbt_late_rcp45_by_eez$abund
  ab_sbt_by_eez$ab_late_rcp85   <- ab_sbt_late_rcp85_by_eez$abund
  
  #### Write abundance-by-EEZ dataframes to file 
  # id <- substr(d$number.asc.file, 1, nchar(d$number.asc.file) - 4)
  # saveRDS(ab_sst_by_eez, paste0("./data/abundance/spp_change/sst/over_eez/", id, ".rds"))
  
  #### Write rasters to file 
  write_rasters <- FALSE
  if(write_rasters){
    # SST 
    raster::writeRaster(ab_delta_sst_mid_rcp45, 
                        paste0("./data/abundance/change/spp_specific/sst/mid_century/rcp45/", d$number.asc.file), 
                        overwrite = TRUE)
    raster::writeRaster(ab_delta_sst_mid_rcp85, 
                        paste0("./data/abundance/change/spp_specific/sst/mid_century/rcp85/", d$number.asc.file), 
                        overwrite = TRUE)
    raster::writeRaster(ab_delta_sst_late_rcp45, 
                        paste0("./data/abundance/change/spp_specific/sst/late_century/rcp45/", d$number.asc.file), 
                        overwrite = TRUE)
    raster::writeRaster(ab_delta_sst_late_rcp85, 
                        paste0("./data/abundance/change/spp_specific/sst/late_century/rcp85/", d$number.asc.file), 
                        overwrite = TRUE)
    # SBT 
    raster::writeRaster(ab_delta_sbt_mid_rcp45, 
                        paste0("./data/abundance/change/spp_specific/sbt/mid_century/rcp45/", d$number.asc.file), 
                        overwrite = TRUE)
    raster::writeRaster(ab_delta_sbt_mid_rcp85, 
                        paste0("./data/abundance/change/spp_specific/sbt/mid_century/rcp85/", d$number.asc.file), 
                        overwrite = TRUE)
    raster::writeRaster(ab_delta_sbt_late_rcp45, 
                        paste0("./data/abundance/change/spp_specific/sbt/late_century/rcp45/", d$number.asc.file), 
                        overwrite = TRUE)
    raster::writeRaster(ab_delta_sbt_late_rcp85, 
                        paste0("./data/abundance/change/spp_specific/sbt/late_century/rcp85/", d$number.asc.file), 
                        overwrite = TRUE)
  }
  #### Return abundance dataframes 
  out <- list(sst = ab_sst_by_eez, sbt = ab_sbt_by_eez)
  return(out)
})
parallel::stopCluster(cl)
# beepr::beep(10)

#### Save abundance by EEZ dataframes
ab_sst_in_eez <- 
  pbapply::pblapply(ab_in_eez_by_spp, function(elm) elm$sst) %>% dplyr::bind_rows()
ab_sbt_in_eez <- 
  pbapply::pblapply(ab_in_eez_by_spp, function(elm) elm$sbt) %>% dplyr::bind_rows()
saveRDS(ab_sst_in_eez, "./data/abundance/change/eez/ab_sst_in_eez.rds")
saveRDS(ab_sbt_in_eez, "./data/abundance/change/eez/ab_sbt_in_eez.rds")


#### End of code. 
##############################
##############################