##############################
##############################
#### process_sensitivity.R 
#### Edward Lavender (el72@st-andrews.ac.uk)

#### This code:
# 1) Calculates species' sensitivity indicies given thermal niche parameters from (a) SST or (b) SBT
# ... mean thermal niche width (STR) over space
# ... SD thermal niche width (STR) over space
# ... mean thermal bias (STI - baseline) over space
# ... SD thermal bias (STI - baseline) over space 

#### Steps preceding this code:
# 1) Define species list (process_spptraits.R)
# 1) Temperature processing (process_temp.R)
# 2) SDMs processing (process_sdm_aquamaps.R)


##############################
##############################
#### Set up 

#### Wipe workspace
rm(list = ls())

#### Global param
root_sst <- "./data/temperature/sst"
root_sbt <- "./data/temperature/sbt"

#### Load data
spptraits <- readRDS("./data/spptraits.rds")
sst_historical <- raster::raster("./data/temperature/sst/historical/historical.asc")
sbt_historical <- raster::raster("./data/temperature/sbt/historical/historical.asc")


##############################
##############################
#### Define raster lists 

#### Method
# 1) For each species, load SDM
# 2) Across the species' range, define, for thermal parameters based on SST and SBT: 
# ... STI
# ... STR 
# ... STI - baseline 
# 3) Stack all of these rasters across species
# 4) Calculate summary statistics (e.g., mean thermal range etc.)

#### Define thermal niche parameters 
spptraits$index <- 1:nrow(spptraits)
spptraits$sst_str <- spptraits$sst_t90 - spptraits$sst_t10
spptraits$sbt_str <- spptraits$sbt_t90 - spptraits$sbt_t10

#### Set up cluster 
cl <- parallel::makeCluster(11L)
parallel::clusterExport(cl, varlist = c("sst_historical", "sbt_historical"))

#### Define a list of 'sensitivity' rasters for each species 
# [This takes ~ 8 minutes with 8 cores.]
sensitivity_by_spp <- 
  pbapply::pblapply(split(spptraits, spptraits$index), cl = cl, function(d){
    # d <- spptraits[2, ]
    sdm <- raster::raster(paste0("./data/sdm_aquamaps/", d$number.asc.file))
    # raster::plot(sdm); raster::plot(is.na(sdm)); raster::plot(sdm == 0); raster::plot(sdm == 1);
    sst_sti <- sbt_sti <- sst_str <- sbt_str <- sst_tb <- sbt_tb <- sdm
    index <- sdm == 1
    sst_sti[index] <- d$sst_t50
    sbt_sti[index] <- d$sbt_t50
    sst_str[index] <- d$sst_str
    sbt_str[index] <- d$sbt_str
    sst_tb[index]  <- d$sst_t50 - sst_historical[index]
    sbt_tb[index]  <- d$sbt_t50 - sbt_historical[index]
    # raster::plot(sst_sti)
    # raster::plot(sbt_sti)
    # raster::plot(sst_str)
    # raster::plot(sbt_str)
    # raster::plot(sst_tb)
    # raster::plot(sbt_tb)
    out <- list(sst_sti = sst_sti, 
                sbt_sti = sbt_sti, 
                sst_str = sst_str, 
                sbt_str = sbt_str, 
                sst_tb = sst_tb, 
                sbt_tb = sbt_tb)
    return(out)
  })


##############################
##############################
#### Process species' rasters 

#### Collect raster list for each metric 
sst_sti_by_spp <- pbapply::pblapply(sensitivity_by_spp, function(elm) elm$sst_sti)
sbt_sti_by_spp <- pbapply::pblapply(sensitivity_by_spp, function(elm) elm$sbt_sti)
sst_str_by_spp <- pbapply::pblapply(sensitivity_by_spp, function(elm) elm$sst_str)
sbt_str_by_spp <- pbapply::pblapply(sensitivity_by_spp, function(elm) elm$sbt_str)
sst_tb_by_spp  <- pbapply::pblapply(sensitivity_by_spp, function(elm) elm$sst_tb)
sbt_tb_by_spp  <- pbapply::pblapply(sensitivity_by_spp, function(elm) elm$sbt_tb)

#### Brick rasters for each metric 
# [This takes ~ 2.5 minutes]
t1_bk <- Sys.time()
sst_sti_bk <- raster::brick(sst_sti_by_spp)
sbt_sti_bk <- raster::brick(sbt_sti_by_spp)
sst_str_bk <- raster::brick(sst_str_by_spp)
sbt_str_bk <- raster::brick(sbt_str_by_spp)
sst_tb_bk  <- raster::brick(sst_tb_by_spp)
sbt_tb_bk  <- raster::brick(sbt_tb_by_spp)
t2_bk <- Sys.time()
difftime(t2_bk, t1_bk)

#### Summarises
# [This takes ~ 2.5 minutes]
t1_summary <- Sys.time()
# Mean 
sst_sti_mean <- raster::calc(sst_sti_bk, mean, na.rm = TRUE)
sbt_sti_mean <- raster::calc(sbt_sti_bk, mean, na.rm = TRUE)
sst_str_mean <- raster::calc(sst_str_bk, mean, na.rm = TRUE)
sbt_str_mean <- raster::calc(sbt_str_bk, mean, na.rm = TRUE)
sst_tb_mean  <- raster::calc(sst_tb_bk, mean, na.rm = TRUE)
sbt_tb_mean  <- raster::calc(sbt_tb_bk, mean, na.rm = TRUE)
# SD
sst_sti_sd <- raster::calc(sst_sti_bk, sd, na.rm = TRUE)
sbt_sti_sd <- raster::calc(sbt_sti_bk, sd, na.rm = TRUE)
sst_str_sd <- raster::calc(sst_str_bk, sd, na.rm = TRUE)
sbt_str_sd <- raster::calc(sbt_str_bk, sd, na.rm = TRUE)
sst_tb_sd  <- raster::calc(sst_tb_bk, sd, na.rm = TRUE)
sbt_tb_sd  <- raster::calc(sbt_tb_bk, sd, na.rm = TRUE)
t2_summary <- Sys.time()
difftime(t2_bk, t1_bk)


##############################
##############################
#### Save files

#### Save bricks 
raster::writeRaster(sst_sti_bk, "./data/sensitivity/sst/sst_sti_bk.tif")
raster::writeRaster(sbt_sti_bk, "./data/sensitivity/sbt/sbt_sti_bk.tif")
raster::writeRaster(sst_str_bk, "./data/sensitivity/sst/sst_str_bk.tif")
raster::writeRaster(sbt_str_bk, "./data/sensitivity/sbt/sbt_str_bk.tif")
raster::writeRaster(sst_tb_bk, "./data/sensitivity/sst/sst_tb_bk.tif")
raster::writeRaster(sbt_tb_bk, "./data/sensitivity/sbt/sbt_tb_bk.tif")

#### Save summaries
# SST, mean and SD thermal range and thermal bias 
raster::writeRaster(sst_sti_mean, "./data/sensitivity/sst/sst_sti_mean.asc")
raster::writeRaster(sst_sti_sd, "./data/sensitivity/sst/sst_sti_sd.asc")
raster::writeRaster(sst_str_mean, "./data/sensitivity/sst/sst_str_mean.asc")
raster::writeRaster(sst_str_sd, "./data/sensitivity/sst/sst_str_sd.asc")
raster::writeRaster(sst_tb_mean, "./data/sensitivity/sst/sst_tb_mean.asc")
raster::writeRaster(sst_tb_sd, "./data/sensitivity/sst/sst_tb_sd.asc")
# SBT, mean and SD thermal range and thermal bias
raster::writeRaster(sbt_sti_mean, "./data/sensitivity/sbt/sbt_sti_mean.asc")
raster::writeRaster(sbt_sti_sd, "./data/sensitivity/sbt/sbt_sti_sd.asc")
raster::writeRaster(sbt_str_mean, "./data/sensitivity/sbt/sbt_str_mean.asc")
raster::writeRaster(sbt_str_sd, "./data/sensitivity/sbt/sbt_str_sd.asc")
raster::writeRaster(sbt_tb_mean, "./data/sensitivity/sbt/sbt_tb_mean.asc")
raster::writeRaster(sbt_tb_sd, "./data/sensitivity/sbt/sbt_tb_sd.asc")


#### End of code. 
##############################
##############################