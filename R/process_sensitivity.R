##############################
##############################
#### process_sensitivity.R 
#### Edward Lavender (el72@st-andrews.ac.uk)

#### This code:
# 1) Calculates species' sensitivity indices given thermal niche parameters from (a) SST or (b) SBT
# ... mean thermal niche width (STR) over space
# ... mean thermal bias (STI - baseline) over space
# ... Variability in thermal niche parameters over space 

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
source("./R/helpers.R")

#### Load data
spptraits <- readRDS("./data/spptraits.rds")
sst_historical <- raster::raster("./data/temperature/sst/historical/historical.asc")
sbt_historical <- raster::raster("./data/temperature/sbt/historical/historical.asc")


##############################
##############################
#### Brick species' rasters 

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
cl <- parallel::makeCluster(12L)
parallel::clusterExport(cl, varlist = c("sst_historical", "sbt_historical"))

#### Define a list of 'sensitivity' rasters for each species 
# [This takes ~ 8 minutes with 8 cores.]
sensitivity_by_spp <- 
  pbapply::pblapply(split(spptraits, spptraits$index), cl = cl, function(d){
    # d <- spptraits[2, ]
    sdm <- raster::raster(paste0("./data/sdm_aquamaps/maps_occ/", d$spp_key_asc))
    # raster::plot(sdm); raster::plot(is.na(sdm)); raster::plot(sdm == 0); raster::plot(sdm > 0);
    sst_sti <- sbt_sti <- sst_str <- sbt_str <- sst_tb <- sbt_tb <- sdm
    index <- sdm > 0
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
parallel::stopCluster(cl)

#### Collect raster list for each metric 
sst_str_by_spp <- pbapply::pblapply(sensitivity_by_spp, function(elm) elm$sst_str)
sbt_str_by_spp <- pbapply::pblapply(sensitivity_by_spp, function(elm) elm$sbt_str)
sst_tb_by_spp  <- pbapply::pblapply(sensitivity_by_spp, function(elm) elm$sst_tb)
sbt_tb_by_spp  <- pbapply::pblapply(sensitivity_by_spp, function(elm) elm$sbt_tb)

#### Subset raster lists for testing speed
quick <- FALSE
if(quick){
  sst_str_by_spp <- sst_str_by_spp[1:10]
  sbt_str_by_spp <- sst_str_by_spp[1:10]
  sst_tb_by_spp  <- sst_str_by_spp[1:10]
  sbt_tb_by_spp  <- sst_str_by_spp[1:10]
}
  
#### Brick rasters for each metric 
# [This takes ~ 6.5 minutes]
t1_bk <- Sys.time()
sst_str_bk <- raster::brick(sst_str_by_spp)
sbt_str_bk <- raster::brick(sbt_str_by_spp)
sst_tb_bk  <- raster::brick(sst_tb_by_spp)
sbt_tb_bk  <- raster::brick(sbt_tb_by_spp)
t2_bk <- Sys.time()
difftime(t2_bk, t1_bk)
beepr::beep(10)


##############################
##############################
#### Summarise bricks 

#### Summarises
# [This takes ~ 7 minutes]
t1_summary <- Sys.time()
# Mean 
sst_str_mean <- raster::calc(sst_str_bk, mean, na.rm = TRUE)
sbt_str_mean <- raster::calc(sbt_str_bk, mean, na.rm = TRUE)
sst_tb_mean  <- raster::calc(sst_tb_bk, mean, na.rm = TRUE)
sbt_tb_mean  <- raster::calc(sbt_tb_bk, mean, na.rm = TRUE)
# IQR
# sst_str_iqr <- raster::calc(sst_str_bk, IQR, na.rm = TRUE)
# sbt_str_iqr <- raster::calc(sbt_str_bk, IQR, na.rm = TRUE)
sst_tb_iqr    <- raster::calc(sst_tb_bk, IQR, na.rm = TRUE)
sbt_tb_iqr    <- raster::calc(sbt_tb_bk, IQR, na.rm = TRUE)
t2_summary    <- Sys.time()
difftime(t2_summary, t1_summary)


##############################
##############################
#### Save files

#### Save bricks 
# (This takes about 8 minutes)
t1_write <- Sys.time()
raster::writeRaster(sst_str_bk, "./data/sensitivity/sst/sst_str_bk.tif", overwrite = TRUE)
raster::writeRaster(sbt_str_bk, "./data/sensitivity/sbt/sbt_str_bk.tif", overwrite = TRUE)
raster::writeRaster(sst_tb_bk, "./data/sensitivity/sst/sst_tb_bk.tif", overwrite = TRUE)
raster::writeRaster(sbt_tb_bk, "./data/sensitivity/sbt/sbt_tb_bk.tif", overwrite = TRUE)
t2_write <- Sys.time()
difftime(t2_write, t1_write)

#### Save summaries
# SST
raster::writeRaster(sst_str_mean, "./data/sensitivity/sst/sst_str_mean.asc", overwrite = TRUE)
# raster::writeRaster(sst_str_iqr, "./data/sensitivity/sst/sst_str_iqr.asc", overwrite = TRUE)
raster::writeRaster(sst_tb_mean, "./data/sensitivity/sst/sst_tb_mean.asc", overwrite = TRUE)
raster::writeRaster(sst_tb_iqr, "./data/sensitivity/sst/sst_tb_iqr.asc", overwrite = TRUE)
# SBT
raster::writeRaster(sbt_str_mean, "./data/sensitivity/sbt/sbt_str_mean.asc", overwrite = TRUE)
# raster::writeRaster(sbt_str_iqr, "./data/sensitivity/sbt/sbt_str_iqr.asc", overwrite = TRUE)
raster::writeRaster(sbt_tb_mean, "./data/sensitivity/sbt/sbt_tb_mean.asc", overwrite = TRUE)
raster::writeRaster(sbt_tb_iqr, "./data/sensitivity/sbt/sbt_tb_iqr.asc", overwrite = TRUE)


#### End of code. 
##############################
##############################