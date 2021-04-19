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
source("./R/helpers.R")

#### Load data
spptraits <- readRDS("./data/spptraits.rds")
spp_bk   <- raster::brick("./data/spatial/species/spp_bk.tif")
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

make_bricks <- FALSE
if(make_bricks){
  
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
      sdm <- raster::raster(paste0("./data/sdm_aquamaps/", d$spp_key_asc))
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
    spp_bk <- raster::subset(spp_bk, 1:10)
  }
  
  #### Brick rasters for each metric 
  # This takes ~ 10 minutes. 
  t1_bk <- Sys.time()
  bk_by_metric <- list(sst_str_by_spp, 
                       sbt_str_by_spp, 
                       sst_tb_by_spp, 
                       sbt_tb_by_spp)
  names(bk_by_metric) <- c("sst_str_bk", "sbt_str_bk", "sst_tb_bk", "sbt_tb_bk")
  cl <- parallel::makeCluster(length(bk_by_metric))
  parallel::clusterExport(cl, varlist = "bk_by_metric")
  pbapply::pblapply(1:length(bk_by_metric), cl = cl, function(i){
    bk <- raster::brick(bk_by_metric[[i]])
    temp_dir <- paste0(substr(names(bk_by_metric)[i], 1, 3), "/")
    raster::writeRaster(bk, paste0("./data/sensitivity/", temp_dir, names(bk_by_metric)[i], ".tif"), overwrite = TRUE)
    return(invisible())
  })
  parallel::stopCluster(cl)
  t2_bk <- Sys.time()
}


##############################
##############################
#### Summarise bricks 

#### Load bricks
sst_str_bk <- raster::brick("./data/sensitivity/sst/sst_str_bk.tif")
sbt_str_bk <- raster::brick("./data/sensitivity/sbt/sbt_str_bk.tif")
sst_tb_bk <- raster::brick("./data/sensitivity/sst/sst_tb_bk.tif")
sbt_tb_bk <- raster::brick("./data/sensitivity/sbt/sbt_tb_bk.tif")

#### Tailor memory options 
# https://strimas.com/post/processing-large-rasters-in-r/ 
raster::canProcessInMemory(sst_str_bk, verbose = TRUE)
raster::rasterOptions(maxmemory = 8*1e9, memfrac = 0.8)

#### Weighted averages 
## sst_str
t1_wt_sst_str <- Sys.time()
sst_str_mean <- raster::weighted.mean(sst_str_bk, spp_bk, na.rm = TRUE)
raster::writeRaster(sst_str_mean, paste0("./data/sensitivity/sst/sst_str_mean.asc"), overwrite = TRUE)
t2_wt_sst_str <- Sys.time()
beepr::beep(10)
difftime(t2_wt_sst_str, t1_wt_sst_str) # Time difference of 10.44203 hours

#### Weighted variability for thermal bias 
# Get thermal bias bricks 
sst_tb_bk <- bk_by_metric$sst_tb
sbt_tb_bk <- bk_by_metric$sbt_tb
# Get weighted IQR 
sst_tb_iqr <- rb_weighted_quantiles(sst_tb_bk, spp_bk, cl = parallel::makeCluster(2L))
sst_tb_iqr <- sst_tb_iqr[[2]] - sst_tb_iqr[[1]]
raster::writeRaster(sst_tb_iqr, "./data/sensitivity/sst/sst_tb_iqr.asc")
sbt_tb_iqr <- rb_weighted_quantiles(sbt_tb_bk, spp_bk)
sbt_tb_iqr <- sbt_tb_iqr[[2]] - sbt_tb_iqr[[1]]
raster::writeRaster(sbt_tb_iqr, "./data/sensitivity/sbt/sbt_tb_iqr.asc")


#### End of code. 
##############################
##############################