##############################
##############################
#### project_abund_2.R 
#### Edward Lavender (el72@st-andrews.ac.uk)

#### This code:
# 1) This code synthesises abundance projections, defining, for each climate scenario
# ... a brick of the outcomes across species, and summary statistics. For each scenario, 
# ... this includes:
# ... ... mean change in abundance in each grid cell
# ... ... sd change in abundance in each grid cell 

#### Code duration: 
# This code makes 2 x 4 raster bricks (SST, SBT, for for climate scenarios)
# ... Each raster brick takes ~ 2 - 3 hours to make, so this code takes a full 24 hours to run. 

#### Steps preceding this code:
# 1) Make abundance projections (project_abund_1.R)


##############################
##############################
#### Set up 

#### Wipe workspace
rm(list = ls())

#### Load data 
spptraits <- readRDS("./data/spptraits.rds")
spptraits$index <- 1:nrow(spptraits)
# spptraits <- spptraits[1:10, ]
# species richness map 
spp_richness <- raster::raster("./data/spatial/map_spp_richness.asc")

#### Global param 
root_ab_delta_sst <- "./data/abundance/change/spp_specific/sst/"
root_ab_delta_sbt <- "./data/abundance/change/spp_specific/sbt/"


##############################
##############################
#### Helpers

#### Load files (in parallel)
load_projections <- function(root = paste0(root_ab_delta_sst, "mid_century/rcp45/"), 
                             cl = parallel::makeCluster(4L), 
                             varlist = "root_ab_delta_sst"){
  parallel::clusterExport(cl = cl, varlist = varlist)
  rls <- pbapply::pblapply(split(spptraits, spptraits$index), cl = cl, function(d){
    raster::raster(paste0(root, d$number.asc.file))
  })
  parallel::stopCluster(cl)
  return(rls)
}

#### Brick files (x) and summarise via mean and SD 
process_projections <- function(x, brick_names = spptraits$number.asc.file, type = "sst", file_name = deparse(substitute(x))){
  # Brick 
  t1 <- Sys.time()
  cat("Bricking rasters...\n")
  x_bk <- raster::brick(x)
  names(x_bk) <- brick_names
  raster::writeRaster(x_bk, paste0("./data/abundance/change/globe/", type, "/bricks/", file_name, "_bk.tif"), 
                      format = "GTiff", options = c("INTERLEAVE=BAND","COMPRESS=LZW"), 
                      overwrite = TRUE)
  t2 <- Sys.time()
  tdiff <- round(as.numeric(difftime(t2, t1, units = "mins")))
  cat("... raster::brick() complete after", tdiff, "minute(s).\n")
  
  ## Global mean 
  t1 <- Sys.time()
  cat("Summarising rasters (mean)...\n")
  x_mean <- raster::calc(x_bk, mean, na.rm = TRUE)
  raster::writeRaster(x_mean, paste0("./data/abundance/change/globe/", type, "/mean/", file_name, "_mean.asc"), 
                      overwrite = TRUE)
  t2 <- Sys.time()
  tdiff <- round(as.numeric(difftime(t2, t1, units = "mins")))
  cat("... raster::calc() complete after", tdiff, "minute(s).\n") 
  
  ## Global SD 
  t1 <- Sys.time()
  cat("Summarising rasters (SD)...\n")
  x_sd <- raster::calc(x_bk, sd, na.rm = TRUE)
  raster::writeRaster(x_sd, paste0("./data/abundance/change/globe/", type, "/sd/", file_name, "_sd.asc"), 
                      overwrite = TRUE)
  t2 <- Sys.time()
  tdiff <- round(as.numeric(difftime(t2, t1, units = "mins")))
  cat("... raster::calc() complete after", tdiff, "minute(s).\n") 
  return(invisible())
}

#### Helper function to calculate the proportion of species expected to decline in each cell
# This could be implemented within process_projections() but it is implemented afterwards for simplicity. 
process_proportions <- function(x, denom, plot = TRUE){
  x_pr <- raster::calc(x, function(x) x < 0)
  x_pr <- sum(x_pr, na.rm = TRUE)/denom
  if(plot) raster::plot(x_pr)
  return(x_pr)
}


##############################
##############################
#### Process rasters 

##############################
#### SST 

#### Mid_century
## RCP 4.5 
ab_delta_sst_mid_rcp45 <- load_projections(root = paste0(root_ab_delta_sst, "mid_century/rcp45/"),
                                           varlist = "root_ab_delta_sst")
process_projections(ab_delta_sst_mid_rcp45)
## RCP 8.5 
ab_delta_sst_mid_rcp85 <- load_projections(root = paste0(root_ab_delta_sst, "mid_century/rcp85/"), 
                                           varlist = "root_ab_delta_sst")
process_projections(ab_delta_sst_mid_rcp85)

#### Late century, RCP 4.5 
## RCP 4.5 
ab_delta_sst_late_rcp45 <- load_projections(root = paste0(root_ab_delta_sst, "late_century/rcp45/"), 
                                            varlist = "root_ab_delta_sst")
process_projections(ab_delta_sst_late_rcp45)
## RCP 8.5 
ab_delta_sst_late_rcp85 <- load_projections(root = paste0(root_ab_delta_sst, "late_century/rcp85/"), 
                                            varlist = "root_ab_delta_sst")
process_projections(ab_delta_sst_late_rcp85)

#### Additional statistics 
## Load bricks 
t1_pr <- Sys.time()
ab_delta_sst_mid_rcp45_bk  <- raster::brick("./data/abundance/change/globe/sst/bricks/ab_delta_sst_mid_rcp45_bk.tif")
ab_delta_sst_late_rcp45_bk <- raster::brick("./data/abundance/change/globe/sst/bricks/ab_delta_sst_late_rcp45_bk.tif")
ab_delta_sst_mid_rcp85_bk  <- raster::brick("./data/abundance/change/globe/sst/bricks/ab_delta_sst_mid_rcp85_bk.tif")
ab_delta_sst_late_rcp85_bk <- raster::brick("./data/abundance/change/globe/sst/bricks/ab_delta_sst_late_rcp85_bk.tif")
## IQR
ab_delta_sst_mid_rcp45_iqr  <- raster::calc(ab_delta_sst_mid_rcp45_bk, IQR, na.rm = TRUE)
ab_delta_sst_late_rcp45_iqr <- raster::calc(ab_delta_sst_late_rcp45_bk, IQR, na.rm = TRUE)
ab_delta_sst_mid_rcp85_iqr  <- raster::calc(ab_delta_sst_mid_rcp85_bk, IQR, na.rm = TRUE)
ab_delta_sst_late_rcp85_iqr <- raster::calc(ab_delta_sst_late_rcp85_bk, IQR, na.rm = TRUE)
## Calculate proportions across each brick 
ab_delta_sst_mid_rcp45_pr  <- process_proportions(ab_delta_sst_mid_rcp45_bk, spp_richness)
ab_delta_sst_late_rcp45_pr <- process_proportions(ab_delta_sst_late_rcp45_bk, spp_richness)
ab_delta_sst_mid_rcp85_pr  <- process_proportions(ab_delta_sst_mid_rcp85_bk, spp_richness)
ab_delta_sst_late_rcp85_pr <- process_proportions(ab_delta_sst_late_rcp85_bk, spp_richness)
## Save rasters 
# IQR
raster::writeRaster(ab_delta_sst_mid_rcp45_iqr, "./data/abundance/change/globe/sst/iqr/ab_delta_sst_mid_rcp45_iqr.asc")
raster::writeRaster(ab_delta_sst_late_rcp45_iqr, "./data/abundance/change/globe/sst/iqr/ab_delta_sst_late_rcp45_iqr.asc")
raster::writeRaster(ab_delta_sst_mid_rcp85_iqr, "./data/abundance/change/globe/sst/iqr/ab_delta_sst_mid_rcp85_iqr.asc")
raster::writeRaster(ab_delta_sst_late_rcp85_iqr, "./data/abundance/change/globe/sst/iqr/ab_delta_sst_late_rcp85_iqr.asc")
# Proportions 
raster::writeRaster(ab_delta_sst_mid_rcp45_pr, "./data/abundance/change/globe/sst/pr/ab_delta_sst_mid_rcp45_pr.asc")
raster::writeRaster(ab_delta_sst_late_rcp45_pr, "./data/abundance/change/globe/sst/pr/ab_delta_sst_late_rcp45_pr.asc")
raster::writeRaster(ab_delta_sst_mid_rcp85_pr, "./data/abundance/change/globe/sst/pr/ab_delta_sst_mid_rcp85_pr.asc")
raster::writeRaster(ab_delta_sst_late_rcp85_pr, "./data/abundance/change/globe/sst/pr/ab_delta_sst_late_rcp85_pr.asc")
t2_pr <- Sys.time()
difftime(t2_pr, t1_pr)


##############################
#### SBT 

#### Mid_century
## RCP 4.5 
ab_delta_sbt_mid_rcp45 <- load_projections(root = paste0(root_ab_delta_sbt, "mid_century/rcp45/"), 
                                           varlist = "root_ab_delta_sbt")
process_projections(ab_delta_sbt_mid_rcp45, type = "sbt")
## RCP 8.5 
ab_delta_sbt_mid_rcp85 <- load_projections(root = paste0(root_ab_delta_sbt, "mid_century/rcp85/"), 
                                           varlist = "root_ab_delta_sbt")
process_projections(ab_delta_sbt_mid_rcp85, type = "sbt")

#### Late century, RCP 4.5 
## RCP 4.5 
ab_delta_sbt_late_rcp45 <- load_projections(root = paste0(root_ab_delta_sbt, "late_century/rcp45/"), 
                                            varlist = "root_ab_delta_sbt")
process_projections(ab_delta_sbt_late_rcp45, type = "sbt")
## RCP 8.5 
ab_delta_sbt_late_rcp85 <- load_projections(root = paste0(root_ab_delta_sbt, "late_century/rcp85/"), 
                                            varlist = "root_ab_delta_sbt")
process_projections(ab_delta_sbt_late_rcp85, type = "sbt")

#### Additional statistics 
## Load bricks 
t1_pr <- Sys.time()
ab_delta_sbt_mid_rcp45_bk  <- raster::brick("./data/abundance/change/globe/sbt/bricks/ab_delta_sbt_mid_rcp45_bk.tif")
ab_delta_sbt_late_rcp45_bk <- raster::brick("./data/abundance/change/globe/sbt/bricks/ab_delta_sbt_late_rcp45_bk.tif")
ab_delta_sbt_mid_rcp85_bk  <- raster::brick("./data/abundance/change/globe/sbt/bricks/ab_delta_sbt_mid_rcp85_bk.tif")
ab_delta_sbt_late_rcp85_bk <- raster::brick("./data/abundance/change/globe/sbt/bricks/ab_delta_sbt_late_rcp85_bk.tif")
## IQR
ab_delta_sbt_mid_rcp45_iqr  <- raster::calc(ab_delta_sbt_mid_rcp45_bk, IQR, na.rm = TRUE)
ab_delta_sbt_late_rcp45_iqr <- raster::calc(ab_delta_sbt_late_rcp45_bk, IQR, na.rm = TRUE)
ab_delta_sbt_mid_rcp85_iqr  <- raster::calc(ab_delta_sbt_mid_rcp85_bk, IQR, na.rm = TRUE)
ab_delta_sbt_late_rcp85_iqr <- raster::calc(ab_delta_sbt_late_rcp85_bk, IQR, na.rm = TRUE)
## Calculate proportions across each brick 
ab_delta_sbt_mid_rcp45_pr  <- process_proportions(ab_delta_sbt_mid_rcp45_bk, spp_richness)
ab_delta_sbt_late_rcp45_pr <- process_proportions(ab_delta_sbt_late_rcp45_bk, spp_richness)
ab_delta_sbt_mid_rcp85_pr  <- process_proportions(ab_delta_sbt_mid_rcp85_bk, spp_richness)
ab_delta_sbt_late_rcp85_pr <- process_proportions(ab_delta_sbt_late_rcp85_bk, spp_richness)
## Save rasters
# IQR
raster::writeRaster(ab_delta_sbt_mid_rcp45_iqr, "./data/abundance/change/globe/sbt/iqr/ab_delta_sbt_mid_rcp45_iqr.asc")
raster::writeRaster(ab_delta_sbt_late_rcp45_iqr, "./data/abundance/change/globe/sbt/iqr/ab_delta_sbt_late_rcp45_iqr.asc")
raster::writeRaster(ab_delta_sbt_mid_rcp85_iqr, "./data/abundance/change/globe/sbt/iqr/ab_delta_sbt_mid_rcp85_iqr.asc")
raster::writeRaster(ab_delta_sbt_late_rcp85_iqr, "./data/abundance/change/globe/sbt/iqr/ab_delta_sbt_late_rcp85_iqr.asc")
# Proportions 
raster::writeRaster(ab_delta_sbt_mid_rcp45_pr, "./data/abundance/change/globe/sbt/pr/ab_delta_sbt_mid_rcp45_pr.asc")
raster::writeRaster(ab_delta_sbt_late_rcp45_pr, "./data/abundance/change/globe/sbt/pr/ab_delta_sbt_late_rcp45_pr.asc")
raster::writeRaster(ab_delta_sbt_mid_rcp85_pr, "./data/abundance/change/globe/sbt/pr/ab_delta_sbt_mid_rcp85_pr.asc")
raster::writeRaster(ab_delta_sbt_late_rcp85_pr, "./data/abundance/change/globe/sbt/pr/ab_delta_sbt_late_rcp85_pr.asc")
t2_pr <- Sys.time()
difftime(t2_pr, t1_pr)
beepr::beep(10)


#### End of code. 
##############################
##############################