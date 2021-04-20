##############################
##############################
#### project_abund_2.R 
#### Edward Lavender (el72@st-andrews.ac.uk)

#### This code:
# 1) This code synthesises abundance projections, defining, for each climate scenario
# ... the mean change in abundance in each grid cell. 

#### Code duration: 
# Synthesising the predictions for each scenario takes about 20 minutes (10 for the mean, 10 for the Pr)
# So the total computation time is 2 (SST, SBT) x 4 (mid, late, RCP 4.5, 8.5) x 20 mins
# This code is not currently set up to run in parallel, but can be opened in multiple RStudio
# ... windows for parallel computations very simply. 

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
spp_richness <- raster::raster("./data/spatial/species/map_spp_richness.asc")

#### Global param 
root_ab_delta_sst <- "./data/abundance/change/spp_specific/sst/"
root_ab_delta_sbt <- "./data/abundance/change/spp_specific/sbt/"
source("./R/helpers.R")


##############################
##############################
#### Helpers

#### Load files (in parallel)
load_projections <- function(root = paste0(root_ab_delta_sst, "mid_century/rcp45/"), 
                             cl = parallel::makeCluster(2L), 
                             varlist = "root_ab_delta_sst"){
  parallel::clusterExport(cl = cl, varlist = varlist)
  rls <- pbapply::pblapply(split(spptraits, spptraits$index), cl = cl, function(d){
    raster::raster(paste0(root, d$spp_key_asc))
  })
  parallel::stopCluster(cl)
  return(rls)
}

#### Summarise projections 
process_projections <- function(x, type = "sst", file_name = deparse(substitute(x)), 
                                make = c("mean", "pr"), denom){
  
  ## Global mean 
  if("mean" %in% make){
    cat("Summarising rasters (mean)...\n")
    x_sum <- calc_sum(x) 
    x_mean <- x_sum/denom
    raster::plot(x_mean, main = paste(file_name, "[mean]"))
    raster::writeRaster(x_mean, paste0("./data/abundance/change/globe/", type, "/mean/", file_name, "_mean.asc"), 
                        overwrite = TRUE)
  }

  ## Proportions 
  if("pr" %in% make){
    cat("Summarising rasters (Pr)...\n")
    # Define list to summarise 
    xp <- pbapply::pblapply(x, function(r){
      rp <- r < 0
      return(rp)
    })
    # Summarise list 
    xp_sum <- calc_sum(xp) 
    xp_pr <- xp_sum/denom
    raster::plot(xp_pr, main = paste(file_name, "[pr]"))
    raster::writeRaster(xp_pr, paste0("./data/abundance/change/globe/", type, "/pr/", file_name, "_pr.asc"), 
                        overwrite = TRUE)
    return(invisible())
  }
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
process_projections(ab_delta_sst_mid_rcp45, denom = spp_richness)
## RCP 8.5 
ab_delta_sst_mid_rcp85 <- load_projections(root = paste0(root_ab_delta_sst, "mid_century/rcp85/"), 
                                           varlist = "root_ab_delta_sst")
process_projections(ab_delta_sst_mid_rcp85, denom = spp_richness)

#### Late century, RCP 4.5 
## RCP 4.5 
ab_delta_sst_late_rcp45 <- load_projections(root = paste0(root_ab_delta_sst, "late_century/rcp45/"), 
                                            varlist = "root_ab_delta_sst")
process_projections(ab_delta_sst_late_rcp45, denom = spp_richness)
## RCP 8.5 
ab_delta_sst_late_rcp85 <- load_projections(root = paste0(root_ab_delta_sst, "late_century/rcp85/"), 
                                            varlist = "root_ab_delta_sst")
process_projections(ab_delta_sst_late_rcp85, denom = spp_richness)


##############################
#### SBT 

#### Mid_century
## RCP 4.5 
ab_delta_sbt_mid_rcp45 <- load_projections(root = paste0(root_ab_delta_sbt, "mid_century/rcp45/"), 
                                           varlist = "root_ab_delta_sbt")
process_projections(ab_delta_sbt_mid_rcp45, type = "sbt", denom = spp_richness)
## RCP 8.5 
ab_delta_sbt_mid_rcp85 <- load_projections(root = paste0(root_ab_delta_sbt, "mid_century/rcp85/"), 
                                           varlist = "root_ab_delta_sbt")
process_projections(ab_delta_sbt_mid_rcp85, type = "sbt", denom = spp_richness)

#### Late century, RCP 4.5 
## RCP 4.5 
ab_delta_sbt_late_rcp45 <- load_projections(root = paste0(root_ab_delta_sbt, "late_century/rcp45/"), 
                                            varlist = "root_ab_delta_sbt")
process_projections(ab_delta_sbt_late_rcp45, type = "sbt", denom = spp_richness)
## RCP 8.5 
ab_delta_sbt_late_rcp85 <- load_projections(root = paste0(root_ab_delta_sbt, "late_century/rcp85/"), 
                                            varlist = "root_ab_delta_sbt")
process_projections(ab_delta_sbt_late_rcp85, type = "sbt", denom = spp_richness)


#### End of code. 
##############################
##############################