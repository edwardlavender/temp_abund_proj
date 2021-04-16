##############################
##############################
#### process_sdm_aquamaps.R
#### Edward Lavender (el72@st-andrews.ac.uk)

#### This code:
# 1) Processes aquamaps SDMs for selected species

#### Steps preceding this code: 
# 1) process_spptraits.R defines an 'intermediate' list of modelled species
# 2) ... For which get_sdm_aquamaps.R gets 'raw' species distributions 
#    ... At which point, this code should be implemented to process species distributions. 
# 3) Manual acquisition and processing of coastline data used to mask predictions (process_coastline.R)


##############################
##############################
#### Set up 

#### Wipe workspace
rm(list = ls())

#### Essential packages
library(magrittr)

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
# If necessary, re-sample to match the resolution of the temperature files
# ... Instead, temperature files have been re-sampled to match the resolution of species distributions. 
# Ensure 0 'occupancy' to NA
# ... This simplifies acquisition of thermal niche parameters
# ... and temperature projections, by ensuring that we only ever use
# ... areas where a species is predicted to be present, and other areas are ignored as NAs. 
# Remove any cells that overlap with the land via raster::mask()
# ... With this method, some 'cells' will overlap with the coastline and land
# ... but these cells this is because of the low resolution of climate models 
# Expand resolution to standard global extent
# ... To match climate outputs. 
# Define threshold for indicating 'presence' e.g., Pr(presence) >= 0.5 as in Klein et al. (2013). 
# ... Instead, we will base our analyses on the probabilities. 
# Plot distribution
# ... For visual checking. 
# Save processed SDM map. 

#### Loop over each species and save processed predictions 
# [This code takes ~ 3 minutes with 11 cores if figures are generated
# ... or ~ 9 minutes if the figures are generated as well.]
spptraits$index <- 1:nrow(spptraits)
blank <- raster::raster(ext = raster::extent(-180, 180, -90, 90), 
                        resolution = 0.5)
cl <- parallel::makeCluster(11L)
parallel::clusterExport(cl = cl, varlist = c("blank", "coastline"))
pbapply::pblapply(split(spptraits, spptraits$index), cl = cl, function(d){
  ## Process distribution 
  # d <- spptraits[1, ]
  r <- raster::raster(paste0("./data-raw/sdm_aquamaps/maps/", d$spp_key_asc))
  # raster::plot(r)
  r[r == 0] <- NA 
  r <- raster::mask(r, coastline, updatevalue = NA)
  r <- raster::resample(r, blank, method = "ngb")
  ## Plot distribution 
  plot <- FALSE
  if(plot){
    tiff(paste0("./fig/sdm_aquamaps/", d$spp_key, ".tiff"), 
         height = 5, width = 5, units = "in", res = 300)
    raster::plot(r, main = d$spp)
    raster::lines(coastline, lwd = 0.25, col = "dimgrey")
    dev.off()
  }
  ## Save asc file 
  raster::writeRaster(r, paste0("data/sdm_aquamaps/", d$spp_key_asc), overwrite = TRUE)
  return(invisible())
}) %>% invisible()
parallel::stopCluster(cl)


##############################
##############################
#### Checks 

#### Manual comparisons of 'raw' and 'processed' maps 
# E.g. 1
raster::plot(raster::raster(paste0("data-raw/sdm_aquamaps/maps/", spptraits$spp_key_asc[1])))
raster::plot(raster::raster(paste0("data/sdm_aquamaps/", spptraits$spp_key_asc[1])))
# E.g. 2
raster::plot(raster::raster(paste0("data-raw/sdm_aquamaps/maps/", spptraits$spp_key_asc[10])))
raster::plot(raster::raster(paste0("data/sdm_aquamaps/", spptraits$spp_key_asc[10])))
# E.g. 3
raster::plot(raster::raster(paste0("data-raw/sdm_aquamaps/maps/", spptraits$spp_key_asc[1000])))
raster::plot(raster::raster(paste0("data/sdm_aquamaps/", spptraits$spp_key_asc[1000])))

#### Compare 'local' and 'online' versions of maps for a subset of species
# ... to check these match up perfectly. 
## Get fishbase numbers
fb <- rfishbase::species(spptraits$spp)
spptraits$spec_code <- fb$SpecCode
## Define helper function to compare local and online maps 
compare_local_and_online_maps <- function(data){
  url_fishbase <- "https://www.fishbase.se/summary/"
  for(i in 1:nrow(data)){
    spptraits_for_spp <- data[i, , drop = FALSE]
    msg <- paste("[", i, "]", spptraits_for_spp$spp)
    print(msg)
    r <- raster::raster(paste0("./data-raw/sdm_aquamaps/maps/", spptraits_for_spp$spp_key_asc))
    raster::plot(r, main = msg)
    raster::lines(coastline)
    browseURL(paste0(url_fishbase, spptraits_for_spp$spec_code))
    readline(prompt = "Press [enter] to continue...")
  }
}
## Compare maps for a subset of species
compare_local_and_online_maps(spptraits[sample(1:nrow(spptraits), 25), ])

#### Acquire presence data via robis (for checking maps)
# 

#### Make maps, with presence data overlaid (for checking maps)
#

#### Qualitative examination of individual maps
# ... look at species distributions
# ... look for 'artefacts' in the distributions (e.g., induced by FAO boundaries) 



#### End of code.
##############################
##############################