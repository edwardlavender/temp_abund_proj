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
# spptraits_pro <- readRDS("./data/spptraits.rds")
# spptraits     <- spptraits[spptraits$spp %in% spptraits_pro$spp, ]
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
# ... We could also expand the extent to the smallest possible area that includes all SDMs
# ... ... for computational speed later, but this approach is simpler. 
# Define threshold for indicating 'presence' e.g., Pr(presence) >= 0.5 as in Klein et al. (2013). 
# ... This is computationally beneficial when it comes to synthesising predictions across species. 
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
  # Make threshold-based map 
  r_occ <- r
  r_occ[r_occ >= 0.5] <- 1
  r_occ[r_occ < 0.5] <- NA
  # raster::plot(r)
  # raster::plot(r_occ)
  ## Plot distribution 
  plot <- FALSE
  if(plot){
    tiff(paste0("./fig/sdm_aquamaps/", d$spp_key, ".tiff"), 
         height = 5, width = 5, units = "in", res = 300)
    raster::plot(r, main = d$spp)
    raster::lines(coastline, lwd = 0.25, col = "dimgrey")
    dev.off()
  }
  ## Save asc files
  raster::writeRaster(r, paste0("data/sdm_aquamaps/maps_pr/", d$spp_key_asc), overwrite = TRUE)
  raster::writeRaster(r_occ, paste0("data/sdm_aquamaps/maps_occ/", d$spp_key_asc), overwrite = TRUE)
  return(invisible())
}) %>% invisible()
parallel::stopCluster(cl)


##############################
##############################
#### Checks 

#### Manual comparisons of 'raw' and 'processed' maps 
# E.g. 1
raster::plot(raster::raster(paste0("data-raw/sdm_aquamaps/maps/", spptraits$spp_key_asc[1])))
raster::plot(raster::raster(paste0("data/sdm_aquamaps/maps_pr/", spptraits$spp_key_asc[1])))
raster::plot(raster::raster(paste0("data/sdm_aquamaps/maps_occ/", spptraits$spp_key_asc[1])))
# E.g. 2
raster::plot(raster::raster(paste0("data-raw/sdm_aquamaps/maps/", spptraits$spp_key_asc[10])))
raster::plot(raster::raster(paste0("data/sdm_aquamaps/maps_occ/", spptraits$spp_key_asc[10])))
# E.g. 3
raster::plot(raster::raster(paste0("data-raw/sdm_aquamaps/maps/", spptraits$spp_key_asc[1000])))
raster::plot(raster::raster(paste0("data/sdm_aquamaps/maps_occ/", spptraits$spp_key_asc[1000])))

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


##############################
#############################
#### Acquire presence data via robis and rgbif (for checking maps)

get_occurences <- FALSE
if(get_occurences){
  #### Acquire gbif data for each species
  # (Use a loop to mitigate poor internet)
  gbif_by_spp <- list()
  lx <- nrow(spptraits)
  for(i in 1:lx){
    print(i)
    # svMisc::progress(i, lx)
    out <- tryCatch(rgbif::occ_data(scientificName = spptraits$spp[i], limit = 1000), 
                    error = function(e) e)
    gbif_by_spp[[i]] <- out
  }
  # Check for errors
  # There is one species for which we had an error
  # We will download the data for this species again: 
  table(unlist(sapply(gbif_by_spp, class)))
  pos_try_again <- which(sapply(gbif_by_spp, function(elm) inherits(elm, "error")))
  gbif_by_spp[[pos_try_again]] 
  gbif_by_spp[[pos_try_again]] <- rgbif::occ_data(scientificName = spptraits$spp[pos_try_again], 
                                                  limit = 1000)
  # Save raw results 
  # saveRDS(gbif_by_spp, "./data-raw/occurrence/gbif_by_spp.rds")
  
  #### Acquire OBIS data for each species
  # (Use a loop to mitigate poor internet)
  obis_by_spp <- list()
  for(i in 1:lx){
    # svMisc::progress(i, lx)
    print(i)
    out <- tryCatch(robis::occurrence(scientificname = spptraits$spp[i]), 
                    error = function(e) e)
    obis_by_spp[[i]] <- out
  }
  # Check for errors: we have got data for all species. 
  table(unlist(sapply(obis_by_spp, class)))
  # Save raw results 
  # saveRDS(obis_by_spp, "./data-raw/occurrence/obis_by_spp.rds")
}

#### Process GBIF data
gbif_by_spp_xy <- pbapply::pblapply(gbif_by_spp, function(gbif_for_spp){
  # gbif_for_spp <- gbif_by_spp[[1]]
  gbif_xy <- NULL
  if(length(gbif_for_spp) > 0){
    if(!is.null(gbif_for_spp$data)){
      if("decimalLongitude" %in% names(gbif_for_spp$data)){
        gbif_xy <- gbif_for_spp$data[, c("species", "decimalLongitude", "decimalLatitude")]
        gbif_xy$source <- "gbif"
        colnames(gbif_xy) <- c("spp", "long", "lat", "source")
        gbif_xy <- gbif_xy[gbif_xy$long != 0 & gbif_xy$lat != 0, ]
      }
    }
  }
  return(gbif_xy)
})
# No records for Centropogon australis
gbif_by_spp[[which(sapply(gbif_by_spp_xy, is.null))]]
gbif_by_spp_xy <- plyr::compact(gbif_by_spp_xy)
names(gbif_by_spp_xy) <- sapply(gbif_by_spp_xy, function(elm) elm$spp[1])
# saveRDS(gbif_by_spp_xy, "./data/occurrence/gbif_by_spp_xy.rds")   

#### Process OBIS data 
obis_by_spp_xy <- pbapply::pblapply(obis_by_spp, function(obis_for_spp){
  # obis_for_spp <- obis_by_spp[[1]]
  obis_xy <- NULL
  if(length(obis_for_spp) > 0){
    if("decimalLongitude" %in% names(obis_for_spp)){
      obis_xy <- obis_for_spp[, c("species", "decimalLongitude", "decimalLatitude")]
      obis_xy$source <- "obis"
      colnames(obis_xy) <- c("spp", "long", "lat", "source")
      obis_xy <- obis_xy[obis_xy$long != 0 & obis_xy$lat != 0, ]
    }
  }
  return(obis_xy)
})
# There are records for all species:
which(sapply(obis_by_spp_xy, is.null))
names(obis_by_spp_xy) <- sapply(obis_by_spp_xy, function(elm) elm$spp[1])
# saveRDS(obis_by_spp_xy, "./data/occurrence/obis_by_spp_xy.rds")   

#### Join GBIF and OBIS datasets for convenience
gbif_xy <- dplyr::bind_rows(gbif_by_spp_xy)
obis_xy <- dplyr::bind_rows(obis_by_spp_xy)
occurence_xy <- rbind(gbif_xy, obis_xy)
occurence_xy$spp <- factor(occurence_xy$spp, levels = unique(spptraits$spp))
occurence_xy$source <- factor(occurence_xy$source)
occurence_xy <- occurence_xy %>% dplyr::arrange(spp, source)
# saveRDS(occurence_xy, "./data/occurrence/occurence_xy.rds")


##############################
#############################
#### Make maps, with presence data overlaid (for checking maps)

#### Convenience function to make maps 
map_sdm <- function(r, coastline, xy_gbif, xy_obis,...){
  prettyGraphics::pretty_map(r, 
                             add_rasters = list(x = r, 
                                                zlim = c(0.95, 1.05), 
                                                col = "black",
                                                plot_method = raster::image),...
                             )
  raster::lines(coastline, col = "darkgreen", lwd = 2)
  if(nrow(xy_gbif) > 0) points(xy_gbif$long, xy_gbif$lat, 
                               pch = 4, col = scales::alpha("royalblue", 0.5), lwd = 1, cex = 0.5)
  if(nrow(xy_obis) > 0) points(xy_obis$long, xy_obis$lat, 
                               pch = 4, col = scales::alpha("darkred", 0.5), lwd = 1, cex = 0.5)
}

#### Make maps 
t1_maps <- Sys.time()
cl <- parallel::makeCluster(8L)
parallel::clusterExport(cl, varlist = c("occurence_xy", "map_sdm"))
pbapply::pblapply(split(spptraits, 1:nrow(spptraits)), function(d){
  # Load rasters 
  # d <- spptraits[1, ]
  r_raw <- raster::raster(paste0("./data-raw/sdm_aquamaps/maps/", d$spp_key_asc))
  r_occ <- raster::raster(paste0("./data/sdm_aquamaps/maps_occ/", d$spp_key_asc))
  r_pr  <- raster::raster(paste0("./data/sdm_aquamaps/maps_occ/", d$spp_key_asc))
  # Obtain occurrence data 
  xy <- occurence_xy %>% dplyr::filter(spp == d$spp)
  xy_gbif <- xy %>% dplyr::filter(source == "gbif")
  xy_obis <- xy %>% dplyr::filter(source == "obis")
  # Crop spatial data for ease of visualisation 
  ext   <- raster::extent(r_raw)
  xlim  <- range(c(xy$long, ext[1:2]))
  ylim  <- range(c(xy$lat, ext[3:4]))
  ext   <- raster::extent(xlim, ylim)
  r_occ <- raster::crop(r_occ, ext)
  r_pr  <- raster::crop(r_pr, ext)
  coast <- raster::crop(coastline, ext)
  # Map rasters
  tiff(paste0("./fig/sdm_aquamaps/", d$spp_key, ".tiff"), 
       height = 4, width = 12, units = "in", res = 300)
  pp <- par(mfrow = c(1, 2), oma = c(2, 2, 2, 2), mar = c(1, 1, 1, 1))
  map_sdm(r_occ, coast, xy_gbif, xy_obis, main = paste(d$spp, "[occ]"))
  map_sdm(r_pr, coast, xy_gbif, xy_obis, main = paste(d$spp, "[pr]"))
  dev.off()
  return(invisible())
}) %>% invisible()
parallel::stopCluster(cl)
t2_maps <- Sys.time()
difftime(t2_maps, t1_maps)

#### Qualitative examination of individual maps
# ... look at species distributions
# ... look for 'artefacts' in the distributions (e.g., induced by FAO boundaries) 



#### End of code.
##############################
##############################