##############################
##############################
#### process_eez.R
#### Edward Lavender (el72@st-andrews.ac.uk)

#### This code:
# 1) Gets EEZ areas and latitudinal mid-points 

#### Steps preceding this code: 
# 1) Definition of global extent/modelling resolution 
# ... e.g., extracted from processed species richness data from process_spptraits.R


##############################
##############################
#### Set up 

#### Wipe workspace
rm(list = ls())

#### Essential packages
library(magrittr)

#### Load data
# Load spp richness for species numbers and area dimensions
spp_richness <- raster::raster("./data/spatial/species/map_spp_richness.asc")
# Load EEZ data 
eez <- sf::read_sf("./data/spatial/eez", "eez_v11") 
eez <- as(eez, "Spatial")


##############################
##############################
#### Implement processing 

#### Obtain EEZ in each cell
# This approach takes ~ 5 minutes. 
# A more accurate approach would loop over each EEZ, 
# ... b/c some cells may contain more than one EEZ. 
# ... But this is much slower. 
# ... The approach below is sufficient here and is consistent with the 
# ... approach used to summarise abundance predictions across EEZs. 
eez_in_cell <- spp_richness
eez_in_cell <- data.frame(raster::rasterToPoints(eez_in_cell))
colnames(eez_in_cell) <- c("x", "y", "z")
eez_in_cell$n_cell <- 1
sp::coordinates(eez_in_cell) <- c("x", "y")
raster::crs(eez_in_cell) <- raster::crs(eez)
eez_in_cell$eez <- sp::over(eez_in_cell, eez)$SOVEREIGN1
eez_in_cell <- data.frame(eez_in_cell)
eez_in_cell$optional <- NULL
head(eez_in_cell)

#### Calculate the number of species in each grid cell 
# Convert species richness to dataframe 
dat_spp_richness <- data.frame(raster::rasterToPoints(spp_richness))
colnames(dat_spp_richness) <- c("x", "y", "n_spp")
# Define keys for matching 
eez_in_cell$key      <- paste0("(", eez_in_cell$x, ",", eez_in_cell$y, ")")
dat_spp_richness$key <- paste0("(", dat_spp_richness$x, ",", dat_spp_richness$y, ")")
# Add the number of species in each cell
eez_in_cell$n_spp <- dat_spp_richness$n_spp[match(eez_in_cell$key, dat_spp_richness$key)]
eez_in_cell$n_spp[is.na(eez_in_cell$n_spp)] <- 0
range(eez_in_cell$n_spp)

#### Summarise the statistics for each EEZ
# ... The total number of cells
# ... ... Note that b/c some cells may contain multiple EEZs, this estimate is approximate. 
# ... The total number of cells with at least one species 
# ... ... with a Pr(presence > 0.5). 
# ... analyse_abund_across_eezs.R adds the total number of species in each EEZ too. 
eez_stats <- 
  eez_in_cell %>% 
  dplyr::filter(!is.na(eez)) %>% 
  dplyr::group_by(eez) %>% 
  dplyr::summarise(n_cell_tot = sum(n_cell), 
                   n_cell_occ = sum(n_spp > 0.5), 
                   pc_cell_occ = n_cell_occ/n_cell_tot * 100)

#### Get latitudinal midpoints for each EEZ
# Here, for speed, we get latitudinal mid-points using byid = TRUE. 
# This returns multiple values for each EEZ, which we then average. 
# A more accurate approach would be to loop over each EEZ and, within that EEZ, 
# ... use byid = TRUE, but that takes hours and this faster approach is sufficient here. 
lat_mid <- rgeos::gCentroid(eez, byid = TRUE)
lat_mid <- data.frame(eez = eez$SOVEREIGN1, mid_point = sp::coordinates(lat_mid)[, 2])
lat_mid <- 
  lat_mid %>% dplyr::group_by(eez) %>% 
  dplyr::summarise(mid_point = mean(mid_point))
eez_stats$lat_mid <- lat_mid$mid_point[match(eez_stats$eez, lat_mid$eez)]
eez_stats$lat_mid_abs <- abs(eez_stats$lat_mid)

#### Associate latitudinal midpoints with colours
n_cols <- 1000
pal <- viridis::plasma
cols <- rev(pal(n_cols))
breaks <- seq(0, 90, length.out = n_cols)
eez_stats$lat_band <- cut(abs(eez_stats$lat_mid), breaks = breaks)
eez_stats$col <- cols[unclass(eez_stats$lat_band)]

#### Dataframe finalisation 
eez_in_cell$key    <- NULL
eez_in_cell$n_cell <- NULL

#### Make a mask for temperature projections
# ... We will focus on areas within EEZs
# ... We will make a raster which defines, for each cell, whether or not it is inside an EEZ
# ... and then use this to mask temperature/abundance projections to focus on coastal areas. 
# ... This is much quicker than trying to mask files
# ... directly using the EEZ data.
cover <- eez_in_cell
cover$z <- ifelse(!is.na(cover$eez), 0, NA)
cover$eez   <- NULL
cover$n_spp <- NULL
cover <- raster::rasterFromXYZ(cover, res = 0.5)
cover <- raster::extend(cover, raster::extent(spp_richness))
raster::plot(cover)

#### Define a vector of the EEZs with the largest marine capture production (FAO, 2018)
eez_25 <- c("Malaysia", 
            "Indonesia", 
            "Thailand",
            "Philippines",
            "India",
            "Myanmar", 
            "Vietnam",
            "United Kingdom",
            "United States",
            "Denmark",
            "Morocco",
            "Spain",
            "Russia",
            "Iceland",
            "Argentina",
            "Norway",
            "Canada",
            "Japan",
            "Mexico",
            "South Korea",
            "Chile",
            "Ecuador",
            "China",
            "Taiwan",
            "Peru")

#### Save results
save <- TRUE
if(save){
  saveRDS(eez_in_cell, "./data/spatial/eez/eez_in_cell.rds")
  saveRDS(eez_stats, "./data/spatial/eez/eez_stats.rds")
  saveRDS(eez_25, "./data/spatial/eez/eez_with_largest_marine_capture_production.rds")
  raster::writeRaster(cover, "./data/spatial/eez/eez_mask.asc", overwrite = TRUE)
}


#### End of code 
##############################
##############################