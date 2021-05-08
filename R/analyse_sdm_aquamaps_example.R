##############################
##############################
#### analyse_sdm_aquamaps_example.R
#### Edward Lavender (el72@st-andrews.ac.uk)

#### This code:
# 1) Analyses an example SDM (for the GitHub README figure). 

#### Steps preceding this code: 
# 1) Define spptraits via process_spptraits.R
# 2) Process SDMs via process_aquamaps_sdm.R


##############################
##############################
#### Set up 

#### Wipe workspace
rm(list = ls())

#### Define helper function to zoom into a plot easily 
crop_from_click <- function(x, plot = TRUE,...){
  cat("Please click four boundary locations on the map and press [Esc] when you are done...")
  dat <- locator()
  dat <- data.frame(x = dat$x, y = dat$y)
  xlim <- range(dat$x)
  ylim <- range(dat$y)
  ext <- raster::extent(xlim[1], xlim[2], ylim[1], ylim[2])
  x_crop <- raster::crop(x, ext)
  if(plot){
    prettyGraphics::pretty_map(x_crop, 
                               add_rasters = list(x = x_crop),...)
  }
  return(x_crop)
}

#### Load data (example SDM)
spptaits <- readRDS("./data/spptraits.rds")
d <- spptraits[spptraits$spp %in% "Acipenser oxyrinchus", ]
r <- raster::raster(paste0("./data/sdm_aquamaps/maps_occ/", d$spp_key_asc))


##############################
##############################
#### Plot example SDM 

#### Quick plot 
raster::plot(r, col = "darkgreen")
raster::lines(coastline)

#### Define zoomed in SDM interactively
r2 <- crop_from_click(r)

#### Plot zoomed in SDM 
raster::plot(r2, col = "darkred")
raster::plot(coastline, col = scales::alpha("skyblue", 0.5), add = TRUE)
land <- flapper::invert_poly(coastline)
raster::plot(r2, col = "darkred", add = TRUE)
raster::plot(land, lwd = 2, col = "white", add = TRUE)
raster::plot(land, lwd = 4, col = scales::alpha("lightgreen", 0.5), add = TRUE)
text(-90, 40, "Land", cex = 5, font = 2)
text(-60, 30, "Sea", cex = 5, font = 2)


#### End of code. 
##############################
##############################