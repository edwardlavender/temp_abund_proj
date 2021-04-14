##############################
##############################
#### analyse_temp.R
#### Edward Lavender (el72@st-andrews.ac.uk)

#### This code:
# 1) Analyses future temperature projections for SST and SBT

#### Steps preceding this code: 
# 1) Process temperature projections via process_temp.R


##############################
##############################
#### Set up 

#### Wipe workspace
rm(list = ls())

#### Essential packages
library(prettyGraphics)

#### Global param
root_sst <- "./data/temperature/sst"
root_sbt <- "./data/temperature/sbt"

#### Load files 
## SST
list.files(root_sst, recursive = TRUE, pattern = ".*asc", full.names = TRUE)
sst_historical <- raster::raster("./data/temperature/sst/historical/historical.asc")
sst_mid_rcp45  <- raster::raster("./data/temperature/sst/mid_century/rcp45/rcp45.asc")
sst_mid_rcp85  <- raster::raster("./data/temperature/sst/mid_century/rcp85/rcp85.asc")
sst_late_rcp45 <- raster::raster("./data/temperature/sst/late_century/rcp45/rcp45.asc")
sst_late_rcp85 <- raster::raster("./data/temperature/sst/late_century/rcp85/rcp85.asc")
## SBT 
list.files(root_sbt, recursive = TRUE, pattern = ".*asc", full.names = TRUE)
sbt_historical <- raster::raster("./data/temperature/sbt/historical/historical.asc")
sbt_mid_rcp45  <- raster::raster("./data/temperature/sbt/mid_century/rcp45/rcp45.asc")
sbt_mid_rcp85  <- raster::raster("./data/temperature/sbt/mid_century/rcp85/rcp85.asc")
sbt_late_rcp45 <- raster::raster("./data/temperature/sbt/late_century/rcp45/rcp45.asc")
sbt_late_rcp85 <- raster::raster("./data/temperature/sbt/late_century/rcp85/rcp85.asc")
## Spatial data
coastline <- readRDS("./data/spatial/coastline/coastline.rds")
land <- flapper::invert_poly(coastline)

#### Helper functions
summarise_raster <- function(x){
  raster::hist(x)
  print(raster::cellStats(x, mean))
  # print(raster::cellStats(x, median))
  print(raster::cellStats(x, range))
  return(invisible())
}

#### Define graphical param 
paa <-  list(side = 1:4, control_axis = list(lwd.ticks = 0, labels = FALSE, lwd = 2))
axis_args <- list(cex.axis = 1.25)

##############################
##############################
#### Historical temperatures

#### Set up plot to save 
tiff("./fig/sst_sbt_historical.tiff", 
     height = 5, width = 6, units = "in", res = 600)
pp <- par(mfrow = c(2, 1), oma = c(0, 0, 0.5, 4), mar = c(0, 0, 0.5, 0.5))

#### Param 
sp <- c(0.93, 0.96, 0.2, 0.8)
line_main <- -0.4
adj_main  <- 0.1
cex_main  <- 1.5
cols <- RColorBrewer::brewer.pal(3, "RdYlBu")
cols <- rev(colorRampPalette(cols)(100))

#### sst_historical
raster::cellStats(sst_historical, range)
pretty_map(sst_historical, 
           pretty_axis_args = paa,
           add_rasters = list(x = sst_historical, smallplot = sp, axis.args = axis_args, zlim = c(-2, 30), col = cols), 
           add_polys = list(x = land))
mtext(side = 3, "A (SST)", adj = adj_main, line = line_main, cex = cex_main, font = 2)

#### sbt_historical
raster::cellStats(sbt_historical, range)
pretty_map(sbt_historical,
           pretty_axis_args = paa,
           add_rasters = list(x = sbt_historical, smallplot = sp, axis.args = axis_args, zlim = c(-2, 30), col = cols), 
           add_polys = list(x = land))
mtext(side = 3, "B (SBT)", adj = adj_main, line = line_main, cex = cex_main, font = 2)
mtext(side = 4, expression(paste("Temperature (", degree, "C)")), line = 2.5, cex = 1.5, outer = TRUE)

#### Save
dev.off()


##############################
##############################
#### SST 

##############################
#### Exploration 

#### Define rasters
sst_1 <- sst_mid_rcp45 - sst_historical
sst_2 <- sst_mid_rcp85 - sst_historical
sst_3 <- (sst_mid_rcp85 - sst_historical) - (sst_mid_rcp45 - sst_historical)
sst_4 <- sst_late_rcp45 - sst_historical
sst_5 <- sst_late_rcp85 - sst_historical
sst_6 <- (sst_late_rcp85 - sst_historical) - (sst_late_rcp45 - sst_historical)
sst_7 <- (sst_late_rcp45 - sst_historical) - (sst_mid_rcp45 - sst_historical)
sst_8 <- (sst_late_rcp85 - sst_historical) - (sst_mid_rcp85 - sst_historical)
sst_n <- raster::stack(sst_1, sst_2, sst_3, sst_4, sst_5, sst_6, sst_7, sst_8)

#### Summary statistics
summarise_raster(sst_historical)
summarise_raster(sst_1)
summarise_raster(sst_2)
summarise_raster(sst_3)
summarise_raster(sst_4)
summarise_raster(sst_5)
summarise_raster(sst_6)
summarise_raster(sst_7)
summarise_raster(sst_8)

#### Quick visualisation 
raster::plot(sst_1)
raster::plot(sst_4)
raster::plot(sst_2)
raster::plot(sst_5)

#### ssplot() method 
# Define names
titles <- c("A, RCP 4.5 [M] - BL", 
            "B, RCP 8.5 [M] - BL", 
            "C, (RCP 8.5 [M] - BL) - (RCP 4.5 [M] - BL)", 
            "D, RCP 4.5 [L] - BL", 
            "E, RCP 8.5 [L] - BL", 
            "F, (RCP 8.5 [L] - BL) - (RCP 4.5 [L] - BL)",
            "G, (RCP 4.5 [L] - BL) - (RCP 4.5 [M] - BL)", 
            "H, (RCP 8.5 [L] - BL) - (RCP 8.5 [M] - BL)")
# Make plot
raster::spplot(sst_n, 
               names.attr = titles, 
               layout = c(3, 3), 
               par.settings = list(layout.widths = list(axis.key.padding = 4))) +
  latticeExtra::layer(sp::sp.polygons(land, col = "black", fill = "white", lwd = 2))


##############################
#### Publication-quality visualisation 

#### Set up plot to save 
tiff("./fig/sst_projections.tiff", 
     height = 6, width = 12, units = "in", res = 600)
pp <- par(mfrow = c(2, 2), oma = c(0, 0, 1, 5), mar = c(0, 1.25, 0, 1))

#### Param
line_main <- -0.9
adj_main  <- 0.05
cex_main  <- 1.5
sp <- c(0.95, 0.96, 0.2, 0.8)
# cols <- RColorBrewer::brewer.pal(3, "YlOrRd")
# cols <- colorRampPalette(cols)(100)
zlim_mid <- zlim_late <- rep(NA, 2) 
zlim_mid[1]  <- min(sapply(list(sst_1, sst_2), raster::cellStats, min))
zlim_mid[2]  <- max(sapply(list(sst_1, sst_2), raster::cellStats, max))
zlim_late[1] <- min(sapply(list(sst_4, sst_5), raster::cellStats, min))
zlim_late[2] <- max(sapply(list(sst_4, sst_5), raster::cellStats, max))
# Define split heat colour schemes 
col_param_mid  <- pretty_cols_split_heat(zlim = zlim_mid, hot = "YlOrRd")
col_param_late <- pretty_cols_split_heat(zlim = zlim_late, hot = "YlOrRd")

#### delta sst_mid_rcp45
raster::cellStats(sst_1, range)
pretty_map(sst_1, 
           pretty_axis_args = paa,
           add_rasters = list(x = sst_1, smallplot = sp, axis.args = axis_args, 
                              zlim = zlim_mid, breaks = col_param_mid$breaks, col = col_param_mid$col), 
           add_polys = list(x = land))
mtext(side = 3, "A (RCP 4.5 [M] - SST)", adj = adj_main, line = line_main, cex = cex_main, font = 2)

#### delta sst_late_rcp45
raster::cellStats(sst_4, range)
pretty_map(sst_4,
           pretty_axis_args = paa,
           add_rasters = list(x = sst_4, smallplot = sp, axis.args = axis_args, 
                              zlim = zlim_late, breaks = col_param_late$breaks, col = col_param_late$col), 
           add_polys = list(x = land))
mtext(side = 3, "B (RCP 4.5 [L] - SST)", adj = adj_main, line = line_main, cex = cex_main, font = 2)

#### delta sst_mid_rcp45
# Define colour scheme split around 0 
raster::cellStats(sst_2, range)
zlim <- raster::cellStats(sst_2, range)
pretty_map(sst_2, 
           pretty_axis_args = paa,
           add_rasters = list(x = sst_2, smallplot = sp, axis.args = axis_args, 
                              zlim = zlim_mid, breaks = col_param_mid$breaks, col = col_param_mid$col), 
           add_polys = list(x = land))
mtext(side = 3, "C (RCP 8.5 [M] - SST)", adj = adj_main, line = line_main, cex = cex_main, font = 2)

#### delta sst_late_rcp85
raster::cellStats(sst_5, range)
pretty_map(sst_5, 
           pretty_axis_args = paa,
           add_rasters = list(x = sst_5, smallplot = sp, axis.args = axis_args, 
                              zlim = zlim_late, breaks = col_param_late$breaks, col = col_param_late$col), 
           add_polys = list(x = land))
mtext(side = 3, "D (RCP 8.5 [L] - SST)", adj = adj_main, line = line_main, cex = cex_main, font = 2)

#### Save plot 
# mtext(side = 1, expression(paste("Longitude (", degree, "C)")), line = 1, cex = 1.5, outer = TRUE)
# mtext(side = 2, expression(paste("Latitude (", degree, "C)")), line = -2, cex = 1.5, outer = TRUE)
mtext(side = 4, expression(paste("Temperature (", degree, "C)")), line = 2.5, cex = 1.5, outer = TRUE)
par(pp)
dev.off()


##############################
##############################
#### SBT

############################## 
#### Exploration 

#### Define rasters
sbt_1 <- sbt_mid_rcp45 - sbt_historical
sbt_2 <- sbt_mid_rcp85 - sbt_historical
sbt_3 <- (sbt_mid_rcp85 - sbt_historical) - (sbt_mid_rcp45 - sbt_historical)
sbt_4 <- sbt_late_rcp45 - sbt_historical
sbt_5 <- sbt_late_rcp85 - sbt_historical
sbt_6 <- (sbt_late_rcp85 - sbt_historical) - (sbt_late_rcp45 - sbt_historical)
sbt_7 <- (sbt_late_rcp45 - sbt_historical) - (sbt_mid_rcp45 - sbt_historical)
sbt_8 <- (sbt_late_rcp85 - sbt_historical) - (sbt_mid_rcp85 - sbt_historical)
sbt_n <- raster::stack(sbt_1, sbt_2, sbt_3, sbt_4, sbt_5, sbt_6, sbt_7, sbt_8)

#### Raster summarise
summarise_raster(sbt_historical)
summarise_raster(sbt_1)
summarise_raster(sbt_2)
summarise_raster(sbt_3)
summarise_raster(sbt_4)
summarise_raster(sbt_5)
summarise_raster(sbt_6)
summarise_raster(sbt_7)
summarise_raster(sbt_8)

#### Quick visualisation
raster::plot(sbt_1)
zlim <- raster::cellStats(sbt_1, range)
col_param <- pretty_cols_split_heat(zlim)
fields::image.plot(sbt_1, zlim = zlim, breaks = col_param$breaks, col = col_param$col)
raster::plot(sbt_4)
raster::plot(sbt_2)
raster::plot(sbt_5)

##############################
#### Publication-quality visualisation 

#### Set up plot to save 
tiff("./fig/sbt_projections.tiff", 
     height = 6, width = 12, units = "in", res = 600)
pp <- par(mfrow = c(2, 2), oma = c(0, 0, 1, 5), mar = c(0, 1.5, 0, 1))

#### Param
# Define limits 
zlim_mid <- zlim_late <- rep(NA, 2) 
zlim_mid[1]  <- min(sapply(list(sbt_1, sbt_2), raster::cellStats, min))
zlim_mid[2]  <- max(sapply(list(sbt_1, sbt_2), raster::cellStats, max))
zlim_late[1] <- min(sapply(list(sbt_4, sbt_5), raster::cellStats, min))
zlim_late[2] <- max(sapply(list(sbt_4, sbt_5), raster::cellStats, max))
# Define split heat colour schemes 
col_param_mid  <- pretty_cols_split_heat(zlim = zlim_mid, hot = "YlOrRd")
col_param_late <- pretty_cols_split_heat(zlim = zlim_late, hot = "YlOrRd")

#### delta sbt_mid_rcp45
raster::cellStats(sbt_1, range)
pretty_map(sbt_1, 
           pretty_axis_args = paa,
           add_rasters = list(x = sbt_1, smallplot = sp, axis.args = axis_args, 
                              zlim = zlim_mid, breaks = col_param_mid$breaks, col = col_param_mid$col), 
           add_polys = list(x = land))
mtext(side = 3, "A (RCP 4.5 [M] - SBT)", adj = adj_main, line = line_main, cex = cex_main, font = 2)

#### delta sbt_late_rcp45
raster::cellStats(sbt_4, range)
pretty_map(sbt_4,
           pretty_axis_args = paa,
           add_rasters = list(x = sbt_4, smallplot = sp, axis.args = axis_args, 
                              zlim = zlim_late, breaks = col_param_late$breaks, col = col_param_late$col), 
           add_polys = list(x = land))
mtext(side = 3, "B (RCP 4.5 [L] - SBT)", adj = adj_main, line = line_main, cex = cex_main, font = 2)

#### delta sbt_mid_rcp45
# Define colour scheme split around 0 
raster::cellStats(sbt_2, range)
zlim <- raster::cellStats(sbt_2, range)
pretty_map(sbt_2, 
           pretty_axis_args = paa,
           add_rasters = list(x = sbt_2, smallplot = sp, axis.args = axis_args, 
                              zlim = zlim_mid, breaks = col_param_mid$breaks, col = col_param_mid$col), 
           add_polys = list(x = land))
mtext(side = 3, "C (RCP 8.5 [M] - SBT)", adj = adj_main, line = line_main, cex = cex_main, font = 2)

#### delta sbt_late_rcp85
raster::cellStats(sbt_5, range)
pretty_map(sbt_5, 
           pretty_axis_args = paa,
           add_rasters = list(x = sbt_5, smallplot = sp, axis.args = axis_args, 
                              zlim = zlim_late, breaks = col_param_late$breaks, col = col_param_late$col), 
           add_polys = list(x = land))
mtext(side = 3, "D (RCP 8.5 [L] - SBT)", adj = adj_main, line = line_main, cex = cex_main, font = 2)

#### Save plot 
# mtext(side = 1, expression(paste("Longitude (", degree, "C)")), line = 1, cex = 1.5, outer = TRUE)
# mtext(side = 2, expression(paste("Latitude (", degree, "C)")), line = -2, cex = 1.5, outer = TRUE)
mtext(side = 4, expression(paste("Temperature (", degree, "C)")), line = 2.5, cex = 1.5, outer = TRUE)
par(pp)
dev.off()


##############################
##############################
#### SST - SBT differences

#### Baseline temperatures
raster::plot(sst_historical)
raster::plot(sbt_historical)
raster::plot(sst_historical - sbt_historical)
summarise_raster(sst_historical - sbt_historical)



#### End of code. 
##############################
##############################