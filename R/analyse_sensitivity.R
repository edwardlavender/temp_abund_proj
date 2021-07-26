##############################
##############################
#### process_sensitivity.R 
#### Edward Lavender (el72@st-andrews.ac.uk)

#### This code:
# 1) Maps species' sensitivity metrics over space. For SST and SBT, this includes:
# ... mean STR
# ... SD STR
# ... Mean thermal bias (STI - SST)
# ... SD thermal bias (STI - SST)

#### Steps preceding this code:
# 1) Define species sensitivity metrics (process_sensitivity.R)


##############################
##############################
#### Set up 

#### Wipe workspace
rm(list = ls())

#### Essential packages
library(prettyGraphics)

#### Global param
root_sensitivity_sst <- "./data/sensitivity/sst/"
root_sensitivity_sbt <- "./data/sensitivity/sbt/"
source("./R/helpers.R")

#### Load files
## Coastline
coastline <- readRDS("./data/spatial/coastline/coastline.rds")
land <- flapper::invert_poly(coastline)
cover <- raster::raster("./data/spatial/eez/eez_mask.asc")
## SST
sst_str_mean <- raster::raster(paste0(root_sensitivity_sst, "sst_str_mean.asc"))
sst_tb_mean  <- raster::raster(paste0(root_sensitivity_sst, "sst_tb_mean.asc"))
sst_tb_iqr   <- raster::raster(paste0(root_sensitivity_sst, "sst_tb_iqr.asc"))
## SBT
sbt_str_mean <- raster::raster(paste0(root_sensitivity_sbt, "sbt_str_mean.asc"))
sbt_tb_mean  <- raster::raster(paste0(root_sensitivity_sbt, "sbt_tb_mean.asc"))
sbt_tb_iqr   <- raster::raster(paste0(root_sensitivity_sbt, "sbt_tb_iqr.asc"))


##############################
##############################
#### Processing 

#### Mask rasters to focus on coastal areas
## SST 
sst_str_mean <- raster::mask(sst_str_mean, cover)
sst_tb_mean  <- raster::mask(sst_tb_mean, cover)
sst_tb_iqr   <- raster::mask(sst_tb_iqr, cover)
## SBT 
sbt_str_mean <- raster::mask(sbt_str_mean, cover)
sbt_tb_mean  <- raster::mask(sbt_tb_mean, cover)
sbt_tb_iqr   <- raster::mask(sbt_tb_iqr, cover)


##############################
##############################
#### Quick visualisation 

### Quick plots
# STR
rasterVis::levelplot(sst_str_mean, margin = list(FUN = mean))
rasterVis::levelplot(sbt_str_mean, margin = list(FUN = mean))
# TB
rasterVis::levelplot(sst_tb_mean, margin = list(FUN = mean))
rasterVis::levelplot(sbt_tb_mean, margin = list(FUN = mean))
rasterVis::levelplot(sst_tb_iqr, margin = list(FUN = mean))
rasterVis::levelplot(sbt_tb_iqr, margin = list(FUN = mean))


##############################
##############################
#### STR and TB

#### Define graphical param 
line_main <- -0.9
adj_main  <- 0.05
cex_main  <- 1.5

#### Set up plot to save 
png("./fig/sensitivity.png", 
     height = 5.5, width = 12, units = "in", res = 600)
pp <- par(mfrow = c(2, 2), oma = c(0, 0, 1, 5), mar = c(0, 2.45, 0, 2.45))
sp <- c(0.99, 1, 0.2, 0.8)

#### Get sensible z limits to faciliate comparisons
# STR
zlim_str <- range(c(raster::cellStats(sst_str_mean, range), 
                    raster::cellStats(sbt_str_mean, range)))
zlim_str
# TB
zlim_tb <- range(c(raster::cellStats(sst_tb_mean, range), 
                   raster::cellStats(sbt_tb_mean, range)))
zlim_tb

#### sst_str_mean
plot_raster(sst_str_mean, zlim = zlim_str, select = 2:8, rev = TRUE, sp = sp)
mtext(side = 3, "A (mean SST STR)", adj = adj_main, line = line_main, cex = cex_main, font = 2)

#### sbt_str_mean
plot_raster(sbt_str_mean, zlim = zlim_str, select = 2:8, rev = TRUE, sp = sp)
mtext(side = 3, "B (mean SBT STR)", adj = adj_main, line = line_main, cex = cex_main, font = 2)

#### sst_tb_mean
plot_raster(sst_tb_mean, 
            zlim = zlim_tb,
            gen_cols = pretty_cols_split_heat,
            scheme_hot = "Blues",  scheme_cold = "YlOrRd", 
            select_hot = 4:8, select_cold = 4:8,
            sp = sp)
mtext(side = 3, "C (mean SST STB)", adj = adj_main, line = line_main, cex = cex_main, font = 2)

#### sbt_tb_mean
plot_raster(sbt_tb_mean, 
            zlim = zlim_tb,
            gen_cols = pretty_cols_split_heat,
            scheme_hot = "Blues", scheme_cold = "YlOrRd", 
            select_hot = 4:8, select_cold = 4:8,
            sp = sp)
mtext(side = 3, "D (mean SBT STB)", adj = adj_main, line = line_main, cex = cex_main, font = 2)

#### Save plot 
# mtext(side = 1, expression(paste("Longitude (", degree, "C)")), line = 1, cex = 1.5, outer = TRUE)
# mtext(side = 2, expression(paste("Latitude (", degree, "C)")), line = -2, cex = 1.5, outer = TRUE)
mtext(side = 4, expression(paste("Temperature (", degree, "C)")), line = 3.5, cex = 1.5, outer = TRUE)
par(pp)
dev.off()


##############################
##############################
#### Variability in thermal bias 

#### Define graphical param 
line_main <- -0.8
adj_main  <- 0.07
cex_main  <- 1.5
zlim_tb_iqr <- range(c(raster::cellStats(sst_tb_iqr, range),
                       raster::cellStats(sst_tb_iqr, range)))

#### Set up plot to save 

png("./fig/sensitivity_tb_iqr.png", 
     height = 5.5, width = 6.75, units = "in", res = 600)
pp <- par(mfrow = c(2, 1), oma = c(0, 0, 1, 5), mar = c(0, 2.1, 0, 2.1))

#### sst_str_mean
plot_raster(sst_tb_iqr, zlim = zlim_tb_iqr, 
            select = 2:8, rev = TRUE, profile_x = c(185, 220), sp = sp)
mtext(side = 3, "A (IQR STI - SST)", adj = adj_main, line = line_main, cex = cex_main, font = 2)

#### sbt_str_mean
plot_raster(sbt_tb_iqr, zlim = zlim_tb_iqr, 
            select = 2:8, rev = TRUE, profile_x = c(185, 220), sp = sp)
mtext(side = 3, "B (IQR STI - SBT)", adj = adj_main, line = line_main, cex = cex_main, font = 2)

#### Save plot 
# mtext(side = 1, expression(paste("Longitude (", degree, "C)")), line = 1, cex = 1.5, outer = TRUE)
# mtext(side = 2, expression(paste("Latitude (", degree, "C)")), line = -2, cex = 1.5, outer = TRUE)
mtext(side = 4, expression(paste("Temperature (", degree, "C)")), line = 3.5, cex = 1.5, outer = TRUE)
par(pp)
dev.off()


#### End of code. 
##############################
##############################