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
## SST
sst_str_mean <- raster::raster(paste0(root_sensitivity_sst, "sst_str_mean.asc"))
sst_str_sd   <- raster::raster(paste0(root_sensitivity_sst, "sst_str_sd.asc"))
sst_tb_mean  <- raster::raster(paste0(root_sensitivity_sst, "sst_tb_mean.asc"))
sst_tb_sd    <- raster::raster(paste0(root_sensitivity_sst, "sst_tb_sd.asc"))
## SBT
sbt_str_mean <- raster::raster(paste0(root_sensitivity_sbt, "sbt_str_mean.asc"))
sbt_str_sd   <- raster::raster(paste0(root_sensitivity_sbt, "sbt_str_sd.asc"))
sbt_tb_mean  <- raster::raster(paste0(root_sensitivity_sbt, "sbt_tb_mean.asc"))
sbt_tb_sd    <- raster::raster(paste0(root_sensitivity_sbt, "sbt_tb_sd.asc"))


##############################
##############################
#### Quick visualisation 

### Quick plots
# STR
rasterVis::levelplot(sst_str_mean, margin = list(FUN = mean))
rasterVis::levelplot(sbt_str_mean, margin = list(FUN = mean))
rasterVis::levelplot(sst_str_sd, margin = list(FUN = mean))
rasterVis::levelplot(sbt_str_sd, margin = list(FUN = mean))
# TB
rasterVis::levelplot(sst_tb_mean, margin = list(FUN = mean))
rasterVis::levelplot(sbt_tb_mean, margin = list(FUN = mean))
rasterVis::levelplot(sst_tb_sd, margin = list(FUN = mean))
rasterVis::levelplot(sbt_tb_sd, margin = list(FUN = mean))

#### Group plots 
raster::plot(raster::stack(sst_str_mean, sst_str_sd, sst_tb_mean, sst_tb_sd))
raster::plot(raster::stack(sbt_str_mean, sbt_str_sd, sbt_tb_mean, sbt_tb_sd))


##############################
##############################
#### Make maps 

##############################
#### SST 

#### Define graphical param 
line_main <- -0.9
adj_main  <- 0.05
cex_main  <- 1.5

#### Set up plot to save 
tiff("./fig/sst_sensitivity.tiff", 
     height = 5.5, width = 12, units = "in", res = 600)
pp <- par(mfrow = c(2, 2), oma = c(0, 0, 1, 5), mar = c(0, 2.1, 0, 2.1))

#### sst_str_mean
plot_raster(sst_str_mean)
mtext(side = 3, "A (mean STR)", adj = adj_main, line = line_main, cex = cex_main, font = 2)

#### sst_str_mean
plot_raster(sst_str_sd)
mtext(side = 3, "B (SD STR)", adj = adj_main, line = line_main, cex = cex_main, font = 2)

#### sst_tb_mean
plot_raster(sst_tb_mean, 
            gen_cols = pretty_cols_split_heat,
            scheme_hot = "Blues", 
            scheme_cold = "YlOrRd")
mtext(side = 3, "C (mean STI - SST)", adj = adj_main, line = line_main, cex = cex_main, font = 2)

#### sst_tb_sd
plot_raster(sst_tb_sd)
mtext(side = 3, "D (SD STI - SST)", adj = adj_main, line = line_main, cex = cex_main, font = 2)

#### Save plot 
# mtext(side = 1, expression(paste("Longitude (", degree, "C)")), line = 1, cex = 1.5, outer = TRUE)
# mtext(side = 2, expression(paste("Latitude (", degree, "C)")), line = -2, cex = 1.5, outer = TRUE)
mtext(side = 4, expression(paste("Temperature (", degree, "C)")), line = 2.75, cex = 1.5, outer = TRUE)
par(pp)
dev.off()


##############################
#### SBT 
# ... [copied from above but 'sst' replaced by 'sbt']

#### Set up plot to save 
tiff("./fig/sbt_sensitivity.tiff", 
     height = 5.5, width = 12, units = "in", res = 600)
pp <- par(mfrow = c(2, 2), oma = c(0, 0, 1, 5), mar = c(0, 2.1, 0, 2.1))

#### sbt_str_mean
plot_raster(sbt_str_mean)
mtext(side = 3, "A (mean STR)", adj = adj_main, line = line_main, cex = cex_main, font = 2)

#### sbt_str_mean
plot_raster(sbt_str_sd)
mtext(side = 3, "B (SD STR)", adj = adj_main, line = line_main, cex = cex_main, font = 2)

#### sbt_tb_mean
plot_raster(sbt_tb_mean, 
            gen_cols = pretty_cols_split_heat,
            scheme_hot = "Blues", 
            scheme_cold = "YlOrRd")
mtext(side = 3, "C (mean STI - SBT)", adj = adj_main, line = line_main, cex = cex_main, font = 2)

#### sbt_tb_sd
plot_raster(sbt_tb_sd)
mtext(side = 3, "D (SD STI - SBT)", adj = adj_main, line = line_main, cex = cex_main, font = 2)

#### Save plot 
# mtext(side = 1, expression(paste("Longitude (", degree, "C)")), line = 1, cex = 1.5, outer = TRUE)
# mtext(side = 2, expression(paste("Latitude (", degree, "C)")), line = -2, cex = 1.5, outer = TRUE)
mtext(side = 4, expression(paste("Temperature (", degree, "C)")), line = 2.5, cex = 1.75, outer = TRUE)
par(pp)
dev.off()



#### End of code. 
##############################
##############################