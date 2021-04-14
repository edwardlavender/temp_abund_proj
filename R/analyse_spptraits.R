##############################
##############################
#### analyse_spptraits.R
#### Edward Lavender (el72@st-andrews.ac.uk)

#### This code:
# Analyses the species selected for modelling.

#### Steps preceding this code:
# 1) process_spptraits.R defines a list of species for modelling 


##############################
##############################
#### Set up 

#### Wipe workspace
rm(list = ls())

#### Essential packages
library(magrittr)
library(prettyGraphics)

#### Load data
spptraits <- readRDS("./data/spptraits.rds")
coastline <- readRDS("./data/spatial/coastline/coastline.rds")


##############################
##############################
#### Table of species 

#### Number of species
nrow(spptraits)

#### Check high-level taxonomic frequencies
table(spptraits$phylum) # (all Chordata)
table(spptraits$class)  # (all Actinopteri)
table(spptraits$order)  

#### Define table 
spp_tbl <- 
  spptraits %>%
  dplyr::select(order, family, spp, occupancy) %>%
  dplyr::arrange(order, family, spp)


##############################
##############################
#### Depth preferences 

#### Get depth quantiles 
qs <- quantile(spptraits$depth_range_shallow,  probs = seq(0, 1, by = 0.05), na.rm = TRUE) 
qd <- quantile(spptraits$depth_range_deep,  probs = seq(0, 1, by = 0.05))    
qs
qd

#### Make figure

## Set up figure to save 
tiff("./fig/depth_quantiles.tiff", height = 5, width = 5, units = "in", res = 600)
cex_axis <- 1.1
cex_lab  <- 1.1

## Define shallow quantiles 
qs <- data.frame(pc = names(qs), depth = qs, row.names = NULL)
qs$pc <- as.character(qs$pc)
qs$pc <- substr(qs$pc, 1, nchar(qs$pc)-1)
qs$pc <- as.numeric(qs$pc)/100

## Define deep quantiles
qd <- data.frame(pc = names(qd), depth = qd, row.names = NULL)
qd$pc <- as.character(qd$pc)
qd$pc <- substr(qd$pc, 1, nchar(qd$pc)-1)
qd$pc <- as.numeric(qd$pc)/100

## Plot 
plot(qs$depth,qs$pc, type = "l", 
     axes = F, xlim = c(0, 50), ylim = c(0, 1), 
     xlab = "", ylab = "", 
     lwd = 2)
lines(qd$depth, qd$pc, lwd = 2, lty = 3)

## Add titles 
mtext(side = 1, "Depth Limit (m)", line = 2.5, cex = cex_lab)
mtext(side = 2, "Quantile", line = 2.5, cex = cex_lab)
axis(side = 1, seq(0, 50, by = 10), pos = 0, cex = cex_axis)
axis(side = 2, seq(0, 1, by = 0.2), pos = 0,  las = 2, cex = cex_axis)

## Add legend
legend(30, 0.3, lty = c(1, 3), lwd = 2, 
       legend = c("Shallow", "Deep"), bty = "n", 
       y.intersp = 1.5, cex = cex_axis)

dev.off()


##############################
##############################
#### Species richness

#### Method used to map species richness
# 1) Load rasters for all species into R
# ... This is the slow step, so we will do this in parallel
# ... We don't have to do this beforehand (e.g., if insufficient memory), but it is faster
# 2) Define a blank map, and update this by adding predictions for each species
# ... We could stack the raster, and then sum, but this is probably slower. 

#### Load rasters into a list 
# ... [This takes < 10 seconds on 8 cores]
cl <- parallel::makeCluster(8L)
spptraits$index <- 1:nrow(spptraits)
maps_by_spp <- pbapply::pblapply(split(spptraits, spptraits$index), cl = cl, function(d){
  raster::raster(paste0("./data/sdm_aquamaps/", d$number.asc.file))
})
parallel::stopCluster(cl)

#### Define blank map
spp_map <- raster::raster(paste0("./data/sdm_aquamaps/", spptraits$number.asc.file[1]))
spp_map <- raster::setValues(spp_map, 0)
raster::plot(spp_map)

#### Update map sequentially, by species
# [~ Time difference of 13 minutes]
t1 <- Sys.time()
for(i in 1:length(maps_by_spp)){
  svMisc::progress(i, length(maps_by_spp))
  spp_map <- sum(spp_map, maps_by_spp[[i]])
}
t2 <- Sys.time()
difftime(t2, t1)
# raster::writeRaster(spp_map, "./data/spatial/map_spp_richness.asc")

#### Plot map
# Set up figure
tiff("./fig/spp_richness_map.tiff", 
     height = 8, width = 10, units = "in", res = 600)
pp <- par(oma = c(1, 1, 1, 5))
cex_axis <- 1.25
cex_lab  <- 1.5
# Create plot 
spp_map[spp_map == 0] <- NA
raster::plot(spp_map, 
             axes = FALSE, box = FALSE,
             # col = rev(viridis::viridis(100)),
             xlim = c(-180, 180), ylim = c(-90, 90), zlim = c(0, 1400), 
             axis.args = list(cex.axis = cex_axis))
land <- flapper::invert_poly(coastline)
raster::plot(sea, col = "white", add = TRUE)
# Add axes 
xat <- seq(-180, 180, by = 60)
yat <- seq(-90, 90, by = 30)
axis(side = 1, xat, pos = -90, cex.axis = cex_axis)    
axis(side = 2, yat, pos = -180, cex.axis= cex_axis, las = TRUE)   
axis(side = 3, xat, labels = FALSE, pos = 90)
axis(side = 4, yat, labels = FALSE, pos = 180)
# Add labels
mtext(side = 1, expression(paste("Longitude (", degree, ")")), line = -2.5, cex = cex_lab)
mtext(side = 2, expression(paste("Latitude (", degree, ")")), line = 3, cex = cex_lab)
mtext(side = 4, "Species richness (count)", line = 6, cex = cex_lab)
# Save 
dev.off()


##############################
##############################
#### Fisheries statistics

# Number of species with fisheries statistics
table(!is.na(spptraits$importance))
# % of species with fisheries statistics 
table(!is.na(spptraits$importance))/nrow(spptraits)*100
table(spptraits$importance)/table(!is.na(spptraits$importance))[2]*100


##############################
##############################
#### Thermal preferences


##############################
#### Basic stats

#### Thermal param 
spptraits$sst_str <- spptraits$sst_t90 - spptraits$sst_t10
spptraits$sbt_str <- spptraits$sbt_t90 - spptraits$sbt_t10

#### Examine ranges 
# T10
utils.add::basic_stats(spptraits$sst_t10)
utils.add::basic_stats(spptraits$sbt_t10)
utils.add::basic_stats(spptraits$sst_t10 - spptraits$sbt_t10)
# T50
utils.add::basic_stats(spptraits$sst_t50)
utils.add::basic_stats(spptraits$sbt_t50)
utils.add::basic_stats(spptraits$sst_t50 - spptraits$sbt_t50)
# T90 
utils.add::basic_stats(spptraits$sst_t90)
utils.add::basic_stats(spptraits$sbt_t90)
utils.add::basic_stats(spptraits$sst_t90 - spptraits$sbt_t90)
# STR 
utils.add::basic_stats(spptraits$sst_str - spptraits$sbt_str)

#### Check species and distributions for species with min/pax param
# sst_t10, sbt_t10 (MIN) --> ISSUE 
spptraits[which.min(spptraits$sst_t10), ]
spptraits[which.min(spptraits$sbt_t10), ]
raster::plot(raster::raster("./data-raw/sdm_aquamaps/13059.asc"))
# sbt_t10 (MAX)
spptraits[which.max(spptraits$sst_t10), ]
spptraits[which.max(spptraits$sbt_t10), ]
raster::plot(raster::raster("./data-raw/sdm_aquamaps/61353.asc")); raster::lines(coastline)
raster::plot(raster::raster("./data-raw/sdm_aquamaps/55845.asc")); raster::lines(coastline)

# sst_t90, sbt_t90 (MIN) --> ISSUE
spptraits[which.min(spptraits$sst_t90), ]
spptraits[which.min(spptraits$sbt_t10), ]
raster::plot(raster::raster("./data-raw/sdm_aquamaps/56311.asc")); raster::lines(coastline)
# sst_t90, sbt_t90 (MAX)
spptraits[which.max(spptraits$sst_t90), ]
spptraits[which.max(spptraits$sbt_t90), ]
raster::plot(raster::raster("./data-raw/sdm_aquamaps/15311.asc")); raster::lines(coastline)
raster::plot(raster::raster("./data-raw/sdm_aquamaps/55845.asc")); raster::lines(coastline)

#### Correlation between SST and SBT thermal niches
# Models
mod_t10 <- lm(sst_t10 ~ sbt_t10, data = spptraits)
mod_t50 <- lm(sst_t50 ~ sbt_t50, data = spptraits)
mod_t90 <- lm(sst_t90 ~ sbt_t90, data = spptraits)
# Summaries
summary(mod_t10)
summary(mod_t50)
summary(mod_t90)


##############################
#### Plot correlations

#### Set up figure
tiff("./fig/thermal_affinity_correlations.tiff", 
     height = 6, width = 6, units = "in", res = 600)
pp <- par(mfrow = c(2, 2), oma = c(2.5, 2.5, 2, 2), mar = c(1, 2, 2, 1))
xlim <- ylim <- c(-5, 30)

#### T10
pretty_plot(spptraits$sbt_t10, spptraits$sst_t10, 
            xlab = "", ylab = "",
            xlim = xlim, ylim = ylim,
            col = scales::alpha("black", 0.75), cex = 0.8)

lines(xlim, ylim)
mtext(side = 3, "A [T10]", adj = 0, font = 2, cex = 1.25)

#### STR
pretty_plot(spptraits$sbt_t50, spptraits$sst_t50, 
            xlab = "", ylab = "",
            xlim = xlim, ylim = ylim,
            col = scales::alpha("black", 0.75), cex = 0.8)
lines(xlim, ylim)
mtext(side = 3, "B [STI]", adj = 0, font = 2, cex = 1.25)

#### T90
pretty_plot(spptraits$sbt_t90, spptraits$sst_t90, 
            xlab = "", ylab = "",
            xlim = xlim, ylim = ylim,
            col = scales::alpha("black", 0.75), cex = 0.8)
lines(xlim, ylim)
mtext(side = 3, "C [T90]", adj = 0, font = 2, cex = 1.25)

#### STR
pretty_plot(spptraits$sbt_str, spptraits$sst_str, 
            xlab = "", ylab = "",
            xlim = xlim, ylim = ylim,
            col = scales::alpha("black", 0.75), cex = 0.8)
lines(xlim, ylim)
mtext(side = 3, "D [STR]", adj = 0, font = 2, cex = 1.25)

#### Save figure 
mtext(side = 1, expression(paste("SBT (", degree, "C)")), line = 1.5, cex = 1.5, outer = TRUE)
mtext(side = 2, expression(paste("SST (", degree, "C)")), line = 0.25, cex = 1.5, outer = TRUE)
par(pp)
dev.off()


#### End of code. 
##############################
##############################
