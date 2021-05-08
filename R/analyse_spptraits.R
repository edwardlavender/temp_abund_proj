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
source("./R/helpers.R")

#### Load data
spptraits <- readRDS("./data/spptraits.rds")
coastline <- readRDS("./data/spatial/coastline/coastline.rds")
eez <- sf::read_sf("./data/spatial/eez", "eez_boundaries_v11") 
eez <- as(eez, "Spatial")


##############################
##############################
#### Table of species 

#### Number of species
nrow(spptraits)

#### Check high-level taxonomic frequencies
table(spptraits$class)  
# Actinopterygii Elasmobranchii 
#    2271             22
table(spptraits$order)  

#### Occur cells 
utils.add::basic_stats(spptraits$occur_cells)

#### Define table 
# Define table
spp_tbl <- 
  spptraits %>%
  dplyr::select(order, family, genus, species_epiphet, occur_cells, sst_t10, sst_t50, sst_t90) %>%
  dplyr::arrange(order, family, genus, species_epiphet)
colnames(spp_tbl) <- c("Order", "Family", "Genus", "Species", "Cells", "T10", "T50", "T90")
# Round thermal affinities
spp_tbl$T10 <- prettyGraphics::add_lagging_point_zero(round(spp_tbl$T10, 2), 2)
spp_tbl$T50 <- prettyGraphics::add_lagging_point_zero(round(spp_tbl$T50, 2), 2)
spp_tbl$T90 <- prettyGraphics::add_lagging_point_zero(round(spp_tbl$T90, 2), 2)
# 'Cells' is the number of unique 0.5 degree cells containing valid occurrences ('occurcells' parameter)
# Aquamaps warns against using maps with fewer than 10 'occurcells'

#### Save table
head(spp_tbl)
# write.table(spp_tbl, "./fig/spp_tbl.txt", sep = ",", row.names = FALSE, quote = FALSE)


##############################
##############################
#### Depth preferences 

#### Lifestyle
get_habits <- FALSE
if(get_habits){
  #### Get fishbase data 
  # Get data 
  fb <- rfishbase::species()
  # Check proportion of species with data 
  table(!is.na(fb$Importance))/nrow(fb)*100
  # Check background importance proportions - is there a bias towards commercial species? 
  table(fb$Importance)/sum(table(fb$Importance))*100
  
  #### Get fishbase data for modelled species
  spptraits$habit <- fb$DemersPelag[match(spptraits$spp, fb$Species)]
  table(is.na(spptraits$habit))
  table(spptraits$habit)
}

#### Get depth quantiles 
table(!is.na(spptraits$depth_range_shallow))
spptraits$depth_mean <- apply(spptraits[, c("depth_range_shallow", "depth_range_deep")], 1, mean)
utils.add::basic_stats(spptraits$depth_mean, na.rm = TRUE)
qs <- quantile(spptraits$depth_range_shallow,  probs = seq(0, 1, by = 0.05), na.rm = TRUE) 
qd <- quantile(spptraits$depth_range_deep,  probs = seq(0, 1, by = 0.05))    
qs
qd

#### Make figure

## Set up figure to save 
png("./fig/depth_quantiles.png", height = 5, width = 5, units = "in", res = 600)
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

##############################
#### Species richness (based on threshold maps)

run <- FALSE
if(run){
  
  ## Define a list of species' rasters 
  # This takes ~ 8 s on 8 cores. 
  cl <- parallel::makeCluster(8L)
  spptraits$index <- 1:nrow(spptraits)
  maps_by_spp <- pbapply::pblapply(split(spptraits, spptraits$index), cl = cl, function(d){
    raster::raster(paste0("./data/sdm_aquamaps/maps_occ/", d$spp_key_asc))
  })
  parallel::stopCluster(cl)
  
  ## Calculate total number of species via calc_sum() 
  # This takes ~ 14 minutes. 
  spp_map <- calc_sum(maps_by_spp)
  
  ## Save
  raster::writeRaster(spp_map, "./data/spatial/species/map_spp_richness.asc", overwrite = TRUE)
  
} else {
  
  ## Load 
  spp_map <- raster::raster("./data/spatial/species/map_spp_richness.asc")
  
}


##############################
#### Species richness (based on probability)

#### Define a map of species richness 
run <- FALSE
if(run){
  
  ## Define a list of species' rasters 
  # This takes ~ 8 s on 8 cores. 
  cl <- parallel::makeCluster(8L)
  spptraits$index <- 1:nrow(spptraits)
  maps_by_spp <- pbapply::pblapply(split(spptraits, spptraits$index), cl = cl, function(d){
    raster::raster(paste0("./data/sdm_aquamaps/maps_pr/", d$spp_key_asc))
  })
  parallel::stopCluster(cl)
  
  ## Brick species' rasters
  # This approach is implemented so that we can, in theory, use the species richness brick
  # ... to weight abundance estimates. 
  # This takes 16.61 hours [...]
  t1_bk <- Sys.time()
  spp_bk <- raster::brick(maps_by_spp)
  raster::writeRaster(spp_bk, "./data/spatial/species/spp_pr_bk.tif", overwrite = TRUE)
  beepr::beep(10)
  t2_bk <- Sys.time()
  difftime(t2_bk, t1_bk)
  
  ## Calculate species richness 
  # This takes 1 minute 
  t1_bk_calc <- Sys.time()
  spp_map <- raster::calc(spp_bk, sum, na.rm = TRUE)
  t2_bk_calc <- Sys.time()
  difftime(t2_bk_calc, t1_bk_calc)
  spp_map[is.na(spp_map)] <- NA
  raster::plot(spp_map)
  raster::writeRaster(spp_map, "./data/spatial/species/map_spp_richness_pr.asc", overwrite = TRUE)
  
} 


##############################
#### Map species richness

#### Plot map
# Set up figure
png("./fig/spp_richness_map.png", 
     height = 8, width = 10, units = "in", res = 600)
pp <- par(oma = c(1, 1, 1, 5))
cex_axis <- 1.25
cex_lab  <- 1.5
# Create plot 
spp_map[spp_map == 0] <- NA
# spp_map <- sqrt(spp_map)
zlim <- raster::cellStats(spp_map, range)
# col_param <- pretty_cols_brewer(zlim = zlim, scheme = "YlOrRd")
raster::plot(spp_map, 
             axes = FALSE, box = FALSE,
             col = viridis::plasma(100),
             xlim = c(-180, 180), ylim = c(-90, 90), zlim = zlim, 
             axis.args = list(cex.axis = cex_axis))
land <- flapper::invert_poly(coastline)
raster::plot(land, col = "white", add = TRUE)
# raster::lines(eez, col = "dimgrey", lwd = 0.75)
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
# % fiszh in each 'importance' group, out of those with statistics 
table(spptraits$importance)/table(!is.na(spptraits$importance))[2]*100
23.5364397 + 1.1947431 + 39.1875747


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
# sst_t10, sbt_t10 (MIN)
pos_1 <- which.min(spptraits$sst_t10)
pos_2 <- which.min(spptraits$sbt_t10)
spptraits[pos_1, ]
spptraits[pos_2, ]
raster::plot(raster::raster(paste0("./data/sdm_aquamaps/maps_occ/", spptraits$spp_key_asc[pos_1]))); raster::lines(coastline)
raster::plot(raster::raster(paste0("./data/sdm_aquamaps/maps_occ/", spptraits$spp_key_asc[pos_2]))); raster::lines(coastline)
# sbt_t10 (MAX)
pos_1 <- which.max(spptraits$sst_t10)
pos_2 <- which.max(spptraits$sbt_t10)
spptraits[pos_1, ]
spptraits[pos_2, ]
raster::plot(raster::raster(paste0("./data/sdm_aquamaps/maps_occ/", spptraits$spp_key_asc[pos_1]))); raster::lines(coastline)
raster::plot(raster::raster(paste0("./data/sdm_aquamaps/maps_occ/", spptraits$spp_key_asc[pos_2]))); raster::lines(coastline)
# sst_t90, sbt_t90 (MIN)
pos_1 <- which.min(spptraits$sst_t90)
pos_2 <- which.min(spptraits$sbt_t90)
spptraits[pos_1, ]
spptraits[pos_2, ]
# sst_t90, sbt_t90 (MAX)
pos_1 <- which.max(spptraits$sst_t90)
pos_2 <- which.max(spptraits$sbt_t90)
spptraits[pos_1, ]
spptraits[pos_2, ]
raster::plot(raster::raster(paste0("./data/sdm_aquamaps/maps_occ/", spptraits$spp_key_asc[pos_1]))); raster::lines(coastline)
raster::plot(raster::raster(paste0("./data/sdm_aquamaps/maps_occ/", spptraits$spp_key_asc[pos_2]))); raster::lines(coastline)

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
#### Plot correlations between weighted and unweighted thermal affinities

#### Method
# Thermal affinities have been derived using weighted quantiles, based on the 
# ... probability of species' presence across their distributions. They have 
# ... also been derived from threshold-based distributions, in which absence/presence
# ... has been assigned based on a 0.5 probability threshold and then un-weighted
# ... quantiles have been used to define thermal affinities across areas of species
# ... 'presence'. Here, we check the difference between these two methods in 
# ... the resultant thermal affinities, focusing on SST. 

#### Basic stats
# SST 
utils.add::basic_stats(spptraits$sst_t10 - spptraits$sst_t10_wt)
utils.add::basic_stats(spptraits$sst_t50 - spptraits$sst_t50_wt)
utils.add::basic_stats(spptraits$sst_t90 - spptraits$sst_t90_wt)
# SBT
utils.add::basic_stats(spptraits$sbt_t10 - spptraits$sbt_t10_wt)
utils.add::basic_stats(spptraits$sbt_t50 - spptraits$sbt_t50_wt)
utils.add::basic_stats(spptraits$sbt_t90 - spptraits$sbt_t90_wt)

#### Plots
pp <- par(mfrow = c(2, 2))
pretty_plot(spptraits$sst_t10, spptraits$sst_t10_wt)
pretty_plot(spptraits$sst_t50, spptraits$sst_t10_wt)
pretty_plot(spptraits$sst_t90, spptraits$sst_t10_wt)
par(pp)


##############################
#### Plot correlations between SST and SBT thermal affinities

#### Set up figure
png("./fig/thermal_affinity_correlations.png", 
     height = 6, width = 6, units = "in", res = 600)
pp <- par(mfrow = c(2, 2), oma = c(2.5, 2.5, 2, 2), mar = c(1, 2, 2, 1))
xlim <- ylim <- c(-5, 30)
line_main <- -0.4

#### T10
pretty_plot(spptraits$sbt_t10, spptraits$sst_t10, 
            xlab = "", ylab = "",
            xlim = xlim, ylim = ylim,
            col = scales::alpha("black", 0.75), cex = 0.8)

lines(xlim, ylim)
mtext(side = 3, expression(paste(bold("A "), "[", T[10], "]")),
      adj = 0, font = 2, cex = 1.25, line = line_main)

#### STR
pretty_plot(spptraits$sbt_t50, spptraits$sst_t50, 
            xlab = "", ylab = "",
            xlim = xlim, ylim = ylim,
            col = scales::alpha("black", 0.75), cex = 0.8)
lines(xlim, ylim)
mtext(side = 3, expression(paste(bold("B "), "[STI]")),
      adj = 0, font = 2, cex = 1.25, line = line_main)

#### T90
pretty_plot(spptraits$sbt_t90, spptraits$sst_t90, 
            xlab = "", ylab = "",
            xlim = xlim, ylim = ylim,
            col = scales::alpha("black", 0.75), cex = 0.8)
lines(xlim, ylim)
mtext(side = 3, expression(paste(bold("C "), "[", T[90], "]")),
      adj = 0, font = 2, cex = 1.25, line = line_main)

#### STR
pretty_plot(spptraits$sbt_str, spptraits$sst_str, 
            xlab = "", ylab = "",
            xlim = xlim, ylim = ylim,
            col = scales::alpha("black", 0.75), cex = 0.8)
lines(xlim, ylim)
mtext(side = 3, expression(paste(bold("D "), "[STR]")),
      adj = 0, font = 2, cex = 1.25, line = line_main)

#### Save figure 
mtext(side = 1, expression(paste("SBT (", degree, "C)")), line = 1.5, cex = 1.5, outer = TRUE)
mtext(side = 2, expression(paste("SST (", degree, "C)")), line = 0.25, cex = 1.5, outer = TRUE)
par(pp)
dev.off()


#### End of code. 
##############################
##############################
