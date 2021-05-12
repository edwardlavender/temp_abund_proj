##############################
##############################
#### analyse_abund_across_eezs.R
#### Edward Lavender (el72@st-andrews.ac.uk)

#### This code:
# 1) Analyses relative abundance projections across EEZs;

#### Steps preceding this code:
# 1) Relative abundance predictions (project_abund_* scripts) for modelled species;


##############################
##############################
#### Set up 

#### Wipe workspace
rm(list = ls())

#### Essential packages
library(magrittr)
library(prettyGraphics)

#### Load data 
## EEZ data
eez <- sf::read_sf("./data/spatial/eez", "eez_v11") 
eez <- as(eez, "Spatial")
eez_stats <- readRDS("./data/spatial/eez/eez_stats.rds")
eez_25    <- readRDS("./data/spatial/eez/eez_with_largest_marine_capture_production.rds")
## Species richness map 
spp_richness <- raster::raster("./data/spatial/species/map_spp_richness.asc")
## Abundance predictions across globe 
# SST mean predictions 
ab_delta_sst_mid_rcp45_mean  <- raster::raster("./data/abundance/change/globe/sst/mean/ab_delta_sst_mid_rcp45_mean.asc")
ab_delta_sst_late_rcp45_mean <- raster::raster("./data/abundance/change/globe/sst/mean/ab_delta_sst_late_rcp45_mean.asc")
ab_delta_sst_mid_rcp85_mean  <- raster::raster("./data/abundance/change/globe/sst/mean/ab_delta_sst_mid_rcp85_mean.asc")
ab_delta_sst_late_rcp85_mean <- raster::raster("./data/abundance/change/globe/sst/mean/ab_delta_sst_late_rcp85_mean.asc")
# SBT mean predictions 
ab_delta_sbt_mid_rcp45_mean  <- raster::raster("./data/abundance/change/globe/sbt/mean/ab_delta_sbt_mid_rcp45_mean.asc")
ab_delta_sbt_late_rcp45_mean <- raster::raster("./data/abundance/change/globe/sbt/mean/ab_delta_sbt_late_rcp45_mean.asc")
ab_delta_sbt_mid_rcp85_mean  <- raster::raster("./data/abundance/change/globe/sbt/mean/ab_delta_sbt_mid_rcp85_mean.asc")
ab_delta_sbt_late_rcp85_mean <- raster::raster("./data/abundance/change/globe/sbt/mean/ab_delta_sbt_late_rcp85_mean.asc")
## Abundance predictions by EEZs
ab_sst_in_eez <- readRDS("./data/abundance/change/eez/ab_sst_in_eez.rds")
ab_sbt_in_eez <- readRDS("./data/abundance/change/eez/ab_sbt_in_eez.rds")

#### Helper function to summarise mean predictions across EEZs
summarise_by_eez <- function(x, eez, level_scenario = "mid_rcp45"){
  x_dat <- data.frame(raster::rasterToPoints(x))
  colnames(x_dat) <- c("x", "y", "z")
  sp::coordinates(x_dat) <- c("x", "y")
  raster::crs(x_dat) <- raster::crs(eez)
  x_dat$eez <- sp::over(x_dat, eez)$SOVEREIGN1
  x_dat$eez <- factor(x_dat$eez)
  x_dat <- data.frame(x_dat)
  x_dat <- x_dat[!is.na(x_dat$eez), ]
  x_dat <- 
    x_dat %>%
    dplyr::group_by(eez) %>%
    dplyr::summarise(avg = mean(z),
                     q25 = quantile(z, 0.25), 
                     q75 = quantile(z, 0.75)) %>%
    dplyr::mutate(scenario = level_scenario)
  
  return(x_dat)
}


##############################
##############################
#### Identify relevant EEZs

### Focus on relevant EEZs:
# We will analyse predictions for EEZs, for which 
# ... we have more than a threshold number of species
# ... we have species in more than a threshold number of grid cells 

#### Define the number of species in each EEZ
eez_spp <- sort(table(ab_sst_in_eez$eez))
eez_spp  <- data.frame(eez = names(eez_spp), n_spp = as.integer(eez_spp))
eez_stats$n_spp <- eez_spp$n_spp[match(eez_stats$eez, eez_spp$eez)]
eez_stats <- eez_stats[!is.na(eez_stats$n_spp), ]

#### Analyse EEZ statistics in terms of species and area coverage
## All EEZs 
# Number of spp in EEZs
utils.add::basic_stats(eez_stats$n_spp)
# % of cells in EEZs with modelled fish
utils.add::basic_stats(eez_stats$pc_cell_occ)
# Results:
# The number of species in EEZs ranges from 2 - 1590 (mean = 408.45 spp per EEZ)
# The areal coverage of EEZs (on grid scale) ranges from 1 - 100 % (mean = 28.18 %)
## 25 EEZs
eez_stats_25 <- eez_stats[which(eez_stats$eez %in% eez_25), ]
utils.add::basic_stats(eez_stats_25$n_spp)
utils.add::basic_stats(eez_stats_25$pc_cell_occ)

#### Get EEZ abbreviations and continents
# EEZ abbreviations 
eez_stats$eez_abb <- countrycode::countrycode(eez_stats$eez, 
                                              origin = "country.name", 
                                              destination = "cldr.short.en")
eez_stats$eez_abb[eez_stats$eez == "Micronesia"] <- "Micronesia"
eez_stats$eez_abb[eez_stats$eez == "Comores"] <- "Comores"
# Continents 
eez_stats$continent <- countrycode::countrycode(eez_stats$eez, 
                                              origin = "country.name", 
                                              destination = "continent")
table(eez_stats$continent)
eez_stats$continent[eez_stats$eez == "Micronesia"] <- "Oceania"
eez_stats$continent[eez_stats$eez == "Comores"] <- "Africa"


##############################
##############################
#### Prepare data
# ... Calculate EEZ means and Prs 

##############################
#### SST 

#### For each scenario, get average, mean changes for each EEZ
# This code takes a few minutes to run, so we'll run it once and 
# ... and thereafter load the summaries directly. 
run <- FALSE
if(run){
  t1 <- Sys.time()
  ab_sst_summary_mid_rcp45  <- summarise_by_eez(ab_delta_sst_mid_rcp45_mean, eez, "mid_rcp45")
  ab_sst_summary_late_rcp45 <- summarise_by_eez(ab_delta_sst_late_rcp45_mean, eez, "late_rcp45")
  ab_sst_summary_mid_rcp85  <- summarise_by_eez(ab_delta_sst_mid_rcp85_mean, eez, "mid_rcp85")
  ab_sst_summary_late_rcp85 <- summarise_by_eez(ab_delta_sst_late_rcp85_mean, eez, "late_rcp85")
  t2 <- Sys.time()
  difftime(t2, t1)
  
  #### Join dataframes across scenarios 
  ab_sst_summary <- rbind(ab_sst_summary_mid_rcp45, 
                          ab_sst_summary_late_rcp45, 
                          ab_sst_summary_mid_rcp85, 
                          ab_sst_summary_late_rcp85)
  
  #### Add proportions
  # Calculate, for each species, the overall change in abundance in each EEZ (under SST)
  # ... -ve values are decreases in abundance
  # ... +ve values are increases in abundance 
  ab_sst_in_eez$delta_mid_rcp45  <- ab_sst_in_eez$ab_mid_rcp45  - ab_sst_in_eez$ab_historical
  ab_sst_in_eez$delta_late_rcp45 <- ab_sst_in_eez$ab_late_rcp45 - ab_sst_in_eez$ab_historical
  ab_sst_in_eez$delta_mid_rcp85  <- ab_sst_in_eez$ab_mid_rcp85  - ab_sst_in_eez$ab_historical
  ab_sst_in_eez$delta_late_rcp85 <- ab_sst_in_eez$ab_late_rcp85 - ab_sst_in_eez$ab_historical
  # Calculate the proportion of species expected to decline in each EEZ (under SST)
  ab_sst_summary_pr <- 
    ab_sst_in_eez %>%
    dplyr::group_by(eez) %>%
    dplyr::summarise(
      mid_rcp45 = sum(delta_mid_rcp45 < 0)/dplyr::n(), 
      late_rcp45 = sum(delta_late_rcp45 < 0)/dplyr::n(),
      mid_rcp85 = sum(delta_mid_rcp85 < 0)/dplyr::n(),
      late_rcp85 = sum(delta_late_rcp85 < 0)/dplyr::n()
    ) %>% tidyr::pivot_longer(!eez, names_to = "scenario", values_to = "pr") 
  # Add to dataframe 
  ab_sst_summary <- dplyr::full_join(ab_sst_summary, ab_sst_summary_pr, by = c("eez", "scenario"))
  
  #### Add relevant information 
  # EEZ abbreviations and continent s
  ab_sst_summary$eez_abb   <- eez_stats$eez_abb[match(ab_sst_summary$eez, eez_stats$eez)]
  ab_sst_summary$continent <- eez_stats$continent[match(ab_sst_summary$eez, eez_stats$eez)]
  # Mid point latitude 
  ab_sst_summary$lat_mid  <- abs(eez_stats$lat_mid)[match(ab_sst_summary$eez, eez_stats$eez)]
  # Scenario
  ab_sst_summary$scenario <- factor(ab_sst_summary$scenario, levels = c("mid_rcp45", "mid_rcp85", "late_rcp45", "late_rcp85"))
  # N spp
  ab_sst_summary$n_spp     <- eez_stats$n_spp[match(ab_sst_summary$eez, eez_stats$eez)]
  ab_sst_summary$area      <- round(eez_stats$pc_cell_occ)[match(ab_sst_summary$eez, eez_stats$eez)]
  ab_sst_summary$eez_label <- paste0(ab_sst_summary$eez_abb, " (", ab_sst_summary$n_spp, ")")
  
  #### Define colour scheme
  # Define n cols 
  n_cols <- 1000
  # Define colour palette 
  pal <- viridis::plasma
  # Define colours 
  cols <- rev(pal(n_cols))
  # Define latitudinal bands 
  # (use Absolute latitudes so that the plot is symmetric about the equator)
  breaks <- seq(0, 90, length.out = n_cols)
  data_legend <- data.frame(x = breaks, col = cols)
  data_legend$col <- as.character(data_legend$col)
  # Add colours to dataframe 
  ab_sst_summary$lat_band <- cut(ab_sst_summary$lat_mid, breaks = breaks)
  # Add colours to dataframe 
  ab_sst_summary$col <- cols[unclass(ab_sst_summary$lat_band)]
  # Adjust colour transparency to distinguish scenarios
  ab_sst_summary <- 
    ab_sst_summary %>%
    dplyr::group_by(eez) %>% 
    dplyr::arrange(scenario) %>%
    dplyr::mutate(trans = c(0.25, 0.5, 0.75, 1), 
                  col = scales::alpha(col, trans)) %>%
    dplyr::arrange(continent, lat_mid, scenario)
  
  #### Save
  saveRDS(ab_sst_summary, "./data/abundance/change/eez/ab_sst_summary.rds")

} else{
  
  #### Load file 
  ab_sst_summary <- readRDS("./data/abundance/change/eez/ab_sst_summary.rds")
}

#### Focus on relevant EEZs
# Relevant EEZs should have more than a threshold number of species
# And more than a threshold areal coverage 
head(eez_stats)
length(unique(eez_stats$eez))
eez_stats <- eez_stats[eez_stats$n_spp > 5, ]
# eez_stats <- eez_stats[eez_stats$pc_cell_occ > 10, ]
length(unique(eez_stats$eez))
ab_sst_summary <- ab_sst_summary[ab_sst_summary$eez %in% unique(eez_stats$eez), ]

#### Define convenience plotting order
ab_sst_summary$continent <- factor(ab_sst_summary$continent, levels = c("Asia", "Africa", "Americas", "Oceania", "Europe"))


##############################
#### SBT
# [copy from above, replace 'sst' for 'sbt']



##############################
##############################
#### Continent-specific plots


##############################
#### Helper plotting function(s) 

#### Barplot of changes across EEZ (avg or pr changes)
barplot_across_eezs <- function(data, type = "avg", 
                                dens = NULL,
                                xlim = c(-1, 1), x_grid = 0.5,
                                arrows_col = NULL, arrows_length = 0.01, arrows_lwd = 0.5, 
                                add_y_labels = TRUE, cex.axis = 1
                                ){
  
  #### Focus on specific continent
  # "Africa"   "Americas" "Asia"     "Europe"   "Oceania" 
  # data <- ab_sst_summary[ab_sst_summary$continent == "Oceania", ]
  
  #### Set plotting window and other param
  # pp <- par(oma = c(2, 2, 2, 2))
  cex_axis <- 1 
  
  #### Define x axis param 
  if(is.null(xlim)){
    xlim <- range(c(data$q25, data$q75))
    xlim <- pretty_seq(xlim)$lim
  }

  #### Blank barplot 
  b1 <- barplot(as.numeric(data[, type, drop = TRUE]), 
                axes = FALSE,
                xlim = xlim,
                space = 0,
                horiz = TRUE,
                col = "white", 
                border = "white")
  
  #### Define y axis param 
  n_bar <- length(which(data$eez == data$eez[1]))
  y_shift <- (b1[2] - b1[1])/2
  dat_y_axis <- data.frame(id = rep(seq(1, length(unique(data$eez))), each = n_bar), 
                           at = as.numeric(b1))
  dat_y_axis <- 
    dat_y_axis %>% 
    dplyr::group_by(id) %>% 
    dplyr::summarise(at = mean(at)) %>% 
    data.frame()
  
  #### Add grid 
  # x_at <- pretty_seq(xlim, lim = xlim, pretty_args = list(n = 10))$at
  # x_at <- c(xlim[1], x_at, xlim[2])
  # x_at <- x_at[!duplicated(x_at)]
  x_at <- seq(xlim[1], xlim[2], by = x_grid)
  y_at <- dat_y_axis$at
  y_at <- c(b1[seq(1, max(b1), by = n_bar)]) - y_shift
  y_at <- c(y_at, max(y_at) + (y_at[2] - y_at[1]))
  add_grid_rect_xy(x = x_at, y = y_at, lty = 1)
  
  #### Overlay barplot 
  if(is.null(dens)) dens <- seq(40, by = 20, length.out = n_bar)
  b1 <- barplot(as.numeric(data[, type, drop = TRUE]), 
                axes = FALSE,
                density = dens,
                xlim = xlim,
                space = 0,
                horiz = TRUE, 
                add = TRUE, 
                col = data$col,
                border = "black", 
                lwd = 0.25)
  
  #### Add error bars 
  if(type == "avg"){
    if(is.null(arrows_col)) arrows_col <- scales::alpha("black", rev(seq(0.9, by = -0.1, length.out = n_bar)))
    add_error_bars(x = as.numeric(b1), 
                   fit = data$avg, 
                   lwr = data$q25, 
                   upr = data$q75, 
                   add_fit = NULL,
                   horiz = TRUE,
                   col = arrows_col,
                   length = arrows_length, 
                   lwd = arrows_lwd
    )
  }

  #### Add axes 
  axis(side = 1, x_at, pos = 0, cex.axis = cex.axis)
  axis(side = 2, at = y_at, labels = FALSE, pos = xlim[1])
  if(add_y_labels) axis(side = 2, at = dat_y_axis$at, labels = unique(data$eez_label), lwd.ticks = 0, pos = xlim[1], las = TRUE, cex.axis = cex.axis)
  return(invisible())
}


##############################
#### Average changes in each EEZ 

#### SST
## Set up plot 
png("./fig/proj_abund_by_eez_mean_sst.png", 
     height = 8, width = 12, units = "in", res = 600)
pp <- par(oma = c(2, 9, 1, 2), mar = c(1, 6, 1, 6))
mat <- matrix(c(1, 2, 3, 4,
                1, 2, 3, 5), nrow = 2, byrow = TRUE)
nf <- layout(mat)
# layout.show(nf)

## Make plots 
lapply(1:length(levels(ab_sst_summary$continent)), function(i){
  # i = 1
  cont <- levels(ab_sst_summary$continent)[i]
  ab_sst_summary_for_cont <- ab_sst_summary[ab_sst_summary$continent == cont, ]
  barplot_across_eezs(ab_sst_summary_for_cont, xlim = c(-1, 1), x_grid = 1)
  pl <- par(las = -0.75)
  mtext(side = 3, paste0(LETTERS[i], " (", cont, ")"), font = 2, cex = 1.5, line = -0.75)
  par(pl)
}) %>% invisible()

## Add colour bar legend 
# Define legend col param (copied from above)
n_cols <- 1000
pal <- viridis::plasma
cols <- rev(pal(n_cols))
breaks <- seq(0, 90, length.out = n_cols)
data_legend <- data.frame(x = breaks, col = cols)
data_legend$col <- as.character(data_legend$col)
# Define legend axis axis param 
axis_legend <- prettyGraphics::pretty_axis(side = 4, lim = list(c(0, 90)), pretty = list(), units = list(15))
axis_legend[[1]]$axis$pos <- 1
axis_legend[[1]]$axis$cex.axis <- 1.5
# Add legend 
pn <- par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
TeachingDemos::subplot(add_colour_bar(data_legend = data_legend,
               pretty_axis_args = axis_legend,
               mtext_args = list(side = 4, 
                                 text = expression(paste("Absolute latitude (", degree, ")")), 
                                 line = 4.25, 
                                 cex = 1.5)), 
               x = 2.2, y = 6, size = c(0.08, 5), vadj = 0, hadj = 0)
par(pn)
## Save 
mtext(side = 1, expression(paste("Mean ", E(Delta ~ CRCA), " in each EEZ, ", E(Delta ~ CRCA[EEZ]))), 
      cex = 1.5, line = 1, outer = TRUE)
mtext(side = 2, "EEZ authority", cex = 1.5, line = 4, outer = TRUE)
dev.off()


#### SBT [copy from above, replace 'sst' for 'sbt']
# 


##############################
#### Proportion of species expected to decline

#### SST
## Set up plot 
png("./fig/proj_abund_by_eez_pr_sst.png", 
     height = 8, width = 12, units = "in", res = 600)
mat <- matrix(c(1, 2, 3, 4,
                1, 2, 3, 5), nrow = 2, byrow = TRUE)
layout(mat)
pp <- par(# mfrow = c(2, 3), 
  oma = c(2, 9, 1, 2), 
  mar = c(1, 6.5, 1, 6.5))
## Make plots 
lapply(1:length(levels(ab_sst_summary$continent)), function(i){
  cont <- levels(ab_sst_summary$continent)[i]
  ab_sst_summary_for_cont <- ab_sst_summary[ab_sst_summary$continent == cont, ]
  barplot_across_eezs(ab_sst_summary_for_cont, type = "pr", xlim = c(0, 1))
  pl <- par(las = -0.75)
  mtext(side = 3, paste0(LETTERS[i], " (", cont, ")"), font = 2, cex = 1.4, line = -0.7)
  par(pl)
}) %>% invisible()
## Add colour bar legend 
axis_legend <- prettyGraphics::pretty_axis(side = 4, lim = list(c(0, 90)), pretty = list(), units = list(15))
axis_legend[[1]]$axis$pos <- 1
axis_legend[[1]]$axis$cex.axis <- 1.5

pn <- par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
TeachingDemos::subplot(add_colour_bar(data_legend = data_legend,
                                      pretty_axis_args = axis_legend,
                                      mtext_args = list(side = 4, 
                                                        text = expression(paste("Absolute latitude (", degree, ")")), 
                                                        line = 4, 
                                                        cex = 1.5)), 
                       x = 1.75, y = 8, size = c(0.08, 5), vadj = 0, hadj = 0)
par(pn)
## Save 
mtext(side = 1, expression("Pr. spp. predicted to decline in CRCA in each EEZ," ~ Pr(Delta ~ CRCA[EEZ] < 0)), 
      cex = 1.5, line = 1, outer = TRUE)
mtext(side = 2, "EEZ authority", cex = 1.5, line = 4, outer = TRUE)
dev.off()

#### SBT [copy from above with 'sst' replaced by 'sbt']
# 


##############################
##############################
#### Visualise changes in the 25 EEZs with the largest marine capture production

#### Define the EEZs with the largest marine capture production, for plotting 
# Isolate 25 biggest EEZs
ab_sst_summary_25 <- 
  ab_sst_summary %>% 
  dplyr::filter(eez %in% eez_25 & scenario %in% c("late_rcp45", "late_rcp85")) %>% 
  dplyr::arrange(eez)
# Define EEZ order based on late RCP 4.5
eez_25_order <- 
  ab_sst_summary_25 %>% 
  dplyr::filter(scenario == "late_rcp45") %>% 
  dplyr::arrange(avg) %>% 
  dplyr::pull(eez) 
eez_25_order <- factor(eez_25_order, unique(eez_25_order))
# Order EEZs data
ab_sst_summary_25 <- 
  ab_sst_summary_25 %>% 
  dplyr::mutate(eez = factor(eez, levels = eez_25_order)) %>% 
  dplyr::arrange(desc(eez), scenario)
data = ab_sst_summary_25

#### Set up plot 
png("./fig/proj_abund_by_eez_25_sst.png", 
     height = 10, width = 10, units = "in", res = 600)

#### Set plotting window 
mat <- matrix(c(1, 1, 2), ncol = 3, byrow = T); mat
layout(mat)
pp <- par(oma = c(5, 20, 0, 12), 
          mar = c(2, 2, 2, 2)) 

#### Mean change in CRCA 
barplot_across_eezs(ab_sst_summary_25, 
                    type = "avg", 
                    x_grid = 0.5,
                    dens = c(40, 100), 
                    arrows_col = c(scales::alpha("black", 0.6), "black"),
                    arrows_length = 0.025, arrows_lwd = 1.5, 
                    cex.axis = 1.75)
mtext(side = 1, expression(E(Delta ~ CRCA[EEZ])), cex = 1.5, line = 2.5)
mtext(side = 2, "EEZ authority", cex = 1.5, line = 15, outer = TRUE)

#### Pr species predicted to decline in CRCA 
barplot_across_eezs(ab_sst_summary_25, 
                    type = "pr", xlim = c(0, 1),
                    x_grid = 0.25,
                    dens = c(40, 100), 
                    arrows_col = c(scales::alpha("black", 0.6), "black"),
                    arrows_length = 0.025, arrows_lwd = 1.5, 
                    cex.axis = 1.75,
                    add_y_labels = FALSE)
mtext(side = 1.5, expression(Pr(Delta ~ CRCA[EEZ] < 0)), cex = 1.5, line = 2.5)

#### Legend and labelling
axis_legend[[1]]$axis$cex.axis <- 1.75
TeachingDemos::subplot(add_colour_bar(data_legend = data_legend,
                                      pretty_axis_args = axis_legend,
                                      mtext_args = list(side = 4, 
                                                        text = expression(paste("Absolute latitude (", degree, ")")), 
                                                        line = 6, 
                                                        cex = 1.5)), 
                       x = 1.125, y = 6, size = c(0.08, 6), vadj = 0, hadj = 0)


#### Save 
dev.off()



#### End of code. 
##############################
##############################