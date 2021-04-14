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
eez_cells <- readRDS("./data/spatial/eez/eez_n_cell_tot.rds")
## Species richness map 
spp_richness <- raster::raster("./data/spatial/map_spp_richness.asc")
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

#### Run the code once
# Thereafter, we can simply load in the necessary data. 
run <- FALSE
if(run){
  
  #### Define a starting dataframe with:
  # ... the number of species in each EEZ (only include EEZs represented by model for speed)
  # ... the number of cells with species in (to be defined below)
  # ... the total number of cells in that EE~
  dat_eez <- sort(table(ab_sst$eez))
  dat_eez  <- data.frame(eez = names(dat_eez), n_spp = as.integer(dat_eez), n_cell_occ = NA, n_cell_tot = NA)
  dat_eez$n_cell_tot <- eez_cells$n_cell_tot[match(dat_eez$eez, eez_cells$eez)]
  
  #### Define species richness as SpatialPointsDataFrame
  # This is used to work out the number of cells in each EEZ
  # ... in which we have made species projections 
  dat_spp_richness <- data.frame(raster::rasterToPoints(spp_richness))
  colnames(dat_spp_richness) <- c("x", "y", "n_spp")
  sp::coordinates(dat_spp_richness) <- c("x", "y")
  raster::crs(dat_spp_richness) <- raster::crs(eez)
  
  #### Determine the number of cells with species in in each EEZ
  ## Set up cluster (optional - does not increase speed much)
  cl <- NULL 
  # cl <- parallel::makeCluster(5L)
  if(!is.null(cl)) parallel::clusterExport(cl, c("eez", "dat_spp_richness"))
  # Define a list with the number of cells with species in in each EEZ
  dat_eez_by_eez <- 
    pbapply::pblapply(split(dat_eez, 1:nrow(dat_eez)), cl = cl, function(d){
      # d <- dat_eez[47, ]
      eez_for_state <- subset(eez, SOVEREIGN1 == d$eez)
      # rgdal::writeOGR(eez_for_state, "/Users/el72/Desktop/tmp", d$eez, driver = "ESRI Shapefile")
      # raster::plot(eez_for_state)
      # Number of cells for that EEZ that contain spp 
      tmp <- sp::over(dat_spp_richness, eez_for_state)$SOVEREIGN1
      d$n_cell_occ <- sum(tmp == d$eez, na.rm = TRUE)
      return(d)
    }) 
  if(!is.null(cl)) parallel::stopCluster(cl)
  
  #### Finalise dataframe 
  dat_eez <- do.call(rbind, dat_eez_by_eez)
  dat_eez$pc_cell_occ <- (dat_eez$n_cell_occ/dat_eez$n_cell_tot)*100
  
  #### Analyse 
  # Number of spp in EEZs
  utils.add::basic_stats(dat_eez$n_spp)
  # min   mean median  max     sd   IQR    MAD
  #   2 410.36    232 1744 457.87 727.5 308.38
  # % of cells in EEZs with modelled fish
  utils.add::basic_stats(dat_eez$pc_cell_occ)
  # min  mean median max   sd  IQR MAD
  # 66.67 97.77    100 100 5.21 2.06   0
  
  #### Results
  # The number of species in EEZs ranges from 2 - 1744 (mean = 410.36 spp per EEZ)
  # The areal coverage of EEZs (on 1 dg scale) ranges from 67 - 100 % (mean = 97.77 %)
  saveRDS(dat_eez, "./data/spatial/eez/eez_stats.rds")
  
} else{
  
  #### Load computed file 
  dat_eez <- readRDS("./data/spatial/eez/eez_stats.rds")
  
}

#### Focus on relevant EEZs
# Relevant EEZs should have more than a threshold number of species
# And more than a threshold areal coverage 
# ... But since areal coverage is greater than 67 % in all cases, we'll ignore this criterion. 
head(dat_eez)
dat_eez <- dat_eez[dat_eez$n_spp > 25, ]
length(unique(dat_eez$eez))

#### Get EEZ abbreviations and continents
# EEZ abbreviations 
dat_eez$eez_abb <- countrycode::countrycode(dat_eez$eez, 
                                              origin = "country.name", 
                                              destination = "cldr.short.en")
dat_eez$eez_abb[dat_eez$eez == "Micronesia"] <- "Micronesia"
dat_eez$eez_abb[dat_eez$eez == "Comores"] <- "Comores"
# Continents 
dat_eez$continent <- countrycode::countrycode(dat_eez$eez, 
                                              origin = "country.name", 
                                              destination = "continent")
table(dat_eez$continent)
dat_eez$continent[dat_eez$eez == "Micronesia"] <- "Oceania"
dat_eez$continent[dat_eez$eez == "Comores"] <- "Africa"


##############################
##############################
#### Prepare data
# ... Calculate EEZ means and Prs 

##############################
#### SST 

#### For each scenario, get average, mean changes for each EEZ
# This code takes ~ 2.5 minutes to run, so we'll run it once and 
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
  
  #### Save
  saveRDS(ab_sst_summary, "./data/abundance/change/eez/ab_sst_summary.rds")
  
} else{
  
  #### Load file 
  ab_sst_summary <- readRDS("./data/abundance/change/eez/ab_sst_summary.rds")
}

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

#### Focus on relevant EEZs
ab_sst_summary <- ab_sst_summary[ab_sst_summary$eez %in% unique(dat_eez$eez), ]

#### Add relevant information 
# EEZ abbreviations and continent s
ab_sst_summary$eez_abb   <- dat_eez$eez_abb[match(ab_sst_summary$eez, dat_eez$eez)]
ab_sst_summary$continent <- dat_eez$continent[match(ab_sst_summary$eez, dat_eez$eez)]
# Mid point latitude 
ab_sst_summary$lat_mid  <- abs(eez_cells$lat_mid)[match(ab_sst_summary$eez, eez_cells$eez)]
# Scenario
ab_sst_summary$scenario <- factor(ab_sst_summary$scenario, levels = c("mid_rcp45", "mid_rcp85", "late_rcp45", "late_rcp85"))
# N spp
ab_sst_summary$n_spp     <- dat_eez$n_spp[match(ab_sst_summary$eez, dat_eez$eez)]
ab_sst_summary$area      <- round(dat_eez$pc_cell_occ)[match(ab_sst_summary$eez, dat_eez$eez)]
ab_sst_summary$eez_label <- paste0(ab_sst_summary$eez_abb, " (", ab_sst_summary$n_spp, ")")

#### Define colour scheme
# Define n cols 
n_cols <- 1000
# Define colour palette 
pal <- viridis::plasma
# Define colours 
cols <- rev(pal(n_cols))
# Define latitudinal bands 
# (use absolute latitudes so that the plot is symmetric about the equator)
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


##############################
#### SBT
# [copied from above, replacing 'sst' for 'sbt']

#### For each scenario, get average, mean changes for each EEZ
# This code takes ~ 5 minutes to run, so we'll run it once and 
# ... and thereafter load the summaries directly. 
run <- FALSE
if(run){
  t1 <- Sys.time()
  ab_sbt_summary_mid_rcp45  <- summarise_by_eez(ab_delta_sbt_mid_rcp45_mean, eez, "mid_rcp45")
  ab_sbt_summary_late_rcp45 <- summarise_by_eez(ab_delta_sbt_late_rcp45_mean, eez, "late_rcp45")
  ab_sbt_summary_mid_rcp85  <- summarise_by_eez(ab_delta_sbt_mid_rcp85_mean, eez, "mid_rcp85")
  ab_sbt_summary_late_rcp85 <- summarise_by_eez(ab_delta_sbt_late_rcp85_mean, eez, "late_rcp85")
  t2 <- Sys.time()
  difftime(t2, t1)
  
  #### Join dataframes across scenarios 
  ab_sbt_summary <- rbind(ab_sbt_summary_mid_rcp45, 
                          ab_sbt_summary_late_rcp45, 
                          ab_sbt_summary_mid_rcp85, 
                          ab_sbt_summary_late_rcp85)
  
  #### Save
  saveRDS(ab_sbt_summary, "./data/abundance/change/eez/ab_sbt_summary.rds")
  
} else{
  
  #### Load file 
  ab_sbt_summary <- readRDS("./data/abundance/change/eez/ab_sbt_summary.rds")
}

#### Add proportions
# Calculate, for each species, the overall change in abundance in each EEZ (under sbt)
# ... -ve values are decreases in abundance
# ... +ve values are increases in abundance 
ab_sbt_in_eez$delta_mid_rcp45  <- ab_sbt_in_eez$ab_mid_rcp45  - ab_sbt_in_eez$ab_historical
ab_sbt_in_eez$delta_late_rcp45 <- ab_sbt_in_eez$ab_late_rcp45 - ab_sbt_in_eez$ab_historical
ab_sbt_in_eez$delta_mid_rcp85  <- ab_sbt_in_eez$ab_mid_rcp85  - ab_sbt_in_eez$ab_historical
ab_sbt_in_eez$delta_late_rcp85 <- ab_sbt_in_eez$ab_late_rcp85 - ab_sbt_in_eez$ab_historical
# Calculate the proportion of species expected to decline in each EEZ (under sbt)
ab_sbt_summary_pr <- 
  ab_sbt_in_eez %>%
  dplyr::group_by(eez) %>%
  dplyr::summarise(
    mid_rcp45 = sum(delta_mid_rcp45 < 0)/dplyr::n(), 
    late_rcp45 = sum(delta_late_rcp45 < 0)/dplyr::n(),
    mid_rcp85 = sum(delta_mid_rcp85 < 0)/dplyr::n(),
    late_rcp85 = sum(delta_late_rcp85 < 0)/dplyr::n()
  ) %>% tidyr::pivot_longer(!eez, names_to = "scenario", values_to = "pr") 
# Add to dataframe 
ab_sbt_summary <- dplyr::full_join(ab_sbt_summary, ab_sbt_summary_pr, by = c("eez", "scenario"))

#### Focus on relevant EEZs
ab_sbt_summary <- ab_sbt_summary[ab_sbt_summary$eez %in% unique(dat_eez$eez), ]

#### Add relevant information 
# EEZ abbreviations and continent s
ab_sbt_summary$eez_abb   <- dat_eez$eez_abb[match(ab_sbt_summary$eez, dat_eez$eez)]
ab_sbt_summary$continent <- dat_eez$continent[match(ab_sbt_summary$eez, dat_eez$eez)]
# Mid point latitude 
ab_sbt_summary$lat_mid  <- abs(eez_cells$lat_mid)[match(ab_sbt_summary$eez, eez_cells$eez)]
# Scenario
ab_sbt_summary$scenario <- factor(ab_sbt_summary$scenario, levels = c("mid_rcp45", "mid_rcp85", "late_rcp45", "late_rcp85"))
# N spp
ab_sbt_summary$n_spp     <- dat_eez$n_spp[match(ab_sbt_summary$eez, dat_eez$eez)]
ab_sbt_summary$area      <- round(dat_eez$pc_cell_occ)[match(ab_sbt_summary$eez, dat_eez$eez)]
ab_sbt_summary$eez_label <- paste0(ab_sbt_summary$eez_abb, " (", ab_sbt_summary$n_spp, ")")

#### Define colour scheme
# Define n cols 
n_cols <- 1000
# Define colour palette 
pal <- viridis::plasma
# Define colours 
cols <- rev(pal(n_cols))
# Define latitudinal bands 
# (use absolute latitudes so that the plot is symmetric about the equator)
breaks <- seq(0, 90, length.out = n_cols)
data_legend <- data.frame(x = breaks, col = cols)
data_legend$col <- as.character(data_legend$col)
# Add colours to dataframe 
ab_sbt_summary$lat_band <- cut(ab_sbt_summary$lat_mid, breaks = breaks)
# Add colours to dataframe 
ab_sbt_summary$col <- cols[unclass(ab_sbt_summary$lat_band)]
# Adjust colour transparency to distinguish scenarios
ab_sbt_summary <- 
  ab_sbt_summary %>%
  dplyr::group_by(eez) %>% 
  dplyr::arrange(scenario) %>%
  dplyr::mutate(trans = c(0.25, 0.5, 0.75, 1), 
                col = scales::alpha(col, trans)) %>%
  dplyr::arrange(continent, lat_mid, scenario)


##############################
##############################
#### Continent-specific plots


##############################
#### Helper plotting function(s) 

#### Barplot of changes across EEZ (avg or pr changes)
barplot_across_eezs <- function(data, type = "avg", xlim = c(-1, 1)){
  
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
  y_shift <- (b1[2] - b1[1])/2
  dat_y_axis <- data.frame(id = rep(seq(1, length(unique(data$eez))), each = 4), 
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
  x_at <- seq(-1, 1, by = 0.5)
  y_at <- dat_y_axis$at
  y_at <- c(b1[seq(1, max(b1), by = 4)]) - y_shift
  y_at <- c(y_at, max(y_at) + (y_at[2] - y_at[1]))
  add_grid_rect_xy(x = x_at, y = y_at, lty = 1)
  
  #### Overlay barplot 
  b1 <- barplot(as.numeric(data[, type, drop = TRUE]), 
                axes = FALSE,
                xlim = xlim,
                space = 0,
                horiz = TRUE, 
                add = TRUE, 
                col = data$col,
                border = "black", 
                lwd = 0.25)
  
  #### Add error bars 
  if(type == "avg"){
    add_error_bars(x = as.numeric(b1), 
                   fit = data$avg, 
                   lwr = data$q25, 
                   upr = data$q75, 
                   add_fit = NULL,
                   horiz = TRUE,
                   col = scales::alpha("black", c(0.6, 0.7, 0.8, 0.9)),
                   length = 0.01, 
                   lwd = 0.5
    )
  }

  #### Add axes 
  # axis(side = 1, x_at, lwd.ticks = 0, labels = FALSE, pos = 0)
  # x_at <- x_at[-c(1, length(x_at))]
  if(type == "avg") x_at <- c(-1, 0, 1) else x_at <- c(0, 0.5, 1)
  axis(side = 1, x_at, pos = 0)
  axis(side = 2, at = y_at, labels = FALSE, pos = xlim[1])
  axis(side = 2, at = dat_y_axis$at, labels = unique(data$eez_label), lwd.ticks = 0, pos = xlim[1], las = TRUE)
  return(invisible())
}


##############################
#### Average changes in each EEZ 

#### SST
## Set up plot 
tiff("./fig/proj_abund_by_eez_mean_sst.tif", 
     height = 8, width = 12, units = "in", res = 600)
pp <- par(oma = c(2, 9, 1, 2), mar = c(1, 6, 1, 6))
mat <- matrix(c(1, 2, 3, 4,
                1, 2, 3, 5), nrow = 2, byrow = TRUE)
nf <- layout(mat)
# layout.show(nf)

## Make plots 
lapply(1:length(unique(ab_sst_summary$continent)), function(i){
  # i = 1
  cont <- unique(ab_sst_summary$continent)[i]
  ab_sst_summary_for_cont <- ab_sst_summary[ab_sst_summary$continent == cont, ]
  barplot_across_eezs(ab_sst_summary_for_cont, xlim = c(-1, 1))
  pl <- par(las = -0.75)
  mtext(side = 3, paste0(LETTERS[i], " (", cont, ")"), font = 2, cex = 1.5, line = -0.75)
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
                                 text = expression(paste("Absolute Latitude (", degree, ")")), 
                                 line = 4, 
                                 cex = 1.5)), 
               x = 2.2, y = 8, size = c(0.08, 5), vadj = 0, hadj = 0)
par(pn)
## Save 
mtext(side = 1, expression(mu[E]), cex = 1.5, line = 0, outer = TRUE)
mtext(side = 2, "EEZ Authority", cex = 1.5, line = 4, outer = TRUE)
dev.off()


#### SBT [copied from above, replacing 'sst' for 'sbt']

#### SST
## Set up plot 
tiff("./fig/proj_abund_by_eez_mean_sbt.tif", 
     height = 8, width = 12, units = "in", res = 600)
pp <- par(oma = c(2, 9, 1, 2), mar = c(1, 6, 1, 6))
mat <- matrix(c(1, 2, 3, 4,
                1, 2, 3, 5), nrow = 2, byrow = TRUE)
nf <- layout(mat)
# layout.show(nf)

## Make plots 
lapply(1:length(unique(ab_sbt_summary$continent)), function(i){
  # i = 1
  cont <- unique(ab_sbt_summary$continent)[i]
  ab_sbt_summary_for_cont <- ab_sbt_summary[ab_sbt_summary$continent == cont, ]
  barplot_across_eezs(ab_sbt_summary_for_cont, xlim = c(-1, 1))
  pl <- par(las = -0.75)
  mtext(side = 3, paste0(LETTERS[i], " (", cont, ")"), font = 2, cex = 1.5, line = -0.75)
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
                                                        text = expression(paste("Absolute Latitude (", degree, ")")), 
                                                        line = 4, 
                                                        cex = 1.5)), 
                       x = 2.2, y = 8, size = c(0.08, 5), vadj = 0, hadj = 0)
par(pn)
## Save 
mtext(side = 1, expression(mu[E]), cex = 1.5, line = 0, outer = TRUE)
mtext(side = 2, "EEZ Authority", cex = 1.5, line = 4, outer = TRUE)
dev.off()


##############################
#### Proportion of species expected to decline

#### SST
## Set up plot 
tiff("./fig/proj_abund_by_eez_pr_sst.tif", 
     height = 8, width = 12, units = "in", res = 600)
mat <- matrix(c(1, 2, 3, 4,
                1, 2, 3, 5), nrow = 2, byrow = TRUE)
layout(mat)
pp <- par(# mfrow = c(2, 3), 
  oma = c(2, 9, 1, 2), 
  mar = c(1, 6.5, 1, 6.5))
## Make plots 
lapply(1:length(unique(ab_sst_summary$continent)), function(i){
  cont <- unique(ab_sst_summary$continent)[i]
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
                                                        text = expression(paste("Absolute Latitude (", degree, ")")), 
                                                        line = 4, 
                                                        cex = 1.5)), 
                       x = 1.75, y = 8, size = c(0.08, 5), vadj = 0, hadj = 0)
par(pn)
## Save 
mtext(side = 1, expression(P[E]), cex = 1.5, line = 0, outer = TRUE)
mtext(side = 2, "EEZ Authority", cex = 1.5, line = 4, outer = TRUE)
dev.off()

#### SBT [copied from above with 'sst' replaced by 'sbt']
## Set up plot 
tiff("./fig/proj_abund_by_eez_pr_sbt.tif", 
     height = 8, width = 12, units = "in", res = 600)
mat <- matrix(c(1, 2, 3, 4,
                1, 2, 3, 5), nrow = 2, byrow = TRUE)
layout(mat)
pp <- par(# mfrow = c(2, 3), 
  oma = c(2, 9, 1, 2), 
  mar = c(1, 6.5, 1, 6.5))
## Make plots 
lapply(1:length(unique(ab_sbt_summary$continent)), function(i){
  cont <- unique(ab_sbt_summary$continent)[i]
  ab_sbt_summary_for_cont <- ab_sbt_summary[ab_sbt_summary$continent == cont, ]
  barplot_across_eezs(ab_sbt_summary_for_cont, type = "pr", xlim = c(0, 1))
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
                                                        text = expression(paste("Absolute Latitude (", degree, ")")), 
                                                        line = 4, 
                                                        cex = 1.5)), 
                       x = 1.75, y = 8, size = c(0.08, 5), vadj = 0, hadj = 0)
par(pn)
## Save 
mtext(side = 1, expression(P[E]), cex = 1.5, line = 0, outer = TRUE)
mtext(side = 2, "EEZ Authority", cex = 1.5, line = 4, outer = TRUE)
dev.off()



#### End of code. 
##############################
##############################