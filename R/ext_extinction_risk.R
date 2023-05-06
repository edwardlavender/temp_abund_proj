##############################
##############################
#### ext_extinction_risk.R 
#### Edward Lavender (el72@st-andrews.ac.uk)

#### This code:
# 1) Projects changes in thermal habitat suitability under different climate scenarios
# ... and synthesises, for each species, the total change in the area 
# ... of thermally suitable habitat (weighted by the thermal suitability score) 
# ... between baseline and future scenarios. 
# 2) The code is an extension of the original temp_abund_project
# ... to support MU's extinction-risk work. 
# ... The number of species 'at risk of extinction' under different
# ... thresholds for the % of habitat loss (80, 95, 100 %) is analysed. 

#### Steps preceding this code:
# 1) This code is closely modelled on project_abund_1.R.


##############################
##############################
#### Set up 

#### Wipe workspace
rm(list = ls())
try(pacman::p_unload("all"), silent = TRUE)
dv::clear()

#### Essential packages
library(ggplot2)
library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)
library(tictoc)

#### Global param
root_sst <- "./data/temperature/sst"
root_sbt <- "./data/temperature/sbt"

#### Load data
## spptraits
spptraits <- readRDS("./data/spptraits.rds")
## SST baselines and projections
sst_historical <- raster::raster("./data/temperature/sst/historical/historical.asc")
sst_mid_rcp45  <- raster::raster("./data/temperature/sst/mid_century/rcp45/rcp45.asc")
sst_mid_rcp85  <- raster::raster("./data/temperature/sst/mid_century/rcp85/rcp85.asc")
sst_late_rcp45 <- raster::raster("./data/temperature/sst/late_century/rcp45/rcp45.asc")
sst_late_rcp85 <- raster::raster("./data/temperature/sst/late_century/rcp85/rcp85.asc")
## SBT baselines and projections 
sbt_historical <- raster::raster("./data/temperature/sbt/historical/historical.asc")
sbt_mid_rcp45  <- raster::raster("./data/temperature/sbt/mid_century/rcp45/rcp45.asc")
sbt_mid_rcp85  <- raster::raster("./data/temperature/sbt/mid_century/rcp85/rcp85.asc")
sbt_late_rcp45 <- raster::raster("./data/temperature/sbt/late_century/rcp45/rcp45.asc")
sbt_late_rcp85 <- raster::raster("./data/temperature/sbt/late_century/rcp85/rcp85.asc")

#### Global parameters
overwrite <- FALSE


##############################
##############################
#### Define helper functions 

#### Define helper function to predict species' abundance
# Predicts relative abundance onto a blank raster, x,
# ... across the range of the species' distribution
# ... defined by a species distribution map (sdm), 
# ... using a raster of temperatures (temperature), 
# ... and species-specific thermal parameters (sti, str_sd). 
# ... This assumes that x, sdm and temperature are defined across
# ... the same extent with the same resolution. 
predict_abund <- function(x, sdm, temperature, sti, str_sd, scale = dnorm(sti, mean = sti, sd = str_sd)){
  x[sdm > 0] <- dnorm(temperature[sdm > 0], mean = sti, sd = str_sd)/scale
  return(x)
}

#### Calculate area of suitable habitat, weighted by thermal suitability score 
#' The function calculates the area of suitable habitat, weighted by the suitability, for a given RasterLayer (`prediction`). 
#' @param area A RasterLayer that defines the area of each grid cell
#' @param prediction A RasterLayer that contains a predicted thermal suitability surface (under a historical or future scenario)
calc_suitable_area <- function(area, prediction, plot = FALSE) {
  # Calculate area grid, weighted by the (non-normalised) predictions
  wta <- area * prediction
  # Calculate the total (suitability-weighted) area
  out <- raster::cellStats(wta, "sum")
  if (plot) {
    pp <- par(mfrow = c(1, 3))
    on.exit(par(pp), add = TRUE)
    a_1 <- round(raster::cellStats(area, "sum"))
    a_2 <- round(raster::cellStats(prediction, "sum"))
    # Plot area (unweighted)
    raster::plot(raster::trim(area), main = paste("area", a_1))
    # Plot prediction (sums to one)
    raster::plot(raster::trim(prediction), main = paste("prediction", a_2))
    # Plot area (weighted by prediction)
    raster::plot(raster::trim(wta), main = paste("wta", round(out)))
  }
  out
}

#### Calculate % change in suitable area
#' @param area A RasterLayer (see above).
#' @param historical A RasterLayer that defines the thermal habitat suitability for a historical scenario. 
#' @param projection A RasterLayer that defines the thermal habitat suitability for a future scenario. 
change <- function(area, historical, projection) {
  historical <- calc_suitable_area(area, historical)
  projection <- calc_suitable_area(area, projection)
  ((historical - projection) / historical) * 100
}

#### Identify the % cells that remain 'hospitable' according to some threshold
hospitable <- function(historical, projection, threshold = 0.05, plot = FALSE) {
  historical <- historical >= threshold
  projection <- projection >= threshold
  historical_hospitable  <- raster::cellStats(historical, "sum")
  projection_hospitable  <- raster::cellStats(projection, "sum")
  change <- ((historical_hospitable - projection_hospitable)/historical_hospitable) * 100
  if (plot) {
    pp <- par(mfrow = c(1, 2))
    on.exit(par(pp), add = TRUE)
    raster::plot(raster::trim(historical), 
                 main = paste0("H (", historical_hospitable, " cells)"))
    raster::plot(raster::trim(projection), 
                 main = paste0("P (", projection_hospitable, " cells [", round(change), " %])"))
  }
  change
}

##############################
##############################
#### Generate predictions

file_out <- "./data/extensions/extinction-risk.rds"

if (overwrite | !file.exists(file_out)) {
  
  #### Set up cluster (~2 s)
  tic()
  cl <- parallel::makeCluster(10L)
  vl <- c("predict_abund", "calc_suitable_area", "change", "hospitable",
          "sst_historical",
          "sst_mid_rcp45",
          "sst_mid_rcp85",
          "sst_late_rcp45",
          "sst_late_rcp85",
          "sbt_historical",
          "sbt_mid_rcp45",
          "sbt_mid_rcp85", 
          "sbt_late_rcp45", 
          "sbt_late_rcp85")
  parallel::clusterExport(cl = cl, varlist = vl)
  toc()
  
  ### Loop over each species and make projections (~6 mins)
  tic()
  spptraits$index <- seq_len(nrow(spptraits))
  # spptraits <- spptraits[1:50, ]
  summary_by_species <- pbapply::pblapply(split(spptraits, spptraits$index), cl = cl, function(d){
    
    #### Define param for projections 
    
    ## Load raster for species 
    # d <- spptraits[spptraits$spp == "Acanthemblemaria castroi", ]
    sdm <- raster::raster(paste0("./data/sdm_aquamaps/maps_occ/",  d$spp_key_asc))
    # Define blank raster
    # ... cells where the species occurs will be replaced with abundance predictions
    # ... cells were the species does not occur will remain as NA 
    blank <- sdm
    blank <- raster::setValues(blank, NA)
    # Compute raster area, using terra (more precise than raster::area())
    area  <- raster::raster(terra::cellSize(terra::rast(sdm), mask = TRUE, unit = "km"))
    if (FALSE) {
      pp <- par(mfrow = c(1, 2))
      raster::plot(raster::trim(sdm > 0), col = "black")
      raster::plot(raster::trim(area > 0), col = "black")
      par(pp)
    }
    
    #### SST projections 
    
    ## Define species thermal niche parameters 
    sst_sti    <- d$sst_t50
    sst_str_sd <- (d$sst_t90 - d$sst_t10)/(2*1.281560031) 
    sst_denom  <- dnorm(sst_sti, mean = sst_sti, sd = sst_str_sd)
    
    ## Predict abundance from baseline temps
    ab_sst_historical <- predict_abund(x = blank,  sdm = sdm, temperature = sst_historical, 
                                       sti = sst_sti, str_sd = sst_str_sd, scale = sst_denom)
    if (FALSE) raster::plot(raster::trim(ab_sst_historical))
    
    ## Predict abundance from mid-century scenarios
    ab_sst_mid_rcp45 <- predict_abund(x = blank,  sdm = sdm, temperature = sst_mid_rcp45, 
                                      sti = sst_sti, str_sd = sst_str_sd, scale = sst_denom)
    ab_sst_mid_rcp85 <- predict_abund(x = blank,  sdm = sdm, temperature = sst_mid_rcp85, 
                                      sti = sst_sti, str_sd = sst_str_sd, scale = sst_denom)
    
    ## Predict abundance from late-century scenarios
    ab_sst_late_rcp45 <- predict_abund(x = blank,  sdm = sdm, temperature = sst_late_rcp45, 
                                       sti = sst_sti, str_sd = sst_str_sd, scale = sst_denom)
    ab_sst_late_rcp85 <- predict_abund(x = blank,  sdm = sdm, temperature = sst_late_rcp85, 
                                       sti = sst_sti, str_sd = sst_str_sd, scale = sst_denom)
    if (FALSE) raster::plot(raster::trim(ab_sst_late_rcp85))
    
    #### SBT projections 
    
    ## Define species thermal niche parameters 
    sbt_sti    <- d$sbt_t50
    sbt_str_sd <- (d$sbt_t90 - d$sbt_t10)/(2*1.281560031) 
    sbt_denom  <- dnorm(sbt_sti, mean = sbt_sti, sd = sbt_str_sd)
    
    ## Predict abundance from baseline temps
    ab_sbt_historical <- predict_abund(x = blank,  sdm = sdm, temperature = sbt_historical, 
                                       sti = sbt_sti, str_sd = sbt_str_sd, scale = sbt_denom)
    
    ## Predict abundance from mid-century scenarios
    ab_sbt_mid_rcp45 <- predict_abund(x = blank,  sdm = sdm, temperature = sbt_mid_rcp45, 
                                      sti = sbt_sti, str_sd = sbt_str_sd, scale = sbt_denom)
    ab_sbt_mid_rcp85 <- predict_abund(x = blank,  sdm = sdm, temperature = sbt_mid_rcp85, 
                                      sti = sbt_sti, str_sd = sbt_str_sd, scale = sbt_denom)
    
    ## Predict abundance from late-century scenarios
    ab_sbt_late_rcp45 <- predict_abund(x = blank,  sdm = sdm, temperature = sbt_late_rcp45, 
                                       sti = sbt_sti, str_sd = sbt_str_sd, scale = sbt_denom)
    ab_sbt_late_rcp85 <- predict_abund(x = blank,  sdm = sdm, temperature = sbt_late_rcp85, 
                                       sti = sbt_sti, str_sd = sbt_str_sd, scale = sbt_denom)
    
    
    #### Calculate statistics
    
    ## Example comparison
    if (FALSE) {
      (a_start <- calc_suitable_area(area, ab_sst_historical, plot = T))
      (a_end <- calc_suitable_area(area, ab_sst_late_rcp85, plot = T))
      ((a_start - a_end)/a_start * 100)
      change(area, ab_sst_historical, ab_sst_late_rcp85)
    }
    
    ## Define baselines 
    baselines <- list(sst = ab_sst_historical, 
                      sbt = ab_sbt_historical)
    ## Define projections 
    projections_by_temperature <- list(sst = list(mid_rcp45 = ab_sst_mid_rcp45,
                                                  mid_rcp85 = ab_sst_mid_rcp85, 
                                                  late_rcp45 = ab_sst_late_rcp45, 
                                                  late_rcp85 = ab_sst_late_rcp85), 
                                       sbt = list(mid_rcp45 = ab_sbt_mid_rcp45,
                                                  mid_rcp85 = ab_sbt_mid_rcp85, 
                                                  late_rcp45 = ab_sbt_late_rcp45, 
                                                  late_rcp85 = ab_sbt_late_rcp85)
                                       )
    
    ## Define a dataframe with % change estimates for each temperature/scenario combination
    # Loop over baselines (SST, SBT)
    ncell <- raster::cellStats(sdm > 0, "sum")
    lapply(seq_len(length(baselines)), function(i) {
      # Pull out baseline-specific information
      baseline    <- baselines[[i]]
      med         <- raster::cellStats(baseline, median)
      temperature <- names(baselines)[i]
      projections <- projections_by_temperature[[temperature]]
      # For each projection, calculate the % change in suitable area
      lapply(seq_len(length(projections)), function(j){
        # Define change in suitable area
        est <- change(area, baseline, projections[[j]])
        hos <- hospitable(baseline, projections[[j]])
        stopifnot(length(est) == 1L)
        stopifnot(length(hos) == 1L)
        # Define dataframe with information
        data.frame(species = d$spp, 
                   ncell = ncell,
                   temperature = temperature, 
                   scenario = names(projections)[j], 
                   median = med,
                   change = est, 
                   hospitable = hos)
      }) |> data.table::rbindlist()
      
    }) |> data.table::rbindlist()
    
  })
  
  if (!is.null(cl)) parallel::stopCluster(cl)
  
  # Tidy data frame 
  out <- rbindlist(summary_by_species)
  out <- 
    out |> 
    mutate(species = factor(species), 
           temperature = factor(temperature, 
                                levels = c("sst", "sbt"), 
                                labels = c("SST", "SBT")), 
           scenario = factor(scenario, 
                             levels = c("mid_rcp45", 
                                        "mid_rcp85", 
                                        "late_rcp45", 
                                        "late_rcp85"), 
                             labels = c("RCP 4.5 (mid-century)", 
                                        "RCP 8.5 (mid-century)", 
                                        "RCP 4.5 (late-century)", 
                                        "RCP 8.5 (late-century)")),
           remaining = 100 - change) |> 
    arrange(species, temperature, scenario) |> 
    select(species, ncell, temperature, scenario, median, change, remaining, hospitable) |>
    as.data.table()
  
  # Save outputs 
  beepr::beep(10)
  saveRDS(out, file_out)
  toc()
  
} else {
  out <- readRDS(file_out)
}


##############################
##############################
#### Analyse results

#### Check structure
# * species
# * temperature
# * scenario
# * % change (positive indicates X % decline in area, negative means increase)
# * % remaining is % remaining suitable habitat
# * % hospitable
str(out)

#### Summarise the % change/% remaining habitat
# There is a wide range in responses
# For some species/scenarios, there is 0 % habitat remaining
# On average, there is 80 % habitat remaining
# Some species see increases in habitat
utils.add::basic_stats(out$change)
utils.add::basic_stats(out$remaining)
utils.add::basic_stats(out$hospitable)
# Check ncell and median suitability for species predicted to lose 95 % of 'hospitable' habitat
utils.add::basic_stats(out$ncell[(100 - out$hospitable) < 5])
utils.add::basic_stats(out$median[(100 - out$hospitable) < 5])
# There are limited difference between scenarios (on average)
# ... But order is expected (changes worse under late century, and RCP85)
out |> 
  group_by(temperature, scenario) |> 
  summarise(q = quantile(remaining, 0.05)) |> 
  arrange(temperature, q)

#### Examine the distribution of % remaining habitat for SST/SBT under different scenarios
# % remaining habitat
png("./fig/extensions/remaining_area_suitable.png", 
    height = 8, width = 5, units = "in", res = 600)
ggplot(out) + 
  geom_histogram(aes(remaining), binwidth = 5) + 
  facet_wrap(~scenario * temperature, ncol = 2) + 
  xlab("The amount of thermally suitable habitat remaining (%)") + 
  ylab("Count")
dev.off()
# As above but for % habitat remaining 'hospitable'
png("./fig/extensions/remaining_area_hospitable.png", 
    height = 8, width = 5, units = "in", res = 600)
ggplot(out) + 
  geom_histogram(aes(100 - hospitable), binwidth = 5) + 
  facet_wrap(~scenario + temperature, ncol = 2) + 
  xlab("The amount of thermally 'hospitable' habitat remaining (%)") + 
  ylab("Count")
dev.off()

#### Number of species projected to be left with selected amounts of habitat
# Define thresholds
thresholds <- setNames(c(0, 5, 20), 
                       c("Threshold: 0 %", "Threshold: 5 %", "Threshold: 20 %"))
# Calculate counts
counts <- 
  lapply(thresholds, function(threshold) {
    out |> 
      group_by(temperature, scenario) |> 
      summarise(n = length(which(remaining <= threshold))) |> 
      filter(n > 0) |>
      ungroup() |> 
      mutate(threshold = factor(threshold, 
                                levels = thresholds, 
                                labels = names(thresholds))) |>
      as.data.table() 
}) |> rbindlist()
# Visualise counts
png("./fig/extensions/extinction-risk.png", 
    height = 4, width = 6, units = "in", res = 600)
ggplot(counts) + 
  geom_bar(aes(x = scenario, y = n, fill = scenario), 
           stat = "identity", colour = "black", lwd = 0.10) + 
  scale_fill_viridis_d(alpha = 0.5, option = "viridis", name = "Scenario") +
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + 
  xlab("Scenario") + ylab("Count") + 
  facet_wrap(~temperature + threshold) + 
  theme(axis.text.x = element_blank())
dev.off()

#### Examine example species
# Look at the species with < 5 % remaining habitat
out[, remaining_hospitable := 100 - hospitable]
out[remaining_hospitable < 5, ]

#### Save table to file
head(out)
nrow(out)
length(unique(out$species)) * 2 * 4
write.csv(out, "./data/extensions/extinction-risk.csv", 
          row.names = FALSE)



#### End of code. 
##############################
##############################