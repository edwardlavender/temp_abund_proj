##############################
##############################
#### analyse_scenarios.R 
#### Edward Lavender (el72@st-andrews.ac.uk)

#### This code:
# 1) Compares the severity of modelled scenarios
# ... mid-century and late-century
# ... RCP 4.5 and RCP 8.5 

#### Steps preceding this code:
# 1) Generate and synthesise predictions 
# ... project_abund_1.R
# ... project_abund_2.R


##############################
##############################
#### Set up 

#### Wipe workspace
rm(list = ls())

#### Essential packages
library(magrittr)
library(prettyGraphics)

#### Global param
root_ab_delta_sst <- "./data/abundance/change/globe/sst/"
root_ab_delta_sbt <- "./data/abundance/change/globe/sbt/"
source("./R/helpers.R")

#### Load data 
# Global grids 
sst_mid_rcp45_mean  <- load_projections(root_ab_delta_sst, "sst", "mid_rcp45", "mean")
sst_late_rcp45_mean <- load_projections(root_ab_delta_sst, "sst", "late_rcp45", "mean")
sst_mid_rcp85_mean  <- load_projections(root_ab_delta_sst, "sst", "mid_rcp85", "mean")
sst_late_rcp85_mean <- load_projections(root_ab_delta_sst, "sst", "late_rcp85", "mean")
# EEZ summaries
eez_stats      <- readRDS("./data/spatial/eez/eez_stats.rds")
ab_sst_in_eez  <- readRDS("./data/abundance/change/eez/ab_sst_in_eez.rds")
ab_sst_summary <- readRDS("./data/abundance/change/eez/ab_sst_summary.rds")


##############################
##############################
#### Data processing 

#### Define the number of species in each EEZ
eez_spp <- sort(table(ab_sst_in_eez$eez))
eez_spp  <- data.frame(eez = names(eez_spp), n_spp = as.integer(eez_spp))
eez_stats$n_spp <- eez_spp$n_spp[match(eez_stats$eez, eez_spp$eez)]
eez_stats <- eez_stats[!is.na(eez_stats$n_spp), ]

#### Focus on relevant EEZs
# Relevant EEZs should have more than a threshold number of species
# And more than a threshold areal coverage 
head(eez_stats)
length(unique(eez_stats$eez))
eez_stats <- eez_stats[eez_stats$n_spp > 5, ]
# eez_stats <- eez_stats[eez_stats$pc_cell_occ > 10, ]
length(unique(eez_stats$eez))
ab_sst_summary <- ab_sst_summary[ab_sst_summary$eez %in% unique(eez_stats$eez), ]


##############################
##############################
#### Mid-century 

#### Global grid
sst_mid_delta <- sst_mid_rcp85_mean - sst_mid_rcp45_mean
pp <- par(mfrow = c(2, 2))
raster::plot(sst_mid_rcp45_mean)
raster::plot(sst_mid_rcp85_mean)
raster::plot(sst_mid_delta)
raster::plot(abs(sst_mid_delta))
par(pp)
# In the mid century, the two scenarios are relatively similar
# ... Declines are marginally are greater in the tropics 
# ... Increases in abundance are greater in temperate regions 
# The absolute differences between the scenarios are small:
utils.add::basic_stats(abs(sst_mid_delta[]), na.rm = TRUE)
quantile(abs(sst_mid_delta[]), probs = c(0.25, 0.5, 0.75), na.rm = TRUE)

#### Across EEZs 
ab_sst_summary_avg <- 
  ab_sst_summary %>% dplyr::select(scenario, eez, avg, pr) %>% 
  tidyr::pivot_wider(names_from = scenario, values_from = c(avg, pr)) %>%
  dplyr::mutate(avg_mid_delta  = avg_mid_rcp85 - avg_mid_rcp45, 
                avg_late_delta = avg_late_rcp85 - avg_late_rcp45)
ab_sst_summary_avg$col   <- ab_sst_summary$col[match(ab_sst_summary_avg$eez, ab_sst_summary$eez)]
ab_sst_summary_avg$n_spp <- ab_sst_summary$n_spp[match(ab_sst_summary_avg$eez, ab_sst_summary$eez)]
utils.add::basic_stats(abs(ab_sst_summary_avg$avg_mid_delta))


##############################
##############################
#### Late-century 

#### Global grid
sst_late_delta <- sst_late_rcp85_mean - sst_late_rcp45_mean
pp <- par(mfrow = c(2, 2))
raster::plot(sst_late_rcp45_mean)
raster::plot(sst_late_rcp85_mean)
raster::plot(sst_late_delta)
raster::plot(abs(sst_late_delta))
par(pp)
# The differences between the scenarios are much more intense in the tropics
utils.add::basic_stats(abs(sst_late_delta[]), na.rm = TRUE)

#### Across EEZs 
utils.add::basic_stats(abs(ab_sst_summary_avg$avg_late_delta))


##############################
##############################
#### Visual comparison (late-century)

### Set up plot 
tiff("./fig/compare_scenarios.tiff", 
     height = 4, width = 8, units = "in", res = 600)
pp <- par(mfrow = c(1, 2), 
          oma = c(2, 3, 2, 5), 
          mar = c(2, 2, 2, 2))

#### Mean changes

## Define colour scheme
ab_sst_summary_avg$lat_mid <- eez_stats$lat_mid_abs[match(ab_sst_summary_avg$eez, eez_stats$eez)]
n_cols <- 1000
pal <- viridis::plasma
cols <- rev(pal(n_cols))
breaks <- seq(0, 90, length.out = n_cols)
data_legend <- data.frame(x = breaks, col = cols)
data_legend$col <- as.character(data_legend$col)
ab_sst_summary_avg$lat_band <- cut(ab_sst_summary_avg$lat_mid, breaks = breaks)
ab_sst_summary_avg$col <- cols[unclass(ab_sst_summary_avg$lat_band)]

## Make plot 
xlim <- ylim <- c(-1, 1)
pretty_plot(ab_sst_summary_avg$avg_late_rcp45, ab_sst_summary_avg$avg_late_rcp85, 
            xlim = xlim, ylim = ylim,
            xlab = "", ylab = "",
            pch = 21,
            bg = ab_sst_summary_avg$col, col = ab_sst_summary_avg$col, 
            cex = ab_sst_summary_avg$n_spp/750)

# Add helper lines
lines(xlim, ylim, lwd = 0.75, lty = 2)
lines(c(0, 0), ylim, lwd = 0.75, lty = 3)
lines(xlim, c(0, 0), lwd = 0.75, lty = 3)

## Add linear model 
# Fit linear model using weighted regression
m1 <- MASS::rlm(avg_late_rcp85 ~ avg_late_rcp45, 
                weights = ab_sst_summary_avg$n_spp, 
                data = ab_sst_summary_avg)
summary(m1)
# Add predictions across the middle 95 % of data 
x <- seq(quantile(ab_sst_summary_avg$avg_late_rcp45, 0.025), 
         quantile(ab_sst_summary_avg$avg_late_rcp45, 0.975), 
         length.out = 100)
preds <- predict(m1, newdata = data.frame(avg_late_rcp45 = x), type = "response", se.fit = TRUE)
preds <- list_CIs(preds)
add_error_envelope(x, preds, 
                   add_fit = list(lwd = 0.75))

## Add titles
mtext(side = 1, expression(paste(E(Delta ~ IRA[EEZ]), " in RCP 4.5 [L]")), line = 2.5, cex = 1.25)
mtext(side = 2, expression(paste(E(Delta ~ IRA[EEZ]), " in RCP 8.5 [L]")), line = 3, cex = 1.25)
mtext(side = 3, font = 2, "A", line = 0.5, adj = 0, cex = 1.25)


#### Proportion of species expected to decline 

## Plot density for rcp4.5 and overlay density for rcp8.5
col_rcp45 <- "grey75"
col_rcp85 <- "grey50"
dens_rcp45 <- density(ab_sst_summary_avg$pr_late_rcp45, from = 0, to = 1)
dens_rcp85 <- density(ab_sst_summary_avg$pr_late_rcp85, from = 0, to = 1)
x <- c(0, 1) 
y <- range(c(dens_rcp45$y, dens_rcp85$y))
axis_ls <- pretty_plot(dens_rcp45,
                       pretty_axis_args = list(x = list(x = x, y = y)),
                       xlab = "", ylab = "",
                       col = col_rcp45, 
                       type = "l", lwd = 2)
lines(dens_rcp85, col = col_rcp85, lwd = 2)

## Add rugs for observed data 
# Use loop because rug isn't vectorised over colour 
ylim <- axis_ls[[2]]$lim
for(i in 1:nrow(ab_sst_summary_avg)){
  rug(ab_sst_summary_avg$pr_late_rcp85[i], 
      pos = ylim[1]+0.01, ticksize = 0.03, col = ab_sst_summary_avg$col[i], lwd = 1.5) 
  rug(ab_sst_summary_avg$pr_late_rcp45[i], 
      pos = ylim[1]+0.1, ticksize = 0.03, col = ab_sst_summary_avg$col[i], lwd = 1.5)
}

## Add legend 
legend(0.01, ylim[2] * 1.075,
       fill = c(col_rcp45, col_rcp85), 
       border =  c(col_rcp45, col_rcp85), 
       legend = c("RCP 4.5 [L]", "RCP 8.5 [L]"), 
       bty = "n", 
       y.intersp = 1.2, 
       cex = 1)

## Add titles
mtext(side = 1, expression(Pr(Delta ~ IRA[EEZ] < 0)), line = 2.5, cex = 1.25)
mtext(side = 2, "Density", cex = 1.25, line = 2.5)
mtext(side = 3, font = 2, "B", cex = 1.25, adj = 0, line = 0.5)

## Add colour bar legend
axis_legend <- prettyGraphics::pretty_axis(side = 4, lim = list(c(0, 90)), pretty = list(), units = list(15))
axis_legend[[1]]$axis$pos <- 1
axis_legend[[1]]$axis$cex.axis <- 1.25

TeachingDemos::subplot(prettyGraphics::add_colour_bar(data_legend = data_legend,
                                                      pretty_axis_args = axis_legend,
                                                      mtext_args = list(side = 4, 
                                                                        text = expression(paste("Absolute Latitude (", degree, ")")), 
                                                                        line = 2.8, 
                                                                        cex = 1.25
                                                      )
),
x = 1.1, y = 0.1, size = c(0.08, 2), vadj = 0, hadj = 0
)

#### Save
dev.off()


#### Compute summary statistics
# Less than 20 % decline 
table(ab_sst_summary_avg$pr_late_rcp45 <= 0.2)
(table(ab_sst_summary_avg$pr_late_rcp45 <= 0.2)/nrow(ab_sst_summary_avg))*100
# Less than 50 % decline 
table(ab_sst_summary_avg$pr_late_rcp45 <= 0.5)
(table(ab_sst_summary_avg$pr_late_rcp45 <= 0.5)/nrow(ab_sst_summary_avg))*100
# More than 80 % decline 
table(ab_sst_summary_avg$pr_late_rcp45 >= 0.8)
(table(ab_sst_summary_avg$pr_late_rcp45 >= 0.8)/nrow(ab_sst_summary_avg))*100


#### End of code. 
##############################
##############################