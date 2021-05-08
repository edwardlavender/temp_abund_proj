##############################
##############################
#### analyse_abund_across_globe_2.R 
#### Edward Lavender (el72@st-andrews.ac.uk)

#### This code:
# 1) Maps abundance predictions over space. 

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
# 

#### Load data 
spptraits <- readRDS("./data/spptraits.rds")
spptraits$index <- 1:nrow(spptraits)
# spptraits <- spptraits[1:10, ]

#### Global param 
root_ab_delta_sst <- "./data/abundance/change/globe/sst/"
root_ab_delta_sbt <- "./data/abundance/change/globe/sbt/"
source("./R/helpers.R")

#### Load files
## Coastline
coastline <- readRDS("./data/spatial/coastline/coastline.rds")
land <- flapper::invert_poly(coastline)
cover <- raster::raster("./data/spatial/eez/eez_mask.asc")
# cover <- NULL
## Species richness
spp_richness <- raster::raster("./data/spatial/species/map_spp_richness.asc")
## SST (mean) 
sst_mid_rcp45_mean  <- load_projections(root_ab_delta_sst, "sst", "mid_rcp45", "mean", mask = cover)
sst_late_rcp45_mean <- load_projections(root_ab_delta_sst, "sst", "late_rcp45", "mean", mask = cover)
sst_mid_rcp85_mean  <- load_projections(root_ab_delta_sst, "sst", "mid_rcp85", "mean", mask = cover)
sst_late_rcp85_mean <- load_projections(root_ab_delta_sst, "sst", "late_rcp85", "mean", mask = cover)
## SST (pr)
sst_mid_rcp45_pr  <- load_projections(root_ab_delta_sst, "sst", "mid_rcp45", "pr", mask = cover)
sst_late_rcp45_pr <- load_projections(root_ab_delta_sst, "sst", "late_rcp45", "pr", mask = cover)
sst_mid_rcp85_pr  <- load_projections(root_ab_delta_sst, "sst", "mid_rcp85", "pr", mask = cover)
sst_late_rcp85_pr <- load_projections(root_ab_delta_sst, "sst", "late_rcp85", "pr", mask = cover)
## SBT (mean)
sbt_mid_rcp45_mean  <- load_projections(root_ab_delta_sbt, "sbt", "mid_rcp45", "mean", mask = cover)
sbt_late_rcp45_mean <- load_projections(root_ab_delta_sbt, "sbt", "late_rcp45", "mean", mask = cover)
sbt_mid_rcp85_mean  <- load_projections(root_ab_delta_sbt, "sbt", "mid_rcp85", "mean", mask = cover)
sbt_late_rcp85_mean <- load_projections(root_ab_delta_sbt, "sbt", "late_rcp85", "mean", mask = cover)


##############################
##############################
#### Helpers

#### Mask projections in areas with fewer than a threshold number of species
# spp_mask <- spp_richness >= 25
# spp_mask[spp_mask == 0] <- NA
# ... We will use the EEZ mask instead. 

#### Wrapper function to plot projections
plot_projections <- function(x, zlim = NULL, type = 1, mask = NULL){
  # Define graphical param 
  line_main <- -0.9
  adj_main  <- 0.05
  cex_main  <- 1.5
  # Define wrapper function for plotting with appropriate colour scheme
  # ... predictions are clearly presented with a split colour scheme
  # ... SDs (which can't be negative) are better presented with a contiguous colour scheme
  if(type == 1){
    plot_x <- function(x)
      plot_raster(x, mask = mask, 
                  zlim = zlim, add_legend = TRUE,
                  gen_cols = prettyGraphics::pretty_cols_split_heat, 
                  scheme_cold = "YlOrRd", scheme_hot = "Greens", 
                  select_cold = 3:8, select_hot = 4:8,
                  profile_x = c(185, 210), profile_y = c(-90, 90))
    } else if(type == 2){
      plot_x <- function(x)
        plot_raster(x, mask = mask,
                    zlim = zlim, add_legend = TRUE,
                    # gen_cols = prettyGraphics::pretty_cols_brewer, pal = viridis::plasma, 
                    gen_cols = prettyGraphics::pretty_cols_split_heat, 
                    split = 0.5,
                    scheme_cold = "Greens", scheme_hot = "YlOrRd", 
                    select_cold = 2:8, select_hot = 2:8,
                    profile_x = c(185, 210))
    } else stop("'type' not supported.")
  # mid_rcp45_mean
  col_param <- plot_x(x[[1]])
  mtext(side = 3, "A (RCP 4.5 [M])", adj = adj_main, line = line_main, cex = cex_main, font = 2)
  # late_rcp45_mean
  plot_x(x[[2]])
  mtext(side = 3, "B (RCP 4.5 [L])", adj = adj_main, line = line_main, cex = cex_main, font = 2)
  # mid_rcp85_mean
  plot_x(x[[3]])
  mtext(side = 3, "C (RCP 8.5 [M])", adj = adj_main, line = line_main, cex = cex_main, font = 2)
  # late_rcp85_mean
  plot_x(x[[4]])
  mtext(side = 3, "D (RCP 8.5 [L])", adj = adj_main, line = line_main, cex = cex_main, font = 2)
  
  # Add legend
  add_legend <- FALSE
  if(add_legend){
    pn <- par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    data_legend <- data.frame(x = col_param$breaks[1:(length(col_param$breaks)-1)], col = col_param$col)
    axis_legend <- prettyGraphics::pretty_axis(side = 4, lim = list(c(-1, 1)), axis = list(at = c(-0.5, 0, 0.5)))
    axis_legend[[1]]$axis$pos <- 1
    axis_legend[[1]]$axis$cex.axis <- 1.25
    TeachingDemos::subplot(prettyGraphics::add_colour_bar(data_legend = data_legend,
                                                          pretty_axis_args = axis_legend,
                                                          mtext_args = list(side = 4, 
                                                                            text = expression(paste(Delta, "CRCA")), 
                                                                            line = 3.5, 
                                                                            cex = 1.25)), 
                           x = 180, y = -77, size = c(0.08, 4), vadj = 0, hadj = 0)
    par(pn)
  }
}


##############################
##############################
#### Make plots 

#### SST (mean)
png("./fig/proj_abund_sst_mean.png", 
     height = 5.5, width = 12, units = "in", res = 600)
pp <- par(mfrow = c(2, 2), oma = c(0, 0, 1, 5), mar = c(0, 2.4, 0, 2.4))
plot_projections(list(sst_mid_rcp45_mean, sst_late_rcp45_mean, sst_mid_rcp45_mean, sst_late_rcp85_mean), zlim = c(-1, 1))
mtext(side = 4, 
      expression(paste("Mean change in CRCA, ", E(Delta ~ CRCA["i,j"]))), 
      cex = 1.5, line = 4, outer = TRUE)
dev.off()

#### SST (pr)
png("./fig/proj_abund_sst_pr.png", 
     height = 5.5, width = 12, units = "in", res = 600)
pp <- par(mfrow = c(2, 2), oma = c(0, 0, 1, 5), mar = c(0, 2.4, 0, 2.4))
plot_projections(list(sst_mid_rcp45_pr, sst_late_rcp45_pr, sst_mid_rcp45_pr, sst_late_rcp85_pr), type = 2)
mtext(side = 4, 
      expression(paste("Pr. spp. w/ declines in CRCA, ", Pr(Delta ~ CRCA["i,j"] < 0))),
      cex = 1.5, line = 4, outer = TRUE)
dev.off()

#### SBT (mean)
png("./fig/proj_abund_sbt_mean.png", 
     height = 5.5, width = 12, units = "in", res = 600)
pp <- par(mfrow = c(2, 2), oma = c(0, 0, 1, 5), mar = c(0, 2.4, 0, 2.4))
plot_projections(list(sbt_mid_rcp45_mean, sbt_late_rcp45_mean, sbt_mid_rcp45_mean, sbt_late_rcp85_mean), zlim = c(-1, 1))
mtext(side = 4, 
      expression(paste("Mean change in CRCA, ", E(Delta ~ CRCA["i,j"]))), 
      cex = 1.5, line = 4, outer = TRUE)
dev.off()


#### End of code. 
##############################
##############################