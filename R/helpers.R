##############################
##############################
#### calc_sum()

#' Calculate the sum of a list of rasters
#' This function calculates the sum of a list of rasters. This approach is designed for situations in which it is extremely slow to brick rasters and then sum rasters in the usual way. Instead, the function loops over each raster in the list, adding it to the previous one, to generate the summary.  
calc_sum <- function(x, verbose = TRUE){
  # Set up verbose
  cat_to_console <- function(..., show = verbose) if(show) cat(paste(..., "\n"))
  if(verbose){
    t_onset <- Sys.time()
    cat_to_console(paste0("calc_mean() called (@ ", t_onset, ")..."))
  }
  # Define starting raster 
  r <- x[[1]]
  lx <- length(x)
  # Update raster 
  for(i in 2:lx){
    svMisc::progress(i, lx)
    r <- sum(r, x[[i]], na.rm = TRUE)
  }
  # Return outputs 
  if(verbose){
    t_end <- Sys.time()
    duration <- difftime(t_end, t_onset, units = "mins")
    cat_to_console(paste0("... calc_mean() call completed (@ ", t_end, ") after ~", round(duration, digits = 2), " minutes."))
  }
  return(r)
}


##############################
##############################
#### load_projections()

#' A helper function to load files 
load_projections <- function(source = root_ab_delta_sst, type = "sst", name = "mid_rcp45", stat = "mean", mask = NULL){
  cat(paste0("Loading '", source, stat, "/ab_delta_", type, "_", name, "_", stat, ".asc'...\n"))
  r <- raster::raster(paste0(source, stat, "/ab_delta_", type, "_", name, "_", stat, ".asc"))
  if(!is.null(mask)) r <- raster::mask(r, mask = mask)
  return(r)
} 


##############################
##############################
#### plot_raster()

#' A helper function for plotting rasters 
plot_raster <- function(x, 
                        mask = NULL, 
                        zlim = NULL, 
                        gen_cols = prettyGraphics::pretty_cols_brewer,
                        profile_x = c(181, 220), profile_y = c(-90, 90),
                        add_legend = TRUE,
                        sp = NULL,...){
  
  paa <-  list(side = 1:4, control_axis = list(lwd.ticks = 0, labels = FALSE, lwd = 2))
  axis_args <- list(cex.axis = 1.25)
  if(is.null(sp)) sp <- c(0.98, 0.99, 0.2, 0.8)
  if(!is.null(mask)) x <- raster::mask(x = x, mask = mask)
  if(is.null(zlim)) zlim <- raster::cellStats(x, range)
  col_param <- gen_cols(zlim,...)
  if(add_legend){
    prettyGraphics::pretty_map(x, 
                               pretty_axis_args = paa,
                               add_rasters = list(x = x, smallplot = sp, axis.args = axis_args, 
                                                  zlim = zlim, breaks = col_param$breaks, col = col_param$col, 
                                                  pretty_axis_args = list()), 
                               add_polys = list(x = land))
  } else {
    prettyGraphics::pretty_map(x, 
                               pretty_axis_args = paa,
                               add_rasters = list(x = x, plot_method = raster::image, legend = FALSE,
                                                  zlim = zlim, breaks = col_param$breaks, col = col_param$col, 
                                                  pretty_axis_args = list()), 
                               add_polys = list(x = land))
  }
  add_profile <- TRUE
  if(add_profile){
    prettyGraphics::add_profile_lat(x = x, 
                                    calc_lwr = NULL, calc_upr = NULL,
                                    ylim = c(-90, 90),
                                    add_fit = list(breaks = col_param$breaks, cols = col_param$col, lwd = 3),
                                    x_at = profile_x, y_at = profile_y,
                                    axes = FALSE)
  }
  return(invisible(col_param))
}



#### End of code. 
##############################
##############################