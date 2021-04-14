##############################
##############################
#### plot_raster()

#' A helper function for plotting rasters 
plot_raster <- function(x, 
                        mask = NULL, 
                        zlim = NULL, 
                        gen_cols = prettyGraphics::pretty_cols_brewer,
                        profile_x = c(181, 220), profile_y = c(-90, 90),
                        add_legend = TRUE,...){
  
  paa <-  list(side = 1:4, control_axis = list(lwd.ticks = 0, labels = FALSE, lwd = 2))
  axis_args <- list(cex.axis = 1.25)
  sp <- c(0.98, 0.99, 0.2, 0.8)
  if(!is.null(mask)) x <- raster::mask(x = x, mask = mask)
  if(is.null(zlim)) zlim <- raster::cellStats(x, range)
  col_param <- gen_cols(zlim,...)
  if(add_legend){
    prettyGraphics::pretty_map(x, 
                               pretty_axis_args = paa,
                               add_rasters = list(x = x, smallplot = sp, axis.args = axis_args, 
                                                  zlim = zlim, breaks = col_param$breaks, col = col_param$col), 
                               add_polys = list(x = land))
  } else {
    prettyGraphics::pretty_map(x, 
                               pretty_axis_args = paa,
                               add_rasters = list(x = x, plot_method = raster::image, legend = FALSE,
                                                  zlim = zlim, breaks = col_param$breaks, col = col_param$col), 
                               add_polys = list(x = land))
  }
  add_profile <- TRUE
  if(add_profile){
    prettyGraphics::add_profile_lat(x = x, 
                                    ylim = c(-90, 90),
                                    add_fit = list(breaks = col_param$breaks, cols = col_param$col, lwd = 3),
                                    x_at = profile_x, y_at = profile_y,
                                    axes = FALSE)
  }
  return(col_param)
}


#### End of code. 
##############################
##############################