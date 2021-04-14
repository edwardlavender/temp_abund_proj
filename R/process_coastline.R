##############################
##############################
#### process_coastline.R
#### Edward Lavender (el72@st-andrews.ac.uk)

#### This code: 
# 1) Processes the 'coastline' data from NE used to mask 
# ... SDMs and temperature projections, so that all are expressed over the same area. 

#### Steps preceding this code
# 1) Get coastline data from Natural Earth (in ./data-raw/)


##############################
##############################
#### Set up 

#### Wipe workspace
rm(list = ls())

#### Load data
coastline <- rgdal::readOGR("./data-raw/spatial/ne_110m_ocean", "ne_110m_ocean")


##############################
##############################
#### Processing

#### Define extent 
coastline@bbox[1] <- -180
coastline@bbox[2] <- -90
coastline@bbox[3] <- 180
coastline@bbox[4] <- 90

#### Save file 
coastline
saveRDS(coastline, "./data/spatial/coastline/coastline.rds")


#### End of code. 
##############################
##############################