##############################
##############################
#### get_sdm_aquamaps.R
#### Edward Lavender (el72@st-andrews.ac.uk)

#### This code:
# 1) Gets Aquamaps species distributions for an intermediate list of species. 

#### Steps preceding this code: 
# 1) An intermediate species list, from process_spptraits.R 


##############################
##############################
#### Set up 

#### Wipe workspace
rm(list = ls())

#### Essential packages
library(magrittr)

#### Load data 
spptraits <- readRDS("./data-raw/spptraits_for_get_sdm_aquamaps.rds")


##############################
##############################
#### Initialise aquamapsdata package 

#### Method
# The aquamapsdata package provides aquamaps data (https://github.com/raquamaps/aquamapsdata)
# To use the package, GPG needs to be installed. The package can be used to download the aquamaps database
# Then package functions can be used to get species' data 
# (including the underlying data e.g., 'occurcells' and the 0.5 dg rasters). 
# Here, we set up the package and download the database. Then we query the database
# ... for each species to get relevant information. 

setup_aquamapsdata <- FALSE
if(setup_aquamapsdata){
  # Test is GPG is installed, which is required for the aquamapsdata package. 
  system("gpg --version")
  # Download and install GPG: https://gpgtools.org
  # Install package 
  remotes::install_github("raquamaps/aquamapsdata", dependencies = TRUE)
  # One-time download of the database:
  # ... source: https://archive.org/download/aquamapsdb/am.db.bz2 -> 
  # ... sink: /var/folders/lx/dhz6yx2n2b7bg93hwz8t97zr0000gq/T/am.db.bz2 -> /Users/el72/Library/Application Support/aquamaps/am.db
  # ... unpacking /var/folders/lx/dhz6yx2n2b7bg93hwz8t97zr0000gq/T/am.db.bz2 to /Users/el72/Library/Application Support/aquamaps/am.db
  aquamapsdata::download_db(force = TRUE)
  # Set database 
  aquamapsdata::default_db("sqlite")
  # <SQLiteConnection>
  #   Path: /Users/el72/Library/Application Support/aquamaps/am.db
  #   Extensions: TRUE
  # Note: a local copy of this database has been stored within this project
  # ... at './data-raw/sdm_aquamaps_db/', as well, for safe-keeping. 
}


##############################
##############################
#### Get species information and maps

#### Register database 
aquamapsdata::default_db("sqlite")

#### Obtain species data
# This takes about 3 minutes
# ... We need to do this in a standard loop (rather than lapply or with foreach) 
# ... because of an environment issue with how aquamapsdata interprets objects. 
t1_query <- Sys.time()
am_by_spp <- list()
for(i in 1:nrow(spptraits)){
  d <- spptraits[i, , drop = FALSE]
  data <- map <- NULL
  data <- aquamapsdata::am_search_exact(Species = d$species, 
                                        Genus = d$genus)
  if(nrow(data) == 1) map <- aquamapsdata::am_raster(data$SpeciesID)
  am_by_spp[[i]] <- list(data = data, map = map)
}
t2_query <- Sys.time()
difftime(t2_query, t1_query)

#### The number of species with data in aquamapsdata package
# This package includes data for 26,399 /33,518 species (see vignette)
# It only includes species with >=10 ‘good cells’ and excludes records of data-poor species
# Here, we have data for most species in our 'intermediate' list but not all species. 
table(sapply(am_by_spp, function(elm) nrow(elm$data)))

#### Get names of species with data 
# Get names 
names(am_by_spp) <- sapply(am_by_spp, function(elm) paste0(elm$data$Genus, "_", elm$data$Species))
# Remove blank species 
am_by_spp <- am_by_spp[-which(names(am_by_spp) == "_")]

#### Define dataframe with species information 
am_tbl <- 
  lapply(am_by_spp, function(elm) elm$data) %>% 
  dplyr::bind_rows() %>%
  data.frame()
# This includes the species information, such as the 
# ... 'number of good cells used to generate species envelope'
am_tbl$OccurCells[1]

#### Define dataframe for map processing 
spptraits_for_process_sdm_aquamaps <- spptraits[spptraits$spp %in% paste(am_tbl$Genus, am_tbl$Species), ]

#### Save these data 
## Save maps 
save <- FALSE
if(save){
  cl <-  parallel::makeCluster(2L)
  parallel::clusterExport(cl, "am_by_spp")
  pbapply::pblapply(1:length(am_by_spp), cl = cl, function(i){
    spp <- names(am_by_spp)[i]
    elm <- am_by_spp[[i]]
    sink <- paste0("./data-raw/sdm_aquamaps/maps/", spp, ".asc")
    raster::writeRaster(elm$map, sink, overwrite = TRUE)
  }) %>% invisible()
  parallel::stopCluster(cl)
  ## Save species table 
  saveRDS(am_tbl, "./data-raw/sdm_aquamaps/tables/am_data.rds")
  message("Data SAVED.")
  # Save list for map processing
  saveRDS(spptraits_for_process_sdm_aquamaps, "./data-raw/spptraits_for_process_sdm_aquamaps.rds")
  
} else message("Data NOT saved.") 


#### End of code. 
##############################
##############################