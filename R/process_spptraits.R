##############################
##############################
#### process_spptraits.R
#### Edward Lavender (el72@st-andrews.ac.uk)

#### This code:
# Defines a list of species, and associated thermal preferences, for modelling. 
# ... Uses an initial list of > 17,000 species
# ... Focuses on species with depth ranges found on fishbase and selects appropriate species (based on depth)
# ... Checks for synonyms
# ... Saves a temporary (reduced) list of species for which aquamaps SDMs are obtained
# ... Using aquamaps distributions, 
# ... ... Focus on species for which there are no obvious artefacts in SDMs
# ... ... Gets thermal niche parameters 
# ... For the final list of species, gets the full taxonomic breakdown 

#### Steps preceding this code:
# 1) Raw, starting list of species for consideration (./data-raw/aqqsstrans.csv)
# 2) Processed, historical climatology predictions (./data/temperature/)


##############################
##############################
#### Set up 

#### Wipe workspace
rm(list = ls())

#### Essential packages
library(magrittr)

#### Load data
# We start with a list of > 17,000 species for consideration 
spptraits <- read.csv("./data-raw/aqqsstrans.csv")
coastline <- readRDS("./data/spatial/coastline/coastline.rds")


##############################
##############################
#### Basic dataframe processing 

#### Define species names
# spp is the full binomial name, genus and species refer to the parts of the name
spptraits <- spptraits[, c("spp", "Genus", "Species")]
colnames(spptraits) <- c("spp", "genus", "species")
# Join species names via "_"
spptraits$spp_key <- apply(spptraits[, c("genus", "species")], 1, function(x) paste0(x, collapse = "_"))
# Define associated file names 
spptraits$spp_key_asc <- paste0(spptraits$spp_key, ".asc")
head(spptraits)

#### Basic information 
# 17,343 species 
nrow(spptraits) 
nrow(spptraits) == length(unique(spptraits$spp))


##############################
##############################
#### Get fishbase data 

#### Define blank dataframe to save output of query to fishbase
fishbase <- rfishbase::species(spptraits$spp)
str(fishbase)

#### Add key fields to spptraits
# spptraits$spp_code            <- fishbase$SpecCode
spptraits$depth_range_shallow <- fishbase$DepthRangeShallow
spptraits$depth_range_deep    <- fishbase$DepthRangeDeep
spptraits$importance          <- fishbase$Importance

#### Focus on species recognised by fishbase
# Number of species without fishbase data (7,365 species)
table(is.na(fishbase$Species))
# Focus on remaining species
# ... This is necessary for filtering based on depth (see below). 
spptraits <- spptraits[which(!is.na(fishbase$Species)), ]


##############################
##############################
#### Focus on species with appropriate depth ranges

#### Fishbase depth data 
# For an explanation of the different depth values given, 
# ... see http://www.fishbase.org/manual/english/fishbasethe_species_table.htm
# We will restrict the species list to those species exclusively found
# ... above 50 m. 

#### Subset based on depth 
depth_threshold <- 50 
spptraits <- spptraits[which(spptraits$depth_range_deep < depth_threshold), ]
nrow(spptraits)


##############################
##############################
#### Check for and drop synonyms 

#### Examine synonyms via fishbase 
# Get full synonyms list
syns_full <- rfishbase::synonyms()

#### Force 'accepted' names 
# Define a list of accepted names 
spptraits$index <- 1:nrow(spptraits)
spptraits$spp_accepted_name <- sapply(split(spptraits, spptraits$index), function(d){
  # For each species, identify all the synonyms 
  # d <- split(spptraits, spptraits$spp)[[1523]]
  syns_for_spp <- syns_full %>% dplyr::filter(synonym == d$spp)
  # If there are synonyms, then check the accepted name and force this if necessary
  if(nrow(syns_for_spp) > 1){
    accepted_name <- syns_for_spp$synonym[syns_for_spp$Status == "accepted name"]
    if(d$spp != accepted_name) d$spp <- accepted_name
  }
  return(d$spp)
})

#### Check and remove any duplicates
table(spptraits$spp != spptraits$spp_accepted_name)
table(duplicated(spptraits$spp_accepted_name))
if(any(duplicated(spptraits$spp_accepted_name))){
  spptraits <- spptraits[!duplicated(spptraits$spp_accepted_name), ]
}

#### Results
# Accepted names used in all cases.
# No synonyms identified. 


##############################
##############################
#### Obtain aquamaps SDMs 

#### Define SDMs
# ... For this intermediate-stage reduced list of species, we will obtain SDMs
# ... because subsequent processing stages require processed SDMs
# Save intermediate spptraits list
# saveRDS(spptraits, "./data-raw/spptraits_for_get_sdm_aquamaps.rds")
# Now get aquamaps data, via 
# ... get_sdm_aquamaps.R

#### Process 'raw' species distributions for this intermediate list of species
# ... via process_sdm_aquamaps.R
# ... For ease of visualisation, below. 

#### Focus on species for which aquamaps data are available
# The aquamapsdata package only provides data for ~20,000 species for which 
# ... the number of 'good' cells is >= 10. 
am_data <- readRDS("./data-raw/sdm_aquamaps/tables/am_data.rds")
am_data$spp <- paste(am_data$Genus, am_data$Species)
table(spptraits$spp %in% am_data$spp)
spptraits <- spptraits[spptraits$spp %in% am_data$spp, ]

#### Qualitatively examine species distributions
# ... see maps created by process_sdm_aquamaps. 

## Define a vector of species that have possible artefacts in their distributions 
spp_to_ck <- c('Prognichthys sealei',
               'Phtheirichthys lineatus', 
               'Parexocoetus mento',
               'Monodactylus argenteus',
               'Hirundichthys albimaculatus',
               'Exocoetus monocirrhus',
               'Euleptorhamphus viridis',
               'Cypselurus naresii',
               'Cololabis adocetus',
               'Cheilopogon spilonotopterus',
               'Cheilopogon cyanopterus',
               'Cheilopogon atrisignis',
               'Allothunnus fallai',
               'Acanthocybium solandri')

## Manually examine these distributions again 
check <- FALSE
if(check){
  lapply(spp_to_ck, function(spp){
    spp <- "Gnatholepis gymnocara"
    spptraits_for_spp <- spptraits[spptraits$spp == spp, ]
    raster::plot(raster::raster(paste0("./data/sdm_aquamaps/", spptraits_for_spp$spp_key_asc)),
                 main = spp)
    raster::lines(coastline)
    readline(prompt = "Press [enter] to continue...")
  }) %>% invisible()
}

## Remove these species 
nrow(spptraits)
spptraits <- spptraits[!(spptraits$spp %in% spp_to_ck), ]
nrow(spptraits)

#### Add 'occurcells' (the number of 'good' aquamaps cells to spptraits)
# (Data from the aquamapsdata package are only provided for species with more than 10 good cells.)
spptraits$occur_cells <- am_data$OccurCells[match(spptraits$spp, am_data$spp)]
utils.add::basic_stats(spptraits$occur_cells)


##############################
##############################
#### Get thermal niche parameters 

#### Method 
# 1) Load historical climatology (for SST or SBT)
# 2) For each species, load species distribution map 
# 3) In historical climatology, focus on areas of predicted presence (everywhere else has been masked to NA)
# 3) Get quantiles of variation in temperature across these areas

#### Get baseline temperatures 
sst_historical <- raster::raster("./data/temperature/sst/historical/historical.asc")
sbt_historical <- raster::raster("./data/temperature/sbt/historical/historical.asc")

#### Examine baseline temperatures 
pp <- par(mfrow = c(1, 2))
raster::plot(sst_historical)
raster::plot(sbt_historical)
par(pp)

#### Get niche quantiles [~16 minutes with 11 cores]
# We will extract thermal affinities based on (a) weighted quantiles across the full predicted distribution
# and (b) un-weighted quantiles across the distribution where the probability of presence is >= 0.5.

## Define cluster 
cl <- parallel::makeCluster(11L)
parallel::clusterExport(cl = cl, varlist = c("sst_historical", "sbt_historical"))

## Get niche quantiles: 
spptraits$index <- 1:nrow(spptraits)
niche_quantiles <- 
  pbapply::pblapply(split(spptraits, spptraits$index), cl = cl, function(d){
    
    ## Load data 
    # d <- spptraits[1, ]
    r <- raster::raster(paste0("./data/sdm_aquamaps/", d$spp_key_asc))
    r_p50 <- r
    r_p50[r_p50 < 0.5] <- NA
    # raster::plot(r)
    
    ## Get temperatures across species range [r is always > 0 following pre-processing]
    sst_historical_for_spp <- raster::mask(sst_historical, r)
    sbt_historical_for_spp <- raster::mask(sbt_historical, r)
    sst_historical_for_spp_p50 <- raster::mask(sst_historical, r_p50)
    sbt_historical_for_spp_p50 <- raster::mask(sbt_historical, r_p50)
    # raster::plot(sst_historical_for_spp)
    # raster::plot(sbt_historical_for_spp)
    # raster::plot(sst_historical_for_spp_p50)
    # raster::plot(sbt_historical_for_spp_p50)
    
    ## Weighted thermal quantiles
    probs <- c(0.1, 0.5, 0.90)
    wts <- r[]
    sst_quant <- as.numeric(Hmisc::wtd.quantile(sst_historical_for_spp[], wts, probs = probs, normwt = FALSE, na.rm = TRUE))
    sbt_quant <- as.numeric(Hmisc::wtd.quantile(sbt_historical_for_spp[], wts, probs = probs, normwt = FALSE, na.rm = TRUE))
    
    ## Un-weighted thermal quantiles based on threshold
    sst_quant_p50 <- as.numeric(raster::quantile(sst_historical_for_spp_p50, probs, na.rm = TRUE))
    sbt_quant_p50 <- as.numeric(raster::quantile(sbt_historical_for_spp_p50, probs, na.rm = TRUE))
    
    ## Return outputs 
    return(c(sst_quant, sst_quant_p50, sbt_quant, sbt_quant_p50))
  })
parallel::stopCluster(cl)
niche_quantiles <- do.call(rbind, niche_quantiles)

#### Check that niche quantiles have been correctly defined for all species:
table(is.na(niche_quantiles))

#### Add quantiles to spptraits
spptraits[, c("sst_t10", "sst_t50", "sst_t90", 
              "sst_t10_p50", "sst_t50_p50", "sst_t90_p50",
              "sbt_t10", "sbt_t50", "sbt_t90", 
              "sbt_t10_p50", "sbt_t50_p50", "sbt_t90_p50"
              )] <- niche_quantiles


##############################
##############################
#### Get taxonomic hierarchy 

fb_taxa <- rfishbase::load_taxa()
index <- match(spptraits$spp, fb_taxa$Species)
spptraits$class  <- fb_taxa$Class[index]
spptraits$order  <- fb_taxa$Order[index]
spptraits$family <- fb_taxa$Family[index]

#### Check that taxonomic levels have been successfully queried
table(is.na(spptraits$class))
table(is.na(spptraits$order))
table(is.na(spptraits$family))


##############################
##############################
#### Final adjustments

spptraits$index <- 1:nrow(spptraits)
spptraits$species_epiphet <- spptraits$species
spptraits <- spptraits[, c("index", "spp", 
                           "occur_cells",
                           "sst_t10", "sst_t50", "sst_t90", 
                           "sst_t10_p50", "sst_t50_p50", "sst_t90_p50", 
                           "sbt_t10", "sbt_t50", "sbt_t90", 
                           "sbt_t10_p50", "sbt_t50_p50", "sbt_t90_p50",
                           "depth_range_shallow", "depth_range_deep", "importance",
                           "species_epiphet", "genus", "family", "order", "class", 
                           "spp_key", "spp_key_asc")
                       ]



##############################
##############################
#### Save file

#### Save list of species for modelling 
save <- FALSE
if(save){
  message("'ssptraits' HAS been saved.")
  saveRDS(spptraits, "./data/spptraits.rds")
} else {
  message("NOTE 'ssptraits' has NOT been saved because save = FALSE.")
}


#### End of code. 
##############################
##############################