##############################
##############################
#### process_spptraits.R
#### Edward Lavender (el72@st-andrews.ac.uk)

#### This code:
# Defines a list of species, and associated thermal preferences, for modelling. 
# ... Uses an initial list of > 17,000 species
# ... Focuses on species with depth ranges found on fishbase and selects appropriate species
# ... Checks for synonyms
# ... Focuses on coastal species 
# ... Saves a temporary (reduced) list of species for which aquamaps SDMs are processed 
# ... ... (see process_sdm_aquamaps)
# ... Using processed species distributions, 
# ... ... checks data quality for any species with few predicted occurrences 
# ... ... Gets thermal niche parameters 
# ... For the final list of species, gets the full taxonomic breakdown 

#### Steps preceding this code:
# 1) Raw, starting list of species for consideration (./data-raw/aqqsstrans.csv)
# 2) Associated species distribution model predictions 
# ... raw data: ./data-raw/sdm_aquamaps/
# ... processing based on intermediate species list: process_sdm_aquamaps.R
# ... processed species distributions: ./data/sdm_aquamaps 
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


##############################
##############################
#### Basic dataframe processing 

#### Tidy names
names(spptraits)[1] <- "row"
names(spptraits)[7] <- "number.asc.file"
spptraits$SpecCode  <- NULL
spptraits$X.x       <- NULL
spptraits$Genus     <- NULL
spptraits$Species   <- NULL
spptraits$spgen     <- NULL

#### Tidy column types
spptraits$number.asc.file <- as.character(spptraits$number.asc.file)
spptraits$spp <- as.character(spptraits$spp)

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
spptraits$spp_code            <- fishbase$SpecCode
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
#### Focus on coastal species 

#### Map all raw species distributions for manual examination
make_maps <- FALSE
if(make_maps){
  spptraits$index <- 1:nrow(spptraits)
  spptraits_by_spp <- split(spptraits, spptraits$index)
  coastline <- rnaturalearth::ne_coastline(scale = "small")
  cl <- parallel::makeCluster(12L)
  parallel::clusterExport(cl, c("coastline", "spptraits_by_spp"))
  pbapply::pblapply(1:length(spptraits_by_spp), cl = cl, function(i){
    d <- spptraits_by_spp[[i]]
    r <- raster::raster(paste0("./data-raw/sdm_aquamaps/", d$number.asc.file))
    png(paste0("./fig/raw_sdm_aquamaps/", i, "_", d$spp_accepted_name, ".png"), 
        height = 5, width = 5, units = "in", res = 300)
    raster::plot(r, main = d$spp_accepted_name)
    raster::lines(coastline)
    dev.off()
    # print(d$spp_accepted_name)
    # readline(prompt = "Press [enter] to continue...")
  })
  parallel::stopCluster(cl)
}

#### Define a list of non-coastal species to be dropped
not_coastal <- c("Acanthocybium solandri", 
                 "Allothunnus fallai", 
                 "Euleptorhamphus viridis", 
                 "Phtheirichthys lineatus", 
                 "Exocoetus monocirrhus", 
                 "Prognichthys sealei", 
                 "Cheilopogon atrisignis", 
                 "Cheilopogon cyanopterus", 
                 "Cheilopogon spilonotopterus", 
                 "Hirundichthys albimaculatus", 
                 "Trimma emeryi", 
                 "Lythrypnus gilberti",
                 "Halichoeres leucurus", 
                 "Gnatholepis gymnocara")
length(not_coastal)
all(not_coastal %in% spptraits$spp)
spptraits <- spptraits[!(spptraits$spp %in% not_coastal), ]


##############################
##############################
#### Check species SDM quality 

#### Define processed SDMs
# ... For this intermediate-stage reduced list of species, we will process SDMs
# ... because subsequent processing stages require processed SDMs
# Save intermediate spptraits list
# saveRDS(spptraits, "./data-raw/spptraits_for_process_sdm_aquamaps.rds")
# Implement processing 
# ... process_sdm_aquamaps.R
# Now continue with processing. 

#### Identify species with predicted presence in < 10 cells
# ... as potentially very data poor for manual querying.
# ... [This code takes ~ 3 minutes with 11 cores.]
check_occupancy <- TRUE
if(check_occupancy){
  spptraits$index <- 1:nrow(spptraits)
  cl <- parallel::makeCluster(11L)
  parallel::clusterExport(cl, "spptraits")
  spptraits$occupancy <- 
    pbapply::pblapply(split(spptraits, spptraits$index), cl = cl, function(d){
      # d <- spptraits[1, ]
      file <- paste0("./data/sdm_aquamaps/", d$number.asc.file)
      r <- raster::raster(file)
      # raster::plot(r)
      f <- raster::freq(r)
      f <- data.frame(f)
      n_occurrence <- f$count[which(f$value == 1)]
      return(n_occurrence)
    })
  parallel::stopCluster(cl)
  table(spptraits$occupancy < 10)
}

#### Results 
# There are < 100 species found with predicted occupancy in < 10 cells
spptraits[spptraits$occupancy < 10, "spp"]
# We will exclude these species 
spptraits <- spptraits[spptraits$occupancy >= 10, ]


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

#### Get niche quantiles [~4 minutes with 11 cores]
# Define cluster 
cl <- parallel::makeCluster(11L)
parallel::clusterExport(cl = cl, varlist = c("sst_historical", "sbt_historical"))
# Get niche quantiles: 
# ... a vector of six numbers for each species
# ... sst_t10, sst_t50, sst_t90, sbt_t10, sbt_t50, sbt_t90
spptraits$index <- 1:nrow(spptraits)
niche_quantiles <- 
  pbapply::pblapply(split(spptraits, spptraits$index), cl = cl, function(d){
    ## Load data 
    # d <- spptraits[1, ]
    r <- raster::raster(paste0("./data/sdm_aquamaps/", d$number.asc.file))
    # raster::plot(r)
    ## Get temperatures across species range [when r == 1, otherwise r == NA due to pre-processing]
    sst_historical_for_spp <- raster::mask(sst_historical, r)
    sbt_historical_for_spp <- raster::mask(sbt_historical, r)
    # raster::plot(sst_historical_for_spp)
    # raster::plot(sbt_historical_for_spp)
    ## Get thermal quantiles 
    sst_quant <- as.numeric(raster::quantile(sst_historical_for_spp, c(0.1, 0.5, 0.90), na.rm = TRUE))
    sbt_quant <- as.numeric(raster::quantile(sbt_historical_for_spp, c(0.1, 0.5, 0.90), na.rm = TRUE))
    return(c(sst_quant, sbt_quant))
  })
parallel::stopCluster(cl)
niche_quantiles <- do.call(rbind, niche_quantiles)

#### Check that niche quantiles have been correctly defined for all species:
table(is.na(niche_quantiles))

#### Add quantiles to spptraits
spptraits[, c("sst_t10", "sst_t50", "sst_t90", "sbt_t10", "sbt_t50", "sbt_t90")] <- niche_quantiles


##############################
##############################
#### Get taxonomic hierarchy 

#### Use a loop to get taxonomic hierarchy
# ... (to mitigate patchy internet connections)
# ... [This takes ~ 45 minutes]
t1_taxise <- Sys.time()
spptraits$phylum  <- NA
spptraits$class  <- NA
spptraits$order  <- NA
spptraits$family <- NA
try_ncbi <- FALSE
try_itis <- FALSE
for(i in 1:nrow(spptraits)){
  print(i)
  hier <- taxize::classification(spptraits$spp, db = "worms")
  hier <- hier[[1]]$name[hier[[1]]$rank %in% c("Phylum", "Class", "Order", "Family")]
  spptraits[i, c("phylum", "class", "order", "family")] <- hier
  if(try_ncbi){
    spptraits[i, c("phylum", "class", "order", "family")] <- 
      taxize::tax_name(spptraits$spp_accepted_name[i], 
                       get = c("phylum", "class", "order", "family"), db = "ncbi", 
                       message = FALSE)[, c("phylum", "class", "order", "family")]
  }
  if(try_itis){
    if(is.na(spptraits$phylum[i])){
      spptraits[i, c("phylum", "class", "order", "family")] <-
        taxize::tax_name(spptraits$spp_accepted_name[i], 
                         get = c("phylum", "class", "order", "family"), db = "itis", 
                         message = FALSE)[, c("phylum", "class", "order", "family")]
    }
  }
}
beepr::beep(10)
t2_taxise <- Sys.time()
difftime(t2_taxise, t1_taxise)
# beepr::beep(10)

#### Check that taxonomic levels have been successfully queried
spptraits[is.na(spptraits$phylum), ]
table(is.na(spptraits$phylum))
table(is.na(spptraits$class))
table(is.na(spptraits$order))
table(is.na(spptraits$family))

taxize::tax_name(spptraits$spp_accepted_name[i], 
                 get = c("phylum", "class", "order", "family"), db = "itis", 
                 message = FALSE)[, c("phylum", "class", "order", "family")]




# Define .csv comprising species for which manual querying is necessary 
save <- FALSE
if(save){
  write.csv(spptraits[is.na(spptraits$phylum), c("spp", "phylum", "class", "order", "family")], 
            "./data-raw/spptraits_for_manual_taxonomy_to_fill.csv", row.names = FALSE, quote = FALSE)
}
# Load .csv with manually defined taxonomy and add to spptraits
manual <- read.csv("./data-raw/spptraits_for_manual_taxonomy_filled.csv")
index <- match(spptraits$spp, manual$spp)
spptraits$phylum
spptraits$class
spptraits$family




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