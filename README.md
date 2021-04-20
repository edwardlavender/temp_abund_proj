Modelling the impacts of climate change on the abundance of
shallow-water marine fish at a global scale
================
Edward Lavender<sup>1,2\*</sup>, Clive J. Fox<sup>1</sup> and Michael T.
Burrows<sup>1</sup>

<!-- README.md is generated from README.Rmd. Please edit that file -->

[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)

<sup>1</sup> The Scottish Association for Marine Science, Scottish
Marine Institute, Dunstaffnage, Oban, Scotland, PA37 1QA  
<sup>2</sup> Current Address: Centre for Research into Ecological and
Environmental Modelling, The Observatory, University of St Andrews,
Fife, Scotland, KY16 9LZ

<sup>\*</sup> This repository is maintained by Edward Lavender
(<el72@st-andrews.ac.uk>).

## Introduction

This repository contains methods, written in `R`, for Lavender, Fox and
Burrows (in review). Modelling the impacts of climate change on the
abundance of shallow-water marine fish at a global scale.

In this work, we developed an approach based on species thermal niches
to project changes in the relative abundance of more than 2,000
shallow-water marine fish species under future temperature change
scenarios.

## Structure

1.  `data-raw` contains raw data for the project.
      - `spatial` contains spatial information.
          - `ne_110m_coastline` contains world coastline data from
            [Natural
            Earth](https://www.naturalearthdata.com/downloads/110m-physical-vectors/).
      - `sdm_aquamaps` contains raw [Aquamaps](https://www.aquamaps.org)
        species distribution model (SDM) predictions (maps) for modelled
        species, obtained via `get_sdm_aquamaps.R` (see below);
      - `temperature` contains sea surface temperature (SST) and sea
        bottom temperature (SBT) global scale CMIP5 ensemble average
        projections for (a) historical (1956-2005) (b) mid-century
        (2006-2055) and (c) late-century (2050-2099) time scales under
        Representative Concentration Pathways (RCPs) 4.5 and 8.5. All
        projections were obtained from [NOAA’s Climate Change Web
        Portal](https://psl.noaa.gov/ipcc/ocn/). For SST, for the
        historical scenario, the
        [HADISST](https://www.metoffice.gov.uk/hadobs/hadisst/) dataset
        is also included for comparison to modelled scenarios.  
        <br />
2.  `data` contains processed data.
      - `spptraits.rds` defines modelled species.
      - `spatial` contains processed spatial data.
          - `coastline` contains world coastline data generated by
            `process_coastline.R`.
              - `eez` contains exclusive economic zone boundary data
                from the [Flanders Marine
                Institute](https://www.marineregions.org/eez.php) and
                datasets defining the number of modelled cells and
                species in each EEZ, generated by `process_eez.R` and
                `analyse_abund_across_eezs.R`.
              - `map_species_richness.asc` is a raster of the number of
                modelled species, generated by `analyse_spptraits.R`.
              - `sdm_aquamaps` contains processed
                [Aquamaps](https://www.aquamaps.org) SDMs for modelled
                species, generated by `process_aquamaps.R`.
      - `temperature` contains processed SST and SBT temperature
        projections for the time periods and scenarios described above,
        generated by `process_temp.R`.
      - `sensitivity` contains processed ‘sensitivity’ rasters for
        modelled species, generated by `process_sensitivity.R`.
      - `abundance` contains abundance predictions for SST and SBT for
        all temperature scenarios, generated by `project_abund.R` and
        `process_abund.R`.  
        <br />
3.  `R` contains scripts for data processing, projections and analysis.
      - `get_*` scripts get raw data, where necessary.  
      - `process_*` scripts implement data processing.
      - `project_*` scripts implement relative abundance projections.
      - `analyse_*` scripts analyse the data and projections, including
        figure creation.  
      - `helper_*` scripts contain helper functions used in multiple
        scripts.  
        <br />
4.  `fig` contains figures.

Note that the `data-raw`, `data` and `fig` directories are not provided
in the online version of this repository.

## Workflow

### Data acquisition via `get_*`

1.  `get_sdm_aquamaps.R` gets [Aquamaps](https://www.aquamaps.org) SDMs
    for modelled species, via the
    [`aquamapsdata`](https://github.com/raquamaps/aquamapsdata) `R`
    package.

2.  Other raw data are acquired manually (see the links above).

### Data processing via `process_*`

1.  `process_coastline.R` processes the raw coastline data:
      - Forces an extent of {-180, 180, -90, 90} to match species’
        distributions and temperature projections;  
        <br />
2.  `process_spptraits.R` defines a list of modelled species:
      - Focuses on the subset of species with depth ranges found on
        [FishBase](http://www.fishbase.org/search.php);
      - Checks for synonyms;
      - Saves a temporary (reduced) list of species for which
        [Aquamaps](https://www.aquamaps.org) SDMs are acquired and then
        processed (see `get_sdm_aquamaps.R` and
        `process_sdm_aquamaps.R`);
      - Using processed species’ distributions:
          - Checks SDMs;
          - Gets thermal niche parameters from processed temperature
            projections;
          - Gets the full taxonomic breakdown;  
            <br />
3.  `process_sdm_aquamaps.R` processes
    [Aquamaps](https://www.aquamaps.org) species’ distributions:
      - Forces an extent of {-180, 180, -90, 90} to match temperature
        projections;  
      - Replaces 0 for predicted occurrence with NA;
      - Masks land using processed coastline;  
        <br />
4.  `process_temp.R` processes temperature projections:
      - Extracts SST and SBT temperature projections from raw files;
      - Forces an extent of {-180, 180, -90, 90} to match species’
        distributions;
      - Re-samples temperatures to the same spatial resolution as
        species’ distributions;
      - Masks land using processed coastline to match species’
        distributions;  
        <br />
5.  `process_eez.R` processes EEZ data:
      - Gets EEZ areas (in units of modelled grid cells);
      - Gets EEZ latitudinal mid-points;  
        <br />
6.  `process_sensitivity.R` derives species’ sensitivity indices given
    thermal niche parameters derived from (a) SST or (b) SBT. This
    includes:
      - Mean thermal niche width (species thermal range: STR) over
        space;
      - Mean thermal bias (species thermal index (STI) - baseline
        temperatures) over space;
      - Variability in thermal bias over space;

### Projections via `project_*()`

1.  `project_abund_1.R` projects changes in species’ relative abundance
    under future climate change scenarios.
2.  `project_abund_2.R` synthesises predictions across species.

### Analyses via `analyse_*()`

1.  `analyse_spptraits.R` analyses modelled species, including:
    
      - taxonomic breakdown;
      - depth ranges;
      - distribution;
      - commercial importance;

2.  `analyse_example_thermal_niche.R` analyses an example thermal niche.

3.  `analyse_temp.R` analyses SST and SBT temperature projections, for
    all scenarios.

4.  `analyse_sensitivity.R` analyses species’ sensitivity indices.

5.  `analyse_abund_across_globe.R` analyses relative abundance
    projections across the globe.

6.  `analyse_abund_across_eezs.R` analyses relative abundance
    projections across EEZs.

7.  `analyse_scenarios.R` analyses the relative severity of projections
    under RCP 4.5 and RCP 8.5.

## Citation

Lavender, Fox and Burrows (in review). Modelling the impacts of climate
change on the abundance of shallow-water marine fish at a global scale.
