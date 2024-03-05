# FieldStationConservation

This repository contains the codes for: "Tropical field stations yield high conservation return on investment"

## Code_RedListData.R
This code extract the list of tetrapod species whose range overlaps with the field stations sites and summarizes the list by number of species and proportion of threatened species per group in the regions. 

The code uses the following data:
- dataframe with field station location coordinates
- IUCN range data for the four tetrapod classes
- IUCN RL data for the four tetrapod classes
- IUCN RL taxonomic data for the four tetrapod classes

## Code_ForestCoverLossMatchingApproach.R 
This code was used to evaluate whether field stations prevent forest cover loss we used a statistical matching approach (MatchIt) to identify control points with similar characteristics for each field station. We used the Global Forest Change index v1.9 (Hansen et al., 2013) to quantify differences in forest cover loss. We collected values associated with deforestation as well as environmental properties using Google Earth Engine. Covariates included the percentage of tree cover, temperature and precipitation UN-Adjusted Population Density, the Global Human modification dataset and road density. For both the field station points and the control points, a circular buffer of 5 km was used to extract the site-specific values.

The code uses the following data:
- dataframe (stored as xlsx) with field station location coordinates and founding year
- the GRIP global roads database: grip4_total_dens_m_km2.asc: https://www.globio.info/download-grip-dataset
- google earth engine sources (listed in the R code)
- worldclim 2.1 data downloaded using the geodata package
