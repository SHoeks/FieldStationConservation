#### Field stations yield high return-on-investment for conservation ####
#-------- authors 
# XXXXX

#-------- description 
# this code was used to evaluate whether field stations prevent forest cover loss
# we used a statistical matching approach (MatchIt) to identify control points with 
# similar characteristics for each field station. 
# we used the Global Forest Change index v1.9 (Hansen et al., 2013) to quantify differences in forest cover loss 
# we collected values associated with deforestation as well as environmental properties using Google Earth Engine. 
# covariates included the percentage of tree cover, temperature and precipitation
# UN-Adjusted Population Density, the Global Human modification dataset and road density. 
# For both the field station points and the control points, a circular buffer of 5 km was used to 
# extract the site-specific values

# the code uses the following data:
# - dataframe (stored as xlsx) with field station location coordinates and founding year
# - the GRIP global roads database: grip4_total_dens_m_km2.asc: https://www.globio.info/download-grip-dataset
# - google earth engine sources
# - worldclim 2.1 data downloaded using the geodata package

#-------- packages
library(MatchIt)
library(sandwich) 
library(ggplot2)
library(sp)
library(tidyverse)
library(rgee)
library(raster)
library(geosphere)
library(sf)
library(openxlsx)
library(geodata)

#-------- init
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set working directory 
rm(list = ls())
ee_Initialize(drive = TRUE)
WGS84Proj<-"+proj=longlat +datum=WGS84 +no_defs"
date = format(Sys.time(), "%Y_%m_%d")

#------- open xlsx with field station data + convert to google earth engine (gee)
fs_data=read.xlsx("../Field Stations-Spatial Data 1.2.xlsx", sheet = 1) # open field station file
fs_coords = data.frame(id=fs_data$Survey.Number,lat=fs_data$Latitude,long=fs_data$Longitude,year=fs_data$Founded.in,Region=fs_data$Region)
fs_coords$year = ifelse(is.na(fs_coords$year),2000,fs_coords$year)
coord_vec = rep(NA,nrow(fs_coords)*2)
coord_vec[seq(1,10000,2)[1:length(fs_coords$long)]] = fs_coords$long
coord_vec[seq(2,10000,2)[1:length(fs_coords$lat)]] = fs_coords$lat
points = ee$Geometry$MultiPoint(coords = coord_vec) 

#------ get protected area (PA) polygons from google earth engine (gee)
PAs = ee$FeatureCollection("WCMC/WDPA/current/polygons") # open data layer connection to PAs
filteredPAs = PAs$filterBounds(points) # filter PAs, only keeping those that contain a field station xy
filteredPAs_sf = ee_as_sf(filteredPAs)
filteredPAs$size()$getInfo() # check filteredPAs 
PA_plot = sf::as_Spatial(sf::st_collection_extract(sf::st_geometry(filteredPAs_sf), "POLYGON")) # find combination points and polygons
fs_coords$poly_id = NA
for(i in 1:length(PA_plot)) {
  idx <- sp::point.in.polygon(point.x = fs_coords$long, point.y = fs_coords$lat, pol.x = raster::geom(PA_plot[i])[,5], pol.y = raster::geom(PA_plot[i])[,6])
  fs_coords$poly_id[which(idx==1)] = i
}

#------- open forest cover layers data layer connections
GFCC_2000 = ee$ImageCollection("NASA/MEASURES/GFCC/TC/v3") %>%
  ee$ImageCollection$filterDate(paste0(2000,'-01-01'), paste0(2000,'-12-31')) %>%
  ee$ImageCollection$map(function(x) x$select("tree_canopy_cover")) %>%
  ee$ImageCollection$mean() %>% ee$ImageCollection$toBands() 
GFCC_2005 = ee$ImageCollection("NASA/MEASURES/GFCC/TC/v3") %>%
  ee$ImageCollection$filterDate(paste0(2005,'-01-01'), paste0(2005,'-12-31')) %>%
  ee$ImageCollection$map(function(x) x$select("tree_canopy_cover")) %>%
  ee$ImageCollection$mean() %>% ee$ImageCollection$toBands() 
GFCC_2010 = ee$ImageCollection("NASA/MEASURES/GFCC/TC/v3") %>%
  ee$ImageCollection$filterDate(paste0(2010,'-01-01'), paste0(2010,'-12-31')) %>%
  ee$ImageCollection$map(function(x) x$select("tree_canopy_cover")) %>%
  ee$ImageCollection$mean() %>% ee$ImageCollection$toBands() 
GFCC_2015 = ee$ImageCollection("NASA/MEASURES/GFCC/TC/v3") %>%
  ee$ImageCollection$filterDate(paste0(2015,'-01-01'), paste0(2015,'-12-31')) %>%
  ee$ImageCollection$map(function(x) x$select("tree_canopy_cover")) %>%
  ee$ImageCollection$mean() %>% ee$ImageCollection$toBands() 

#------- open forest cover loss layers data layer connections
HansenGFC= ee$Image("UMD/hansen/global_forest_change_2021_v1_9") # forest loss layers Hansen
HansenGFCLoss = HansenGFC$select("loss")

#------- open population data layer connections
PopDens = ee$ImageCollection("CIESIN/GPWv411/GPW_UNWPP-Adjusted_Population_Density") %>%
  ee$ImageCollection$filterDate(paste0(2020,'-01-01'), paste0(2020,'2015-12-31')) %>%
  ee$ImageCollection$map(function(x) x$select("unwpp-adjusted_population_density")) %>%
  ee$ImageCollection$mean() %>% ee$ImageCollection$toBands() # 
PopDens$bandNames()$getInfo()

#------- global human modification data connection
GlobalHumanModification= ee$ImageCollection("CSP/HM/GlobalHumanModification")  %>%
  ee$ImageCollection$map(function(x) x$select("gHM")) %>%
  ee$ImageCollection$toBands() 

#------- temperature data connection
ERA5_TEMP1980 = ee$ImageCollection("ECMWF/ERA5/MONTHLY") %>%
  ee$ImageCollection$filterDate(paste0(1980,'-01-01'), paste0(1980,'-12-31')) %>%
  ee$ImageCollection$map(function(x) x$select("mean_2m_air_temperature")) %>%
  ee$ImageCollection$mean()
ERA5_TEMP2000 = ee$ImageCollection("ECMWF/ERA5/MONTHLY") %>%
  ee$ImageCollection$filterDate(paste0(2000,'-01-01'), paste0(2000,'-12-31')) %>%
  ee$ImageCollection$map(function(x) x$select("mean_2m_air_temperature")) %>%
  ee$ImageCollection$mean()
ERA5_TEMP2020 = ee$ImageCollection("ECMWF/ERA5/MONTHLY") %>%
  ee$ImageCollection$filterDate(paste0(2020,'-01-01'), paste0(2020,'-12-31')) %>%
  ee$ImageCollection$map(function(x) x$select("mean_2m_air_temperature")) %>%
  ee$ImageCollection$mean()

#------- load road density, external data, not from gee
RoadDensity = raster("../ExternalData/GRIP4_density_total/grip4_total_dens_m_km2.asc")
RoadDensity = log10(RoadDensity)

#------- download and load tavg climate overall, external data, not from gee
res = 2.5
tavgDir = '../ExternalData/wc2.1_10m/'
tavgFiles = str_pad(1:12,2,'left',pad="0") %>% 
  paste0(tavgDir,'wc2.1_10m_tavg_',.,'.tif')
if(!all(file.exists(paste0(tavgFiles)))){
  AnnTemp = geodata::worldclim_global('tavg', res=res, path='../ExternalData/')
}else{
  AnnTemp = terra::rast(tavgFiles)
}
AnnTemp = mean(AnnTemp)
AnnTemp = raster(AnnTemp)

#------- download and load prec climate overall, external data, not from gee
res = 2.5
precDir = '../ExternalData/wc2.1_10m/'
precFiles = str_pad(1:12,2,'left',pad="0") %>% 
  paste0(tavgDir,'wc2.1_10m_prec_',.,'.tif')
if(!all(file.exists(paste0(precFiles)))){
  AnnPrec = geodata::worldclim_global('prec', res=res, path='../ExternalData/')
}else{
  AnnPrec = terra::rast(precFiles)
}
AnnPrec = mean(AnnPrec)
AnnPrec = raster(AnnPrec)

#------- temperature seasonality (not used in analysis)
SeasTemp = AnnTemp
SeasTemp[] = 0

#-------  merge bands into one image with multiple bands
GFCC_comp = GFCC_2000$addBands(GFCC_2005)
GFCC_comp = GFCC_comp$addBands(GFCC_2010)
GFCC_comp = GFCC_comp$addBands(GFCC_2015)
GFCC_comp = GFCC_comp$addBands(PopDens)
GFCC_comp = GFCC_comp$addBands(GlobalHumanModification)
GFCC_comp = GFCC_comp$addBands(ERA5_TEMP1980)
GFCC_comp = GFCC_comp$addBands(ERA5_TEMP2000)
GFCC_comp = GFCC_comp$addBands(ERA5_TEMP2020)
GFCC_comp$bandNames()$getInfo()

#-------- extact stats from gee layers
df = data.frame() # create df to fill
it = nrow(fs_coords)
makeplot=FALSE
for(fs_n in 1:it){
  
  print(paste0("progress: ",fs_n,"/",it)) # show progress
  n_poly=fs_coords$poly_id[fs_n] # get matching PA 
  fs_p = cbind(fs_coords$long[fs_n],fs_coords$lat[fs_n]) |> SpatialPoints() # convert fs point to spatial
  sp::proj4string(fs_p) = CRS(WGS84Proj) # set projection
  sample_area_poly = buffer(fs_p, width = 50*1e3) # create sample area buffer (50km)
  sample_area_poly_excl = buffer(fs_p, width = 5*1e3) # create sample area exclusion buffer (5km)
  ext = extent(sample_area_poly) # check sample area size
  xsample = runif(2000, ext[1], ext[2]) # sample x
  ysample = runif(2000, ext[3], ext[4]) # sample y
  sample_p = SpatialPoints(cbind(xsample,ysample)) # make sample data spatial
  sp::proj4string(sample_p) = CRS(WGS84Proj) # set projection sampled points
  if(makeplot){
    plot(sample_area_poly)
    plot(sample_area_poly_excl,add=TRUE,col="red")
    plot(sample_p,add=TRUE,col="grey")
  }
  sample_p = sample_p[which(sp::over(sample_p, sample_area_poly)==1)] # only keep points in sample_area_poly within 50 km buffer
  if(makeplot) plot(sample_p,add=TRUE,col="black")
  idx = as.vector(which(sp::over(sample_p, sample_area_poly_excl)==1))
  sample_p = sample_p[!1:length(sample_p)%in%idx] # remove points close to field station 5 km buffer
  if(makeplot) plot(sample_p,add=TRUE,col="blue")
  sample_p = sample_p[sample(1:length(sample_p),50)] # sample 50 points from remaining valid points
  
  # crop external rasters to increase speed
  extt = extent(fs_coords$long[fs_n]-3,fs_coords$long[fs_n]+3,fs_coords$lat[fs_n]-3,fs_coords$lat[fs_n]+3)
  AnnTemp_crop = crop(AnnTemp,extt)
  SeasTemp_crop = crop(SeasTemp,extt)
  annPrec_crop = crop(annPrec,extt) 
  RoadDensity_crop = crop(RoadDensity,extt) 
  
  # create 5km buffer around field station point
  fs_buffer = buffer(fs_p, width = 5000) |> st_as_sf() |> sf_as_ee() 
  
  # extract data for field station buffered area 
  GFCC_stats_fs = GFCC_comp$reduceRegion(
    reducer = ee$Reducer$mean(),
    geometry = fs_buffer$geometry(),
    scale = 30,
    maxPixels = 1e9)$getInfo() |>
    unlist()
  
  # get hansen fc loss
  HansenFCoverloss = HansenGFCLoss$reduceRegion(
    reducer = ee$Reducer$sum(),
    geometry = fs_buffer$geometry(),
    scale = 30,
    maxPixels = 1e9)$getInfo() |>
    unlist()
  
  # extract data from external cropped rasters
  tmp_buf = buffer(fs_p, width = 5000) # 5km buffer 
  RoadDensity_extracted = RoadDensity_crop |> mask(tmp_buf, inverse=FALSE) |> values()
  RoadDensity_extracted = mean(RoadDensity_extracted[is.finite(RoadDensity_extracted)],na.rm=TRUE)
  if(is.na(RoadDensity_extracted)){
    tmp_buf_l = buffer(fs_p, width = 10000)
    RoadDensity_extracted = RoadDensity_crop |> mask(tmp_buf_l, inverse=FALSE) |> values()
    RoadDensity_extracted = mean(RoadDensity_extracted[is.finite(RoadDensity_extracted)],na.rm=TRUE)
  }
  AnnTemp_cell  = AnnTemp_crop  |> mask(tmp_buf,inverse=FALSE) |> cellStats(mean)
  SeasTemp_cell = SeasTemp_crop |> mask(tmp_buf,inverse=FALSE) |> cellStats(mean)
  annPrec_cell  = annPrec_crop  |> mask(tmp_buf,inverse=FALSE) |> cellStats(mean)
  
  # check if point is in PA
  withinPA = 0
  if(!is.na(n_poly)) {
    if(!is.na(sp::over(fs_p, PA_plot[n_poly]))) withinPA = 1
  }
  
  # create data
  d_fs = data.frame(ID = fs_coords$id[fs_n],
                    Region = fs_coords$Region[fs_n],
                    FieldStation = 1, 
                    withinPA = withinPA,
                    FoundingYear = fs_coords$year[fs_n],
                    x = fs_coords$long[fs_n],
                    y = fs_coords$lat[fs_n],
                    HansenFCoverloss = HansenFCoverloss,
                    GFCC_2000 = GFCC_stats_fs[1], 
                    GFCC_2005 = GFCC_stats_fs[2], 
                    GFCC_2010 = GFCC_stats_fs[3], 
                    GFCC_2015 = GFCC_stats_fs[4],
                    population_density = GFCC_stats_fs[5],
                    GlobalHumanModification = GFCC_stats_fs[6],
                    Temp1980 = GFCC_stats_fs[7],
                    Temp2000 = GFCC_stats_fs[8],
                    Temp2020 = GFCC_stats_fs[9],
                    AnnTemp  = AnnTemp_cell,
                    SeasTemp = SeasTemp_cell,
                    annPrec  = annPrec_cell,
                    RoadDensity = RoadDensity_extracted,
                    row.names = 1)
  
  # add data to df
  df = rbind(df,d_fs)
  
  # loop over 25 sampled points
  for(i in 1:length(sample_p)){
    
    # add same size buffer to sampled points
    sample_buffer = buffer(sample_p[i], width = 5000) |> st_as_sf() |> sf_as_ee() 
    
    # extract data for sampled points buffered areas
    GFCC_stats_sampled = GFCC_comp$reduceRegion(
      reducer = ee$Reducer$mean(),
      geometry = sample_buffer$geometry(),
      scale = 30,
      maxPixels = 1e9)$getInfo() |>
      unlist()
    
    # get hansen fc loss
    HansenFCoverloss = HansenGFCLoss$reduceRegion(
      reducer = ee$Reducer$sum(),
      geometry = sample_buffer$geometry(),
      scale = 30,
      maxPixels = 1e9)$getInfo() |>
      unlist()
    
    # extract data from external cropped rasters
    tmp_buf = buffer(sample_p[i], width = 5000)
    RoadDensity_extracted = RoadDensity_crop |> mask(tmp_buf, inverse=FALSE) |> values()
    RoadDensity_extracted = mean(RoadDensity_extracted[is.finite(RoadDensity_extracted)],na.rm=TRUE)
    if(is.na(RoadDensity_extracted)){
      tmp_buf_l = buffer(fs_p, width = 10000)
      RoadDensity_extracted = RoadDensity_crop |> mask(tmp_buf_l, inverse=FALSE) |> values()
      RoadDensity_extracted = mean(RoadDensity_extracted[is.finite(RoadDensity_extracted)],na.rm=TRUE)
    }
    AnnTemp_cell  = AnnTemp_crop  |> mask(tmp_buf,inverse=FALSE) |> cellStats(mean)
    SeasTemp_cell = SeasTemp_crop |> mask(tmp_buf,inverse=FALSE) |> cellStats(mean)
    annPrec_cell  = annPrec_crop  |> mask(tmp_buf,inverse=FALSE) |> cellStats(mean)
    
    # check if point is in PA
    withinPA = 0
    if(!is.na(n_poly)) {
      if(!is.na(sp::over(sample_p[i], PA_plot[n_poly]))) withinPA = 1
    }
    
    # create data
    d_sample = data.frame(ID = fs_coords$id[fs_n],
                          Region = fs_coords$Region[fs_n],
                          FieldStation = 0, 
                          withinPA = withinPA,
                          FoundingYear = fs_coords$year[fs_n],
                          x = sample_p[i]@coords[1],
                          y = sample_p[i]@coords[2],
                          HansenFCoverloss = HansenFCoverloss, 
                          GFCC_2000 = GFCC_stats_sampled[1], 
                          GFCC_2005 = GFCC_stats_sampled[2], 
                          GFCC_2010 = GFCC_stats_sampled[3], 
                          GFCC_2015 = GFCC_stats_sampled[4], 
                          population_density = GFCC_stats_sampled[5],
                          GlobalHumanModification = GFCC_stats_sampled[6],
                          Temp1980 = GFCC_stats_sampled[7],
                          Temp2000 = GFCC_stats_sampled[8],
                          Temp2020 = GFCC_stats_sampled[9],
                          AnnTemp  = AnnTemp_cell,
                          SeasTemp = SeasTemp_cell,
                          annPrec  = annPrec_cell,
                          RoadDensity = RoadDensity_extracted,
                          row.names = 1)
    
    # add data to df
    df = rbind(df,d_sample)
    
  }
  
}
head(df)

#-------- write temporary output
write.csv(df,file=paste0("../ExtractedData_",date,".csv"))
save.image(file=paste0("../ExtractedData_",date,".RData"))

#-------- open temp data
#load(file=paste0("ExtractedData_",date,".RData"))

#-------- set forest cover at founding year
df$ForestCovFoundYr = NA
colOffset = grep("GFCC_2000",names(df))-1
for(i in 1:nrow(df)){
  colIdx = c(df$FoundingYear[i]-2000,df$FoundingYear[i]-2005,
    df$FoundingYear[i]-2010,df$FoundingYear[i]-2015) |>
    abs() |> which.min()
  df$ForestCovFoundYr[i] = df[i,colOffset+colIdx]
}
df$RoadDensity[is.nan(df$RoadDensity)] = 0
df$ID = as.numeric(df$ID)
df = df[complete.cases(df),]
df$PrecForestCoverLoss =  ( (df$HansenFCoverloss / 90174) * 100 )
nrow(df)

#-------- check resulting df
head(df)
length(unique(df$ID))
length(which(df$FieldStation==1))

#-------- matching field station points to control points
m.out1 = matchit(FieldStation ~ ForestCovFoundYr + 
                   RoadDensity + 
                   population_density + 
                   GlobalHumanModification + 
                   AnnTemp + 
                   annPrec,
            data = df, 
            method = "nearest", 
            exact='ID', 
            ratio = 5)

#-------- checking balance after matching
summary(m.out1, un = FALSE)
par(mfrow=c(1,1))
plot(summary(m.out1), xlim=c(0,1))
plot(m.out1, type = "density", interactive = FALSE, which.xs = c("ForestCovFoundYr","RoadDensity","population_density","GlobalHumanModification","AnnTemp","annPrec"))

#-------- extract matched data from matching object
m.data1 = match.data(m.out1)
length(unique(m.data1$ID))
nrow(m.data1)

#-------- check mean % change in field station points vs control points
m.data1$Region[m.data1$Region=="Madagascar"] = "Africa"
dr = aggregate(PrecForestCoverLoss~Region+FieldStation,m.data1,mean)
dr = dr[order(dr$Region),]
dg = aggregate(PrecForestCoverLoss~FieldStation,m.data1,mean)

#-------- check absolute difference between control points and field station
c.data1 = m.data1[m.data1$FieldStation==0,]
fs.data1 = m.data1[m.data1$FieldStation==1,]
fs.data1$FieldAbsEffect = NA

#-------- average of control points compared to treatment point
for(i in 1:nrow(fs.data1)){
  cdatasub = c.data1[fs.data1$ID[i]==c.data1$ID,]
  fs.data1$FieldAbsEffect[i] = fs.data1$PrecForestCoverLoss[i] - mean(cdatasub$PrecForestCoverLoss)
}

#-------- prepare for histogram plots
min = -15
max = 15
fs.data1$CoverlossPlot = fs.data1$FieldAbsEffect
fs.data1$CoverlossPlot[fs.data1$Coverloss<min] = min
fs.data1$CoverlossPlot[fs.data1$Coverloss>max] = max
par(mfrow=c(1,4))
xlims = c(-20,20)
ylims = c(0,50)
hist(c(fs.data1$CoverlossPlot[fs.data1$Region=="Africa"]), 
     breaks=seq(-20,20,5)-2.5, xlim=xlims, xaxs="i",yaxs="i", ylim=ylims,
     main = "Africa", xlab="Forest cover change %",ylab="count")
hist(c(fs.data1$CoverlossPlot[fs.data1$Region=="Americas"]), 
     breaks=seq(-20,20,5)-2.5, xlim=xlims, xaxs="i",yaxs="i", ylim=ylims,
     main = "Americas", xlab="Forest cover change %",ylab="count")
hist(c(fs.data1$CoverlossPlot[fs.data1$Region=="Asia"]),
     breaks=seq(-20,20,5)-2.5, xlim=xlims, xaxs="i",yaxs="i", ylim=ylims,
     main = "Asia", xlab="Forest cover change %",ylab="count")
hist(c(fs.data1$CoverlossPlot),breaks=c(-15,-9,-3,3,9,15),
     xlim=xlims,ylim=c(0,100), xaxs="i",yaxs="i",
     main = "All", xlab="Forest cover change %",ylab="count")

#-------- create density plots
all = ggplot(fs.data1, aes(x=FieldAbsEffect)) + 
  geom_density(aes(y = ..count..),color="black", fill=rgb(0.1,0.1,0.9,0.5)) + 
  theme_classic() + 
  geom_vline(aes(xintercept=0), color="black",linetype="dotted",linewidth=1) +
  xlab("Relative forest cover loss (%)") + 
  scale_x_continuous(expand = c(0, 0), limits=c(-20,20)) + scale_y_continuous(expand = c(0, 0)) +
  ggtitle("All")
am = ggplot(fs.data1[fs.data1$Region=="Americas",], aes(x=FieldAbsEffect)) + 
  geom_density(aes(y = ..count..),color="black", fill=rgb(0.1,0.1,0.9,0.5)) + 
  theme_classic() + 
  geom_vline(aes(xintercept=0), color="black",linetype="dotted",linewidth=1) +
  xlab("Relative forest cover loss (%)") + 
  scale_x_continuous(expand = c(0, 0),limits=c(-20,20)) + scale_y_continuous(expand = c(0, 0))+
  ggtitle("Americas")
af = ggplot(fs.data1[fs.data1$Region=="Africa",], aes(x=FieldAbsEffect)) + 
  geom_density(aes(y = ..count..),color="black", fill=rgb(0.1,0.1,0.9,0.5)) + 
  theme_classic() + 
  geom_vline(aes(xintercept=0), color="black",linetype="dotted",linewidth=1) +
  xlab("Relative forest cover loss (%)") + 
  scale_x_continuous(expand = c(0, 0),limits=c(-20,20)) + scale_y_continuous(expand = c(0, 0))+
  ggtitle("Africa")
as = ggplot(fs.data1[fs.data1$Region=="Asia",], aes(x=FieldAbsEffect)) + 
  geom_density(aes(y = ..count..),color="black", fill=rgb(0.1,0.1,0.9,0.5)) + 
  theme_classic() + 
  geom_vline(aes(xintercept=0), color="black",linetype="dotted",linewidth=1) +
  xlab("Relative forest cover loss (%)") + 
  scale_x_continuous(expand = c(0, 0),limits=c(-20,20)) + scale_y_continuous(expand = c(0, 0))+
  ggtitle("Asia")
ggpubr::ggarrange(all,am,af,as,ncol=2,nrow=2)

#--------- get matched data stats
se = function(x) sqrt(var(x) / length(x))
df_stats = data.frame(subset=c("Asia","Africa","Americas","All"),mean=rep(NA,4))
df_stats$mean[1] = mean(fs.data1$FieldAbsEffect[fs.data1$Region=="Asia"])
df_stats$mean[2] = mean(fs.data1$FieldAbsEffect[fs.data1$Region=="Africa"])
df_stats$mean[3] = mean(fs.data1$FieldAbsEffect[fs.data1$Region=="Americas"])
df_stats$mean[4] = mean(fs.data1$FieldAbsEffect)
df_stats$median[1] = median(fs.data1$FieldAbsEffect[fs.data1$Region=="Asia"])
df_stats$median[2] = median(fs.data1$FieldAbsEffect[fs.data1$Region=="Africa"])
df_stats$median[3] = median(fs.data1$FieldAbsEffect[fs.data1$Region=="Americas"])
df_stats$median[4] = median(fs.data1$FieldAbsEffect)
df_stats$se[1] = se(fs.data1$FieldAbsEffect[fs.data1$Region=="Asia"])
df_stats$se[2] = se(fs.data1$FieldAbsEffect[fs.data1$Region=="Africa"])
df_stats$se[3] = se(fs.data1$FieldAbsEffect[fs.data1$Region=="Americas"])
df_stats$se[4] = se(fs.data1$FieldAbsEffect)
df_stats$sd[1] = sd(fs.data1$FieldAbsEffect[fs.data1$Region=="Asia"])
df_stats$sd[2] = sd(fs.data1$FieldAbsEffect[fs.data1$Region=="Africa"])
df_stats$sd[3] = sd(fs.data1$FieldAbsEffect[fs.data1$Region=="Americas"])
df_stats$sd[4] = sd(fs.data1$FieldAbsEffect)
t1 = t.test(fs.data1$FieldAbsEffect[fs.data1$Region=="Asia"],mu=0)
t2 = t.test(fs.data1$FieldAbsEffect[fs.data1$Region=="Africa"],mu=0)
t3 = t.test(fs.data1$FieldAbsEffect[fs.data1$Region=="Americas"],mu=0)
t4 = t.test(fs.data1$FieldAbsEffect,mu=0)
df_stats$pvalue[1] = t1$p.value
df_stats$pvalue[2] = t2$p.value
df_stats$pvalue[3] = t3$p.value
df_stats$pvalue[4] = t4$p.value
df_stats$ci95_lower[1] = t1$conf.int[1]
df_stats$ci95_lower[2] = t2$conf.int[1]
df_stats$ci95_lower[3] = t3$conf.int[1]
df_stats$ci95_lower[4] = t4$conf.int[1]
df_stats$ci95_upper[1] = t1$conf.int[2]
df_stats$ci95_upper[2] = t2$conf.int[2]
df_stats$ci95_upper[3] = t3$conf.int[2]
df_stats$ci95_upper[4] = t4$conf.int[2]
dd = fs.data1$FieldAbsEffect[fs.data1$Region=="Asia"]
df_stats$PercDataUnderZero[1] = (length(which(dd<0))/length(dd))*100
dd = fs.data1$FieldAbsEffect[fs.data1$Region=="Africa"]
df_stats$PercDataUnderZero[2] = (length(which(dd<0))/length(dd))*100
dd = fs.data1$FieldAbsEffect[fs.data1$Region=="Americas"]
df_stats$PercDataUnderZero[3] = (length(which(dd<0))/length(dd))*100
dd = fs.data1$FieldAbsEffect
df_stats$PercDataUnderZero[4] = (length(which(dd<0))/length(dd))*100
