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
library(rgee)
ee_Initialize(drive = TRUE)
library(MatchIt)
library(sandwich) 
library(ggplot2)
library(sp)
library(tidyverse)
library(raster)
library(geosphere)
library(sf)
library(openxlsx)
library(geodata)
library(geojsonio)
library(exactextractr)

#-------- vars
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set working directory 
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
pa_within = ee$Image$constant(1)$pixelArea()
pa_outside = ee$Image$constant(1)$pixelArea()
pa = ee$FeatureCollection("WCMC/WDPA/current/polygons")
pa_img_mask = ee$Image$constant(1)$clip(pa)$mask()
pa_within = pa_within$updateMask(pa_img_mask)
pa_outside = pa_outside$updateMask(pa_img_mask$Not())
pa_img = pa_within
pa_img = pa_img$addBands(pa_outside)
pa_img = pa_img$rename(c("within_pa","outside_pa"))
pa_img$bandNames()$getInfo()
# m=Map$addLayer(pa_img$select(0), list(min = 0, max = 1, palette = "green"), name = "pa_img")
# m=m+Map$addLayer(pa_img$select(1), list(min = 0, max = 1, palette = "red"), name = "pa_img")
# m

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
gfcc_img = GFCC_2000
gfcc_img = gfcc_img$addBands(GFCC_2005)
gfcc_img = gfcc_img$addBands(GFCC_2010)
gfcc_img = gfcc_img$addBands(GFCC_2015)
gfcc_img = gfcc_img$rename(paste0("GFCC_",c(2000,2005,2010,2015)))
gfcc_img$bandNames()$getInfo()

#------- open forest cover loss layers data layer connections
HansenGFC= ee$Image("UMD/hansen/global_forest_change_2022_v1_10") # forest loss layers Hansen
HansenGFCLoss = HansenGFC$select("loss")
HansenGFCLosstmp = HansenGFCLoss$pixelArea()
HansenGFCLoss = HansenGFCLosstmp$updateMask(HansenGFCLoss$eq(1))
HansenGFCLossYear = HansenGFC$select("lossyear")
floss_img = HansenGFCLoss
for(i in 0:22) {
  tmp = HansenGFCLossYear$gte(i) 
  tmp2 = tmp$pixelArea()
  tmp = tmp2$updateMask(tmp$eq(1))
  floss_img = floss_img$addBands(tmp)
}
floss_img = floss_img$addBands(HansenGFCLoss$pixelArea())
floss_img$bandNames()$getInfo()
floss_img = floss_img$rename(c("losstotal",paste0("loss_from_",2000:2022),"totalarea"))
floss_img$bandNames()$getInfo()
# m=Map$addLayer(floss_img$select(22), list(min = 0, max = 1, palette = "green"))
# m=m+Map$addLayer(floss_img$select(23), list(min = 0, max = 1, palette = "red"))
# m

#------- open population data layer connections
popdens_img = ee$ImageCollection("CIESIN/GPWv411/GPW_UNWPP-Adjusted_Population_Density") %>%
  ee$ImageCollection$filterDate(paste0(2000,'-01-01'), paste0(2000,'-12-31')) %>%
  ee$ImageCollection$map(function(x) x$select("unwpp-adjusted_population_density")) %>%
  ee$ImageCollection$mean() %>% ee$ImageCollection$toBands()
for(yy in c(2005, 2010, 2015, 2020)) {
  tmp = ee$ImageCollection("CIESIN/GPWv411/GPW_UNWPP-Adjusted_Population_Density") %>%
    ee$ImageCollection$filterDate(paste0(yy,'-01-01'), paste0(yy,'-12-31')) %>%
    ee$ImageCollection$map(function(x) x$select("unwpp-adjusted_population_density")) %>%
    ee$ImageCollection$mean() %>% ee$ImageCollection$toBands()
  popdens_img = popdens_img$addBands(tmp)
}
popdens_img = popdens_img$rename(paste0("popDens_",c(2000,2005,2010,2015,2020)))
popdens_img$bandNames()$getInfo()
# m=Map$addLayer(popdens_img$select(0), list(min = 0, max = 1000, palette = c("lightgreen","red")))
# m=m+Map$addLayer(popdens_img$select(4), list(min = 0, max = 1000, palette = c("lightgreen","red")))
# m

#------- global human modification data connection
globalmod_img = ee$ImageCollection("CSP/HM/GlobalHumanModification")  %>%
  ee$ImageCollection$map(function(x) x$select("gHM")) %>%
  ee$ImageCollection$toBands() 

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

#-------- extract stats from gee layers
df = data.frame() # create df to fill
it = nrow(fs_coords)
makeplot=FALSE
fs_n = 1
for(fs_n in 1:it){
  
  print(paste0("progress: ",fs_n,"/",it)) # show progress
  fs_p = cbind(fs_coords$long[fs_n],fs_coords$lat[fs_n]) |> SpatialPoints() # convert fs point to spatial
  sp::proj4string(fs_p) = CRS(WGS84Proj) # set projection
  sample_area_poly = buffer(fs_p, width = 50*1e3) # create sample area buffer (50km)
  #sample_area_poly_excl = buffer(fs_p, width = 5*1e3) # create sample area exclusion buffer (5km)
  sample_area_poly_excl = buffer(fs_p, width = 10*1e3) # create sample area exclusion buffer (10km)
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
  AnnPrec_crop = crop(AnnPrec,extt) 
  RoadDensity_crop = crop(RoadDensity,extt) 
  
  # create 5km buffer around field station and corrosponding control points
  p = st_as_sf(buffer(fs_p, width = 5000))
  p$pt_type = "Fieldstation"
  pl = list(); for(xx in 1:50) pl[[xx]] = st_as_sf(buffer(sample_p[xx], width = 5000))
  pl = do.call(rbind,pl)
  pl$pt_type = "ControlPoint"
  p = rbind(p,pl)
  p$idx = 1:nrow(p)
  p$ID = fs_coords$id[fs_n]
  p_rgee = sf_as_ee(p) 
  
  # rgee plot 
  # m=Map$addLayer(gfcc_img$select(0), list(min = 0, max = 100, palette = c("yellow","green")))
  # m=m+Map$addLayer(p_rgee)
  # m
  
  # extract rgee data for buffered points
  ### pa_img ###
  pa_extract = ee_extract(x = pa_img, y = p_rgee, fun = ee$Reducer$sum(), scale = 1000, tileScale = 1)
  pa_extract$outside_pa = NULL
  pa_extract$within_pa[pa_extract$within_pa>0] = 1
  pa_extract = pa_extract[order(pa_extract$idx),]
  ### gfcc_img ###
  gfcc_extract = ee_extract(x = gfcc_img, y = p_rgee, fun = ee$Reducer$mean(), scale = 1000, tileScale = 1)
  gfcc_extract = gfcc_extract[order(gfcc_extract$idx),]
  ### floss_img ###
  floss_extract = ee_extract(x = floss_img, y = p_rgee, fun = ee$Reducer$sum(), scale = 30, tileScale = 1)
  floss_extract = floss_extract[order(floss_extract$idx),]
  ### popdens_img ###
  popdens_extract = ee_extract(x = popdens_img, y = p_rgee, fun = ee$Reducer$mean(), scale = 1000, tileScale = 1)
  popdens_extract = popdens_extract[order(popdens_extract$idx),]
  ### globalmod_img ###
  globalmod_extract = ee_extract(x = globalmod_img, y = p_rgee, fun = ee$Reducer$mean(), scale = 1000, tileScale = 1)
  globalmod_extract = globalmod_extract[order(globalmod_extract$idx),]

  # extract local data for buffered points
  RoadDensity_crop = rast(RoadDensity_crop)
  RoadDensity_crop = ifel(is.na(RoadDensity_crop),NA,RoadDensity_crop)
  RoadDensity_crop = ifel(is.infinite(RoadDensity_crop),NA,RoadDensity_crop)
  road_extract = data.frame(pt_type=p$pt_type,ID=fs_coords$id[fs_n])
  road_extract$roaddensity = exact_extract(RoadDensity_crop,p,function(values, coverage_fraction){mean(values, na.rm=TRUE)})
  temp_extract = data.frame(pt_type=p$pt_type,ID=fs_coords$id[fs_n])
  temp_extract$annualTemp = exact_extract(rast(AnnTemp_crop),p,function(values, coverage_fraction){mean(values, na.rm=TRUE)})
  prec_extract = data.frame(pt_type=p$pt_type,ID=fs_coords$id[fs_n])
  prec_extract$annualTempPrec = exact_extract(rast(AnnPrec_crop),p,function(values, coverage_fraction){mean(values, na.rm=TRUE)})

  # merge extracted data
  data_extract = cbind(pa_extract,gfcc_extract,floss_extract,popdens_extract,globalmod_extract,road_extract,temp_extract,prec_extract)
  data_extract = data_extract[,unique(names(data_extract))]
  data_extract$lat = fs_coords[fs_n,"lat"]
  data_extract$long = fs_coords[fs_n,"long"]
  data_extract$year = fs_coords[fs_n,"year"]
  data_extract$Region = fs_coords[fs_n,"Region"]
  
  # get actual locations sampling points
  xy = st_coordinates(st_centroid(p))
  data_extract$long2 = xy[,2]
  data_extract$lat2 = xy[,1]
  
  # extract closest time dep variables
  fs_year = fs_coords[fs_n,"year"]
  var = "GFCC_"
  cols = names(data_extract)[grep(var,names(data_extract))]
  years_var = str_remove(cols,var) |> as.numeric()
  colSelect = paste0(var,years_var[which.min(abs(fs_year-years_var))])
  data_extract[,paste0(var,"closestYear")] = data_extract[,colSelect]
  data_extract = data_extract[,!names(data_extract)%in%cols]
  var = "loss_from_"
  cols = names(data_extract)[grep(var,names(data_extract))]
  years_var = str_remove(cols,var) |> as.numeric()
  colSelect = paste0(var,years_var[which.min(abs(fs_year-years_var))])
  data_extract[,paste0(var,"closestYear")] = data_extract[,colSelect]
  data_extract = data_extract[,!names(data_extract)%in%cols]

  # add data to df
  df = rbind(df,data_extract)
    
}
head(df)


#-------- write temporary output
# write.csv(df,file=paste0("../ExtractedData_",date,".csv"))
# save.image(file=paste0("../ExtractedData_",date,".RData"))

#-------- open temp data
# load(file=paste0("../ExtractedData_","2023_11_21",".RData"))

#-------- set forest cover at founding year
df$ForestCovFoundYr = df$GFCC_closestYear
df$RoadDensity = df$roaddensity
df$RoadDensity[is.nan(df$RoadDensity)] = 0
df$roaddensity = NULL
df$ID = as.numeric(df$ID)
df = df[complete.cases(df),]
df$PrecForestCoverLoss =  (df$losstotal / df$totalarea)*100
df$FieldStation = ifelse(df$pt_type=="Fieldstation",1,0)
df$population_density = df$popDens_2020
df$GlobalHumanModification = df$X2016_gHM
df$withinPA = df$within_pa
nrow(df)
head(df)

#-------- check resulting df
head(df)
length(unique(df$ID))
length(which(df$FieldStation==1))

#-------- matching field station points to control points
m.out1 = matchit(FieldStation ~ ForestCovFoundYr + 
                   RoadDensity + 
                   population_density + 
                   GlobalHumanModification + 
                   annualTemp + 
                   withinPA + 
                   annualTempPrec,
            data = df, 
            method = "nearest", 
            exact='ID', 
            ratio = 5)

#-------- checking balance after matching
summary(m.out1, un = FALSE)
par(mfrow=c(1,1))
plot(summary(m.out1), xlim=c(0,1))
plot(m.out1, type = "density", interactive = FALSE, which.xs = c("ForestCovFoundYr","RoadDensity","population_density","GlobalHumanModification","annualTemp","annualTempPrec"))

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
fs.data1$FieldAbsEffect_weighted = NA

#-------- average of control points compared to treatment point
for(i in 1:nrow(fs.data1)){
  cdatasub = c.data1[fs.data1$ID[i]==c.data1$ID,]
  fs.data1$FieldAbsEffect[i] = fs.data1$PrecForestCoverLoss[i] - mean(cdatasub$PrecForestCoverLoss)
  fs.data1$FieldAbsEffect_weighted[i] = fs.data1$PrecForestCoverLoss[i] - weighted.mean(cdatasub$PrecForestCoverLoss,w=(1/cdatasub$distance))
}

#-------- pair-wise diff between control points and treatment point
fs.data1[,paste0("FieldAbsEffect_c",1:5)] = NA
for(i in 1:nrow(fs.data1)){
  cdatasub = c.data1[fs.data1$ID[i]==c.data1$ID,]
  fs.data1[i,paste0("FieldAbsEffect_c",1:5)] = fs.data1$PrecForestCoverLoss[i] - cdatasub$PrecForestCoverLoss
}
mean(c(mean(fs.data1$FieldAbsEffect_c1),
mean(fs.data1$FieldAbsEffect_c2),
mean(fs.data1$FieldAbsEffect_c3),
mean(fs.data1$FieldAbsEffect_c4),
mean(fs.data1$FieldAbsEffect_c5)))

#-------- prepare for histogram plots, average comparison
min = -15
max = 15
fs.data1$CoverlossPlot = fs.data1$FieldAbsEffect_weighted
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

#-------- prepare for histogram plots, pair-wise comparison
min = -15
max = 15
fs.data1_pair = data.frame(Region=rep(fs.data1$Region,5),
                           ID=rep(fs.data1$ID,5),
                           FieldAbsEffect=c(fs.data1$FieldAbsEffect_c1,
                                            fs.data1$FieldAbsEffect_c2,
                                            fs.data1$FieldAbsEffect_c3,
                                            fs.data1$FieldAbsEffect_c4,
                                            fs.data1$FieldAbsEffect_c5))
fs.data1_pair$CoverlossPlot = fs.data1_pair$FieldAbsEffect
fs.data1_pair$CoverlossPlot[fs.data1_pair$Coverloss<min] = min
fs.data1_pair$CoverlossPlot[fs.data1_pair$Coverloss>max] = max
par(mfrow=c(1,4))
xlims = c(-20,20)
ylims = c(0,200)
hist(c(fs.data1_pair$CoverlossPlot[fs.data1_pair$Region=="Africa"]), 
     breaks=seq(-20,20,5)-2.5, xlim=xlims, xaxs="i",yaxs="i", ylim=c(0,300),
     main = "Africa", xlab="Forest cover change %",ylab="count")
hist(c(fs.data1_pair$CoverlossPlot[fs.data1_pair$Region=="Americas"]), 
     breaks=seq(-20,20,5)-2.5, xlim=xlims, xaxs="i",yaxs="i", ylim=ylims,
     main = "Americas", xlab="Forest cover change %",ylab="count")
hist(c(fs.data1_pair$CoverlossPlot[fs.data1_pair$Region=="Asia"]),
     breaks=seq(-20,20,5)-2.5, xlim=xlims, xaxs="i",yaxs="i", ylim=ylims,
     main = "Asia", xlab="Forest cover change %",ylab="count")
hist(c(fs.data1_pair$CoverlossPlot),breaks=c(-15,-9,-3,3,9,15),
     xlim=xlims,ylim=c(0,500), xaxs="i",yaxs="i",
     main = "All", xlab="Forest cover change %",ylab="count")

#-------- create density plots
all = ggplot(fs.data1, aes(x=CoverlossPlot)) + 
  geom_density(aes(y = ..count..),color="black", fill=rgb(0.1,0.1,0.9,0.5)) + 
  theme_classic() + 
  geom_vline(aes(xintercept=0), color="black",linetype="dotted",linewidth=1) +
  xlab("Relative forest cover loss (%)") + 
  scale_x_continuous(expand = c(0, 0), limits=c(-20,20)) + scale_y_continuous(expand = c(0, 0)) +
  ggtitle("All")
am = ggplot(fs.data1[fs.data1$Region=="Americas",], aes(x=CoverlossPlot)) + 
  geom_density(aes(y = ..count..),color="black", fill=rgb(0.1,0.1,0.9,0.5)) + 
  theme_classic() + 
  geom_vline(aes(xintercept=0), color="black",linetype="dotted",linewidth=1) +
  xlab("Relative forest cover loss (%)") + 
  scale_x_continuous(expand = c(0, 0),limits=c(-20,20)) + scale_y_continuous(expand = c(0, 0))+
  ggtitle("Americas")
af = ggplot(fs.data1[fs.data1$Region=="Africa",], aes(x=CoverlossPlot)) + 
  geom_density(aes(y = ..count..),color="black", fill=rgb(0.1,0.1,0.9,0.5)) + 
  theme_classic() + 
  geom_vline(aes(xintercept=0), color="black",linetype="dotted",linewidth=1) +
  xlab("Relative forest cover loss (%)") + 
  scale_x_continuous(expand = c(0, 0),limits=c(-20,20)) + scale_y_continuous(expand = c(0, 0))+
  ggtitle("Africa")
as = ggplot(fs.data1[fs.data1$Region=="Asia",], aes(x=CoverlossPlot)) + 
  geom_density(aes(y = ..count..),color="black", fill=rgb(0.1,0.1,0.9,0.5)) + 
  theme_classic() + 
  geom_vline(aes(xintercept=0), color="black",linetype="dotted",linewidth=1) +
  xlab("Relative forest cover loss (%)") + 
  scale_x_continuous(expand = c(0, 0),limits=c(-20,20)) + scale_y_continuous(expand = c(0, 0))+
  ggtitle("Asia")
out = ggpubr::ggarrange(all,am,af,as,ncol=2,nrow=2)
file_out = glue::glue('out_df_{format(Sys.time(), "%Y_%m_%d_%H%M%S")}.pdf')
pdf(file_out,width=6,height=4)
print(out)
dev.off()

#-------- create density plots pair-wise
all = ggplot(fs.data1_pair, aes(x=FieldAbsEffect)) + 
  geom_density(aes(y = ..count..),color="black", fill=rgb(0.1,0.1,0.9,0.5)) + 
  theme_classic() + 
  geom_vline(aes(xintercept=0), color="black",linetype="dotted",linewidth=1) +
  xlab("Relative forest cover loss (%)") + 
  scale_x_continuous(expand = c(0, 0), limits=c(-20,20)) + scale_y_continuous(expand = c(0, 0)) +
  ggtitle("All")
am = ggplot(fs.data1_pair[fs.data1_pair$Region=="Americas",], aes(x=FieldAbsEffect)) + 
  geom_density(aes(y = ..count..),color="black", fill=rgb(0.1,0.1,0.9,0.5)) + 
  theme_classic() + 
  geom_vline(aes(xintercept=0), color="black",linetype="dotted",linewidth=1) +
  xlab("Relative forest cover loss (%)") + 
  scale_x_continuous(expand = c(0, 0),limits=c(-20,20)) + scale_y_continuous(expand = c(0, 0))+
  ggtitle("Americas")
af = ggplot(fs.data1_pair[fs.data1_pair$Region=="Africa",], aes(x=FieldAbsEffect)) + 
  geom_density(aes(y = ..count..),color="black", fill=rgb(0.1,0.1,0.9,0.5)) + 
  theme_classic() + 
  geom_vline(aes(xintercept=0), color="black",linetype="dotted",linewidth=1) +
  xlab("Relative forest cover loss (%)") + 
  scale_x_continuous(expand = c(0, 0),limits=c(-20,20)) + scale_y_continuous(expand = c(0, 0))+
  ggtitle("Africa")
as = ggplot(fs.data1_pair[fs.data1_pair$Region=="Asia",], aes(x=FieldAbsEffect)) + 
  geom_density(aes(y = ..count..),color="black", fill=rgb(0.1,0.1,0.9,0.5)) + 
  theme_classic() + 
  geom_vline(aes(xintercept=0), color="black",linetype="dotted",linewidth=1) +
  xlab("Relative forest cover loss (%)") + 
  scale_x_continuous(expand = c(0, 0),limits=c(-20,20)) + scale_y_continuous(expand = c(0, 0))+
  ggtitle("Asia")
ggpubr::ggarrange(all,am,af,as,ncol=2,nrow=2)

#--------- get matched data stats
se = function(x) sqrt(var(x) / length(x))
mode = function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
df_stats = data.frame(subset=c("Asia","Africa","Americas","All"),mean=rep(NA,4))
df_stats$mean[1] = mean(fs.data1$FieldAbsEffect_weighted[fs.data1$Region=="Asia"])
df_stats$mean[2] = mean(fs.data1$FieldAbsEffect_weighted[fs.data1$Region=="Africa"])
df_stats$mean[3] = mean(fs.data1$FieldAbsEffect_weighted[fs.data1$Region=="Americas"])
df_stats$mean[4] = mean(fs.data1$FieldAbsEffect_weighted)
df_stats$median[1] = median(fs.data1$FieldAbsEffect_weighted[fs.data1$Region=="Asia"])
df_stats$median[2] = median(fs.data1$FieldAbsEffect_weighted[fs.data1$Region=="Africa"])
df_stats$median[3] = median(fs.data1$FieldAbsEffect_weighted[fs.data1$Region=="Americas"])
df_stats$median[4] = median(fs.data1$FieldAbsEffect_weighted)
df_stats$mode[1] = mode(fs.data1$FieldAbsEffect_weighted[fs.data1$Region=="Asia"])
df_stats$mode[2] = mode(fs.data1$FieldAbsEffect_weighted[fs.data1$Region=="Africa"])
df_stats$mode[3] = mode(fs.data1$FieldAbsEffect_weighted[fs.data1$Region=="Americas"])
df_stats$mode[4] = mode(fs.data1$FieldAbsEffect_weighted)
df_stats$se[1] = se(fs.data1$FieldAbsEffect_weighted[fs.data1$Region=="Asia"])
df_stats$se[2] = se(fs.data1$FieldAbsEffect_weighted[fs.data1$Region=="Africa"])
df_stats$se[3] = se(fs.data1$FieldAbsEffect_weighted[fs.data1$Region=="Americas"])
df_stats$se[4] = se(fs.data1$FieldAbsEffect_weighted)
df_stats$sd[1] = sd(fs.data1$FieldAbsEffect_weighted[fs.data1$Region=="Asia"])
df_stats$sd[2] = sd(fs.data1$FieldAbsEffect_weighted[fs.data1$Region=="Africa"])
df_stats$sd[3] = sd(fs.data1$FieldAbsEffect_weighted[fs.data1$Region=="Americas"])
df_stats$sd[4] = sd(fs.data1$FieldAbsEffect_weighted)
t1 = t.test(fs.data1$FieldAbsEffect_weighted[fs.data1$Region=="Asia"],mu=0)
t2 = t.test(fs.data1$FieldAbsEffect_weighted[fs.data1$Region=="Africa"],mu=0)
t3 = t.test(fs.data1$FieldAbsEffect_weighted[fs.data1$Region=="Americas"],mu=0)
t4 = t.test(fs.data1$FieldAbsEffect_weighted,mu=0)
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
dd = fs.data1$FieldAbsEffect_weighted[fs.data1$Region=="Asia"]
length(dd)
df_stats$PercDataUnderZero[1] = (length(which(dd<0))/length(dd))*100
dd = fs.data1$FieldAbsEffect_weighted[fs.data1$Region=="Africa"]
length(dd)
df_stats$PercDataUnderZero[2] = (length(which(dd<0))/length(dd))*100
dd = fs.data1$FieldAbsEffect_weighted[fs.data1$Region=="Americas"]
length(dd)
df_stats$PercDataUnderZero[3] = (length(which(dd<0))/length(dd))*100
dd = fs.data1$FieldAbsEffect_weighted
length(dd)
df_stats$PercDataUnderZero[4] = (length(which(dd<0))/length(dd))*100

file_out = glue::glue('out_df_{format(Sys.time(), "%Y_%m_%d_%H%M%S")}.csv')
write.csv(df_stats,file_out)

dd_all = fs.data1
mean(dd_all$FieldAbsEffect_weighted)
dd_asia = fs.data1[fs.data1$Region=="Asia",]
dd_africa = fs.data1[fs.data1$Region=="Africa",]
dd_americas = fs.data1[fs.data1$Region=="Americas",]

