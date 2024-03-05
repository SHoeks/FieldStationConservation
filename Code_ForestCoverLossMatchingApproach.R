#### Field stations yield high return-on-investment for conservation ####
#-------- authors 
# Selwyn Hoeks

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
library(terra)
library(tidyterra)
library(marginaleffects)

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

#--------- download forest cover, only for plotting
select_plot_idx = 45
select_plot_year = fs_coords$year[select_plot_idx]
forest_plot_idx = which.min(abs(select_plot_year-c(2000,2005,2010,2015)))
if(file.exists("../ForestCoverRast.tif")){
  r = rast("../ForestCoverRast.tif")
}else{
  r = ee_as_rast(gfcc_img$select(forest_plot_idx),scale=10000)
  terra::writeRaster(r,"../ForestCoverRast.tif")
  r = rast("../ForestCoverRast.tif")
}

#--------- global map plot
colors <- colorRampPalette(colors = c("#D0E7D2", "#B0D9B1", "#79AC78", "#618264"),bias = 1.5)
fs_p = cbind(fs_coords$long,fs_coords$lat) |> SpatialPoints() |> st_as_sf()
st_crs(fs_p) = crs(r)
globalMap = ggplot() + 
  geom_spatraster(
    data = r,
    aes(fill = GFCC_2010)) + # just the name of the layer
  geom_sf(data = fs_p, # lines on top of raster
          fill = NA,
          col = "red",
          size = 0.5) +
  scale_fill_gradientn(
    limits = c(0,100),
    colors = colors(100),
    na.value = "transparent") +
  theme_void()

#--------- prepare zoomed in plot
fs_n = 45
fs_p = cbind(fs_coords$long[fs_n],fs_coords$lat[fs_n]) |> SpatialPoints() # convert fs point to spatial
sp::proj4string(fs_p) = CRS(WGS84Proj) # set projection
sample_area_poly = buffer(fs_p, width = 50*1e3) # create sample area buffer (50km)
sample_area_poly_excl = buffer(fs_p, width = 5*1e3) # create sample area exclusion buffer (5km)
ext = extent(sample_area_poly) # check sample area size
xsample = runif(2000, ext[1], ext[2]) # sample x
ysample = runif(2000, ext[3], ext[4]) # sample y
crop_geometry = function(min_long,max_long,min_lat,max_lat){
  return(ee$Geometry$Rectangle(coords = c(min_long,min_lat,max_long,max_lat),proj="EPSG:4326",geodesic=FALSE))
}
cropgoem = crop_geometry(ext[1]-0.3,ext[2]+0.3,ext[3]-0.3,ext[4]+0.3)
r = ee_as_rast(gfcc_img$select(forest_plot_idx),region=cropgoem,via="drive",container="rgee_backup",scale=500,maxPixels=1e+09,skipEmptyTiles=TRUE)
plot(r)
sample_p = SpatialPoints(cbind(xsample,ysample)) # make sample data spatial
sample_p2 = sample_p[which(sp::over(sample_p, sample_area_poly)==1)] # only keep points in sample_area_poly within 50 km buffer
idx = as.vector(which(sp::over(sample_p2, sample_area_poly_excl)==1))
sample_p3 = sample_p2[!1:length(sample_p2)%in%idx] # remove points close to field station 5 km buffer
sample_p4 = sample_p3[sample(1:length(sample_p3),50)] # sample 50 points from remaining valid points
sp::proj4string(sample_p) = CRS(WGS84Proj) # set projection sampled points

#------------ run match, example only, final match is more elaborate (more co-variates etc)
controls = st_as_sf(sample_p4)
fieldst = st_as_sf(fs_p)
st_data = rbind(fieldst,controls)
example_fc = exact_extract(r,st_buffer(st_data,5000),fun="mean")
dfplot=data.frame(FieldStation=c(1,rep(0,50)),ForestCovFoundYr=example_fc,
                  x=st_coordinates(st_data)[,1],y=st_coordinates(st_data)[,2])
mexample = matchit(FieldStation ~ ForestCovFoundYr,
          # RoadDensity + 
          # population_density + 
          # GlobalHumanModification + 
          # annualTemp + 
          # withinPA + 
          # annualTempPrec,
        data = dfplot, 
        method = "nearest", 
        ratio = 10)
mexample = match.data(mexample)
selectP = cbind(mexample$x,mexample$y) |> SpatialPoints() |> st_as_sf()
st_crs(selectP) = crs(r)

#------------ make zoomed plot
LocalMap1 = ggplot() + 
  geom_spatraster(
    data = r,
    aes(fill = GFCC_2010)) + # just the name of the layer
  geom_sf(data = st_as_sf(fs_p), # lines on top of raster
          fill = NA,
          col = "blue",
          size = 2) +
  geom_sf(data = st_as_sf(sample_area_poly), # lines on top of raster
          fill = NA,
          col = "blue",
          size = 10) +
  geom_sf(data = st_as_sf(sample_area_poly_excl), # lines on top of raster
          fill = NA,
          col = "orange",
          size = 10) +
  geom_sf(data = st_as_sf(sample_p), # lines on top of raster
          fill = NA,
          col = "grey",
          size = 0.5) +
  geom_sf(data = st_as_sf(sample_p2), # lines on top of raster
          fill = NA,
          col = "black",
          size = 0.5) +
  geom_sf(data = st_as_sf(sample_p3), # lines on top of raster
          fill = NA,
          col = "purple",
          size = 0.5) +
  geom_sf(data = st_as_sf(sample_p4), # lines on top of raster
          fill = NA,
          col = "red",
          size = 2) +
  geom_sf(data = st_as_sf(selectP), # lines on top of raster
          fill = NA,
          col = "green",
          size = 2) +
  scale_fill_gradientn(
    limits = c(0,100),
    colors = colors(100),
    na.value = "transparent") +
  theme_void()


#-------- extract stats from gee layers
df = data.frame() # create df to fill
it = nrow(fs_coords)
makeplot = FALSE
for(fs_n in 1:it){
  
  print(paste0("progress: ",fs_n,"/",it)) # show progress
  fs_p = cbind(fs_coords$long[fs_n],fs_coords$lat[fs_n]) |> SpatialPoints() # convert fs point to spatial
  sp::proj4string(fs_p) = CRS(WGS84Proj) # set projection
  sample_area_poly = buffer(fs_p, width = 50*1e3) # create sample area buffer (50km)
  sample_area_poly_excl = buffer(fs_p, width = 5*1e3) # create sample area exclusion buffer (5km)
  #sample_area_poly_excl = buffer(fs_p, width = 10*1e3) # create sample area exclusion buffer (10km)
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

#------------ useful when continuing from here
if(!any(grepl("MatchIt",(.packages())))){
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
  library(terra)
  library(tidyterra)
  library(marginaleffects)
  library(flextable)
}

#-------- open temp data
if(length(grep("data_extract",ls()))<1){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set working directory
  load(file=paste0("../ExtractedData_","2023_11_21",".RData"))
}

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
            ratio = 10)

#-------- checking balance after matching
summary(m.out1, un = FALSE)
par(mfrow=c(1,1))
plot(summary(m.out1), xlim=c(0,1))
plot(m.out1, type = "density", interactive = FALSE, which.xs = c("ForestCovFoundYr","RoadDensity","population_density","GlobalHumanModification","annualTemp","annualTempPrec"))

#-------- extract matched data from matching object
m.data1 = match.data(m.out1)
length(unique(m.data1$ID))
nrow(m.data1)

#-------- estimating the treatment effect
m.data1$treat = as.factor(m.data1$FieldStation)
m.data1$Region[m.data1$Region=="Madagascar"] = "Africa"

fit = lm(PrecForestCoverLoss ~ treat * 
        ( ForestCovFoundYr + RoadDensity + population_density + 
          GlobalHumanModification + annualTemp + withinPA + annualTempPrec), 
          data = m.data1, weights = weights)
avg_comparisons(fit, variables = "treat", vcov = ~subclass + ID,
                newdata = m.data1, wts = "weights")

#-------- check mean % change in field station points vs control points
dr = aggregate(PrecForestCoverLoss~Region+FieldStation,m.data1,mean)
dr = dr[order(dr$Region),]
dg = aggregate(PrecForestCoverLoss~FieldStation,m.data1,mean)
dg$Region = "All"
da = rbind(dg,dr)
print(da)

#-------- calculate relative differences, stats for main text
relRegion = data.frame(region=unique(dr$Region),relPrec=rep(NA,3))
relRegion$relPrec = (1-(dr$PrecForestCoverLoss[dr$FieldStation==1]/dr$PrecForestCoverLoss[dr$FieldStation==0]))*100
relAll = data.frame(region='all',relPrec=NA)
relAll$relPrec = (1-(dg$PrecForestCoverLoss[dg$FieldStation==1]/dg$PrecForestCoverLoss[dg$FieldStation==0]))*100
rel = rbind(relAll,relRegion)
relft = flextable(rel)
print(relft, preview = 'docx')
asciify(rel)

#-------- check absolute difference between control points and field station
c.data1 = m.data1[m.data1$FieldStation==0,]
fs.data1 = m.data1[m.data1$FieldStation==1,]
fs.data1$FieldAbsEffect = NA
fs.data1$FieldAbsEffect_weighted = NA
fs.data1$PrecForestCoverLossweightedmeanControl = NA

#-------- average of control points compared to treatment point
for(i in 1:nrow(fs.data1)){
  cdatasub = c.data1[fs.data1$ID[i]==c.data1$ID,]
  fs.data1$FieldAbsEffect[i] = fs.data1$PrecForestCoverLoss[i] - mean(cdatasub$PrecForestCoverLoss)
  fs.data1$PrecForestCoverLossweightedmeanControl[i] = weighted.mean(cdatasub$PrecForestCoverLoss,w=(1/cdatasub$distance))
  fs.data1$FieldAbsEffect_weighted[i] = fs.data1$PrecForestCoverLoss[i] - weighted.mean(cdatasub$PrecForestCoverLoss,w=(1/cdatasub$distance))
}
fs.data1 = fs.data1[!duplicated(fs.data1$ID),]

#-------- perform t-test, compare field station to weighted means control points
# tdf = data.frame(fsval=fs.data1$PrecForestCoverLoss,
#                  FieldAbsEffect_weighted=fs.data1$FieldAbsEffect_weighted,
#                  ctrlval=fs.data1$PrecForestCoverLoss-fs.data1$FieldAbsEffect_weighted,
#                  region=fs.data1$Region)
# varequal = TRUE
# t1=t.test(x=tdf$fsval,y=tdf$ctrlval,paired = TRUE, var.equal = varequal) # not significant, no differences 
# t2=t.test(x=tdf$fsval[tdf$region=="Africa"],y=tdf$ctrlval[tdf$region=="Africa"],paired = TRUE, var.equal = varequal) # not significant, no differences 
# t3=t.test(x=tdf$fsval[tdf$region=="Americas"],y=tdf$ctrlval[tdf$region=="Americas"],paired = TRUE, var.equal = varequal) # not significant, no differences 
# t4=t.test(x=tdf$fsval[tdf$region=="Asia"],y=tdf$ctrlval[tdf$region=="Asia"],paired = TRUE, var.equal = varequal) # not significant, no differences 
# tdf2 = data.frame(pval=c(t1$p.value,t2$p.value,t3$p.value,t4$p.value),region=c("All","Africa","Americas","Asia"))
# print(tdf2)

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

#-------- create density plots
all = ggplot(fs.data1, aes(x=CoverlossPlot)) + 
  geom_density(aes(y = ..count..),color="black", fill=rgb(0.1,0.1,0.9,0.5)) + 
  theme_classic() + 
  geom_vline(aes(xintercept=0), color="black",linetype="dotdash",linewidth=0.75) +
  xlab("Relative Forest Cover Loss (%)") + 
  ylab("# of Field Stations") + 
  scale_x_continuous(expand = c(0, 0), limits=c(-20,20)) + 
  scale_y_continuous(expand = c(0, 0),limits=c(0,20))+
  ggtitle("Global") + 
  theme(plot.title = element_text(hjust = 0.5))
am = ggplot(fs.data1[fs.data1$Region=="Americas",], aes(x=CoverlossPlot)) + 
  geom_density(aes(y = ..count..),color="black", fill=rgb(0.1,0.1,0.9,0.5)) + 
  theme_classic() + 
  geom_vline(aes(xintercept=0), color="black",linetype="dotdash",linewidth=0.75) +
  xlab("Relative Forest Cover Loss (%)") + 
  ylab("") + 
  scale_x_continuous(expand = c(0, 0),limits=c(-20,20)) + 
  scale_y_continuous(expand = c(0, 0),limits=c(0,20))+
  ggtitle("Americas") + 
  theme(plot.title = element_text(hjust = 0.5))
af = ggplot(fs.data1[fs.data1$Region=="Africa",], aes(x=CoverlossPlot)) + 
  geom_density(aes(y = ..count..),color="black", fill=rgb(0.1,0.1,0.9,0.5)) + 
  theme_classic() + 
  geom_vline(aes(xintercept=0), color="black",linetype="dotdash",linewidth=0.75) +
  xlab("Relative Forest Cover Loss (%)") + 
  ylab("") + 
  scale_x_continuous(expand = c(0, 0),limits=c(-20,20)) + 
  scale_y_continuous(expand = c(0, 0),limits=c(0,20))+
  ggtitle("Africa") +
  theme(plot.title = element_text(hjust = 0.5))
as = ggplot(fs.data1[fs.data1$Region=="Asia",], aes(x=CoverlossPlot)) + 
  geom_density(aes(y = ..count..),color="black", fill=rgb(0.1,0.1,0.9,0.5)) + 
  theme_classic() + 
  geom_vline(aes(xintercept=0), color="black",linetype="dotdash",linewidth=0.75) +
  xlab("Relative Forest Cover Loss (%)") + 
  ylab("") + 
  scale_x_continuous(expand = c(0, 0),limits=c(-20,20)) + 
  scale_y_continuous(expand = c(0, 0),limits=c(0,20))+
  ggtitle("Asia") +
  theme(plot.title = element_text(hjust = 0.5))
out = ggpubr::ggarrange(all,am,af,as,ncol=4,nrow=1)
print(out)

#--------- get number of matched points below 0
dfout = data.frame(subset=c("Asia","Africa","Americas","All"),PercDataUnderZero=rep(NA,4),NPoints=rep(NA,4))
dd = fs.data1$FieldAbsEffect_weighted[fs.data1$Region=="Asia"]
dfout$PercDataUnderZero[1] = (length(which(dd<0))/length(dd))*100
dfout$NPoints[1] = length(dd)
dd = fs.data1$FieldAbsEffect_weighted[fs.data1$Region=="Africa"]
dfout$PercDataUnderZero[2] = (length(which(dd<0))/length(dd))*100
dfout$NPoints[2] = length(dd)
dd = fs.data1$FieldAbsEffect_weighted[fs.data1$Region=="Americas"]
dfout$PercDataUnderZero[3] = (length(which(dd<0))/length(dd))*100
dfout$NPoints[3] = length(dd)
dd = fs.data1$FieldAbsEffect_weighted
dfout$PercDataUnderZero[4] = (length(which(dd<0))/length(dd))*100
dfout$NPoints[4] = length(dd)
print(dfout)

#--------- data for map
mapdata = fs.data1[,c("ID","lat","long","CoverlossPlot","Region")]
mapdataNA = df[df$FieldStation==1,]
mapdataNA = mapdataNA[,c("ID","lat","long","Region")]
mapdataNA$CoverlossPlot = NA
mapdataNA$CoverlossPlot[match(mapdata$ID,mapdataNA$ID)] = mapdata$CoverlossPlot
mapdataNA$n = 1:nrow(mapdataNA)
mapdataNA_sf = st_as_sf(mapdataNA, coords = c("long", "lat"), crs = WGS84Proj)
View(mapdataNA)

# --------- create map
ggmap = ggplot() +
  geom_sf(data = world.sf,lwd=0.2,bg="#D7D8D7") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        plot.margin = unit(c(0,0,0,0), unit = "cm")) + 
  geom_sf(data = mapdataNA_sf, aes(color = CoverlossPlot), fill = "black", size=3) + 
  scale_color_distiller(palette = "RdYlBu", direction = -1) +
  xlim(c(-120,120)) + ylim(c(-55,40)) 
ggmap




