#### Field stations yield high return-on-investment for conservation ####
#Author: Luca Santini

#This code extract the list of tetrapod species whose range overlaps with the field stations sites
#and summarizes the list by number of species and proportion of threatened species per group in the regions
#the code uses the following data:
# - dataframe with field station location coordinates
# - IUCN range data for the four tetrapod classes
# - IUCN RL data for the four tetrapod classes
# - IUCN RL taxonomic data for the four tetrapod classes

library(raster)
library(maptools)
library(sf)
library(sp)
library(ggplot2)

TAXONOMIC_CLASSES=c('AMPHIBIANS', 'REPTILES', 'BIRDS', 'MAMMALS')
DIR='' #path to the files

#load field stations data
data<-read.csv('Field Stations-Spatial Data 1.2.csv')
data<-data[complete.cases(data$Longitude),]

for (tc in TAXONOMIC_CLASSES) {

TAXONOMIC_CLASS=TAXONOMIC_CLASSES[tc]

#load range data for the class
MR<-st_read(paste0(DIR,TAXONOMIC_CLASS,'.shp'))#, proj4string=CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')

#load taxonomic and Red List data for the class
RL<-read.csv(paste0(DIR, TAXONOMIC_CLASS,'_June2022.csv'))
tax<-read.csv(paste(DIR, TAXONOMIC_CLASS,'_taxonomy.csv'))

ids<-unique(MR$id_no) #extrat ids

#generate grid
res=0.5
mat<-matrix(0, ncol=360/res, nrow=180/res)
grid<-raster(mat, xmn=-180, xmx=180, ymn=-90, ymx=90, crs=CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'))

xy<-data.frame(x=data$Longitude, y=data$Latitude)
xy<-xy[!is.na(xy$x),]
coordinates(xy)<-~x+y
xy<-st_as_sf(xy, coords = c("x","y"), )
st_crs(xy)<-st_crs(MR)
xy$id<-1:nrow(xy)

out<-matrix(NA, ncol=length(ids), nrow=nrow(xy))
colnames(out)<-paste0('id_', ids)

#GENERATE MATRIX WITH FIELD STATION LOCATIONS (COORDINATES) IN ROWS AND SPECIES AS COLUMNS (1 if present, NA otherwise)

for (i in 1:length(ids)) {
	id=ids[i]
	print(paste(i, 'out of', length(ids)))
		range<-MR[MR$id_no==id,]

		int<-try(st_intersection(xy, range), silent=TRUE)

		if(class(int)[1]=='try-error') { #se da errore vai di raster
		
		if(is.na(extent(range)[1])){print('skip geometry issue'); next}
		if(length(range)==0) {next}
		grid2<-grid
		grid2<-crop(grid2, extent(range))
		range<-try(rasterize(range, grid2, getCover=TRUE), silent=TRUE)
		range[range>0]<-1
		range[range==0]<-NA
		if(class(range)=='try-error'){print('skip geometry issue'); next}
		range[!is.na(range)]<-1
		out[,paste0('id_', id)]<-extract(range, xy)
		next
		}
		out[int$id,paste0('id_', id)]<-1

		}

out2<-cbind(st_coordinates(xy),out)

#LIST SPECIES PRESENT PER FIELD STATION LOCATION

out3<-data.frame()

for (i in 1:nrow(out2)){
	print(i)
	xy<-out2[i, c('X','Y')]
	sps<-out2[i, 3:ncol(out2)]
	sps<-names(sps)[!is.na(sps)]
	if(length(sps)==0){next}
	sps<-as.numeric(gsub('id_','',sps))
	tmp<-cbind(xy[1], xy[2], sps)
	out3<-rbind(out3,tmp)
}
colnames(out3)[1:2]<-c('x','y')

#ADD TAXONOMIC INFORMATION AND RL INFORMATIONS

RL<-merge(RL[,c('internalTaxonId', 'scientificName', 'redlistCategory')], tax[,c('internalTaxonId', 'orderName', 'familyName', 'genusName')], by='internalTaxonId', all.x=TRUE)   

out3<-merge(out3, RL, by.x='sps', by.y='internalTaxonId', all.x=TRUE) 

write.csv(out3, DIR, TAXONOMIC_CLASS, 'FieldStations2_Vect.csv', row.names=FALSE)

}


#SUMMARY FIGURE

rept<-read.csv(paste0(DIR,'ReptilesFieldStations2_Vect.csv'))
rept$Group<-'Reptiles'
rept<-rept[!rept$familyName %in% c('CHELONIIDAE', 'DERMOCHELYIDAE'),]
rept<-rept[!rept$genusName %in% c('Aipysurus', 'Emydocephalus', 'Ephalophis', 'Hydrelaps', 'Hydrophis', 'Laticauda', 'Parahydrophis'),]
amph<-read.csv(paste0(DIR,'AmphibiansFieldStations2_Vect.csv'))
amph$Group<-'Amphibians'
mamm<-read.csv(paste0(DIR,'MammalsFieldStations2_Vect.csv'))
mamm$Group<-'Mammals'
mamm$Group<-ifelse(mamm$orderName=='PRIMATES', 'Primates', mamm$Group)
birds<-read.csv(paste0(DIR,'BirdsFieldStations2_Vect.csv'))
birds$Group<-'Birds'

data<-rbind(rept, amph, mamm, birds)#mamm

data$RL<-ifelse(data$redlistCategory=='Least Concern', 'LC', 
	ifelse(data$redlistCategory=='Near Threatened', 'NT',
		ifelse(data$redlistCategory=='Vulnerable', 'VU',
			ifelse(data$redlistCategory=='Endangered', 'EN',
				ifelse(data$redlistCategory=='Critically Endangered', 'CR', 'DD')))))
data$RL<-as.factor(data$RL)
data$RL <- factor(data$RL, levels = c('LC','NT','VU','EN','CR', 'DD'))


data$Region<-ifelse(data$x< -30, 'Neotropics',
	ifelse(data$x> -30 & data$x< 55, 'Africa', 
			ifelse(data$x> 55, 'Asia', NA)))#)

#make table for ggplot
data$ID<-paste0(data$scientificName, '_', data$Region) 
data<-data[!duplicated(data$ID),] 

data$N<-1
dt<-aggregate(N~RL+Region+Group, FUN='sum', data=data)

dt$Group <- factor(dt$Group, levels = c('Amphibians','Reptiles','Birds','Mammals', 'Primates'))

#PERCENTAGES OF THREATENED
dt$regG<-paste0(dt$Region,'_',dt$Group)

regG<-unique(dt$regG)

dt2<-data.frame()
for (i in 1:length(regG)){
	regGtmp<-regG[i]
	regGds<-dt[dt$regG==regGtmp,]
	nT<-regGds[regGds$RL %in% c('LC','NT'),]
	T<-regGds[!regGds$RL %in% c('LC','NT'),]
	T$Perc<-100*T$N/sum(nT$N, T$N)
	T$N<-sum(T$N) 
	T$N[2:nrow(T)]<-NA
	T$PercTot<-sum(T$Perc)
	T$PercTot[2:nrow(T)]<-NA
	dt2<-rbind(dt2, T)
}

dt2$RL <- factor(dt2$RL, levels = c('DD','VU','EN','CR'))
dt2$Group <- factor(dt2$Group, levels = c('Amphibians','Reptiles','Birds','Mammals', 'Primates'))

pdf(paste0(DIR, 'PercentageRLperContinent3.pdf'), width=12, height=4)
ggplot(dt2, aes(y=Perc, x=Group, fill=RL)) + 
    geom_bar(position="stack", stat="identity") + facet_wrap( ~ Region, nrow = 1) +
    ylab('% of threatened and data deficient species') + xlab('') +
    geom_text(aes(y=PercTot+5, label=N)) +
 	scale_fill_manual(values=c('VU'='#fee08b', 'EN'='#f46d43', 'CR'='#a50026', 'DD'='#bababa')) +#'LC'='#1a9850', 'NT'='#a6d96a', 
    theme_classic()
dev.off()
