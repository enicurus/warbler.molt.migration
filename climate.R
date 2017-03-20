#####This will project points onto winter and summer ranges for Parulidae , extract climate values at those
#####points from BIOCLIM layers, and then take the means of those values for breeding and winter months
##### for both summer and winter range by species
#####Maybe export some maps and values too


library(GISTools)
library(maptools)
library(raster)
library(geosphere)
library(sp)
library(rgdal)
library(wesanderson)



map.dir = '~/Dropbox/Warbler.Molt.Migration/Parulidae_shapefiles/'
files = list.files(map.dir,pattern='*.shp')

taxa = as.character(sapply(files,function(x) paste(strsplit(x,'_')[[1]][1],'_',strsplit(x,'_')[[1]][2],sep='')))


###making world map for plotting to make sure first few look right
world<-readShapeSpatial("~/Downloads/landsea_mask/landsea_mask.shp")
ex = c(-170,-30,-60,89) #new world
nwmap = crop(world,extent(ex))


####Import BioClim layers
BClim = getData("worldclim", var="bio", res=2.5)
BCLim = crop(BClim,ex)

####Import WorldClim layers
alt<-getData("worldclim",var="alt",res=2.5);alt<-crop(alt,extent(ex))
precip<-getData("worldclim",var="prec",res=2.5);precip<-crop(precip,extent(ex))
tmin<-getData("worldclim",var="tmin",res=2.5);tmin<-crop(tmin,extent(ex))
tmax<-getData("worldclim",var="tmax",res=2.5);tmax<-crop(tmax,extent(ex))
tavg<-getData("worldclim",var="tmean",res=2.5);tavg<-crop(tavg,extent(ex))

####Import solar radiation layers
Sol<-brick("/Users/ryanterrill/Dropbox/Warbler.Molt.Migration/solarRaster.grd")

nwmap<-readShapeSpatial("~/Dropbox/Warbler.Molt.Migration/nwmap.shp")
 

#Read BTNW - for some reason there is a problem with just this sp.

shapefile<-readShapeSpatial("~/Dropbox/Warbler.Molt.Migration/Parulidae_shapefiles/Dendroica_virens_9106_NS.shp")
shp1<-readShapeSpatial("~/Dropbox/Warbler.Molt.Migration/Parulidae_shapefiles/Dendroica_townsendi_9104_NS.shp")

###assign layers to precipitation layers to test
layers<-precip

###this function takes the mean of a climate value over 10,000 points within the range of a species
### migratory species split into breeding and winter range, points only taken from those months for mig and res species


clim_extract<-function(shapefile,layers){
	projection(shapefile) = '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'
	meanVals<-list()
	for(f in unique(shapefile$SEASONAL)){
		if(f==1){b_points<-spsample(shapefile[shapefile$SEASONAL==1,],10000,type="random")
				b_clim_points<-extract(layers,b_points)
				w_clim_points<-b_clim_points
				meanVals$breedMean<-mean(b_clim_points[,5:8],na.rm=TRUE)
				meanVals$winMean<-mean(w_clim_points[,c(1,2,11,12)],na.rm=TRUE)}
		else if(f==2){b_points<-spsample(shapefile[shapefile$SEASONAL==2,],10000,type="random")
				b_clim_points<-extract(layers,b_points)
				meanVals$breedMean<-mean(b_clim_points[,5:8],na.rm=TRUE)}
		else if(f==3){w_points<-spsample(shapefile[shapefile$SEASONAL==3,],10000,type="random")
				w_clim_points<-extract(layers,w_points)
				meanVals$winMean<-mean(w_clim_points[,c(1,2,11,12)],na.rm=TRUE)}
		}
		return(meanVals)
		}
		
		

alt_extract<-function(shapefile,layers){
	projection(shapefile) = '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'
	meanVals<-list()
	for(f in unique(shapefile$SEASONAL)){
		if(f==1){b_points<-spsample(shapefile[shapefile$SEASONAL==1,],10000,type="random")
				b_clim_points<-extract(layers,b_points)
				w_clim_points<-b_clim_points
				meanVals$breedMean<-mean(b_clim_points,na.rm=TRUE)
				meanVals$winMean<-mean(w_clim_points,na.rm=TRUE)}
		else if(f==2){b_points<-spsample(shapefile[shapefile$SEASONAL==2,],10000,type="random")
				b_clim_points<-extract(layers,b_points)
				meanVals$breedMean<-mean(b_clim_points,na.rm=TRUE)}
		else if(f==3){w_points<-spsample(shapefile[shapefile$SEASONAL==3,],10000,type="random")
				w_clim_points<-extract(layers,w_points)
				meanVals$winMean<-mean(w_clim_points,na.rm=TRUE)}
		}
		return(meanVals)
		}
			
	

####Ok great this works
####Now I want to set up a matrix for all species, with breeding and winter values for each climate layer
###So names in rownames, and pairs for each climate variable
###Then write a loop that automatically reads in each map and writes values to the matrix


taxa = as.character(sapply(files,function(x) paste(strsplit(x,'_')[[1]][1],'_',strsplit(x,'_')[[1]][2],sep='')))


climMatrix<-matrix(nrow=length(taxa),ncol=10)
rownames(climMatrix)<-taxa
colnames(climMatrix)<-c("precip_Breed","precip_Winter","minTemp_Breed","minTemp_Winter","maxTemp_Breed","maxTemp_Winter","avgTemp_Breed","avgTemp_Winter","elevation_Breed","elevation_Winter")

map.dir = '~/Dropbox/Warbler.Molt.Migration/Parulidae_shapefiles/'
data = data.frame(matrix(nrow=length(taxa),ncol=length(colnames)))
data$Shapefile = files


precipVals<-list()
for(i in 1:nrow(climMatrix)){
	path = paste(map.dir,data$Shapefile[i],sep='')
	shapefile = readShapeSpatial(path)
	cat("Species",i,"/",nrow(climMatrix),rownames(climMatrix)[i],"precipitation","\n")
	precipVals<-clim_extract(shapefile,precip)
	climMatrix[i,1]<-precipVals$breedMean
	climMatrix[i,2]<-precipVals$winMean
	}
	
	
for(i in 1:nrow(climMatrix)){
	path = paste(map.dir,data$Shapefile[i],sep='')
	shapefile = readShapeSpatial(path)
	cat("Species",i,"/",nrow(climMatrix),rownames(climMatrix)[i],"minTemp","\n")
	Vals<-clim_extract(shapefile,tmin)
	climMatrix[i,3]<-Vals$breedMean
	climMatrix[i,4]<-Vals$winMean
	}	
	
for(i in 1:nrow(climMatrix)){
	path = paste(map.dir,data$Shapefile[i],sep='')
	shapefile = readShapeSpatial(path)
	cat("Species",i,"/",nrow(climMatrix),rownames(climMatrix)[i],"maxTemp","\n")
	Vals<-clim_extract(shapefile,tmax)
	climMatrix[i,5]<-Vals$breedMean
	climMatrix[i,6]<-Vals$winMean
	}	
	
for(i in 1:nrow(climMatrix)){
	path = paste(map.dir,data$Shapefile[i],sep='')
	shapefile = readShapeSpatial(path)
	cat("Species",i,"/",nrow(climMatrix),rownames(climMatrix)[i],"avgTemp","\n")
	Vals<-clim_extract(shapefile,tavg)
	climMatrix[i,7]<-Vals$breedMean
	climMatrix[i,8]<-Vals$winMean
	}	
	
		
for(i in 1:nrow(climMatrix)){
	path = paste(map.dir,data$Shapefile[i],sep='')
	shapefile = readShapeSpatial(path)
	cat("Species",i,"/",nrow(climMatrix),rownames(climMatrix)[i],"alt","\n")
	Vals<-alt_extract(shapefile,alt)
	climMatrix[i,9]<-Vals$breedMean
	climMatrix[i,10]<-Vals$winMean
	}	
		
write.csv(climMatrix,"~/Dropbox/Warbler.Molt.Migration/climate.csv")		
	
btnw<-clim_extract(BTNW,precip)	


###Centroids for Breeding and Witner ranges:

migData<-read.table("/Users/ryanterrill/Dropbox/Warbler.Molt.Migration/migratory_distance.txt")

centroids<-data.frame(migData$centroid_lat_breed,migData$centroid_lat_winter,migData$centroid_lat_resid)
colnames(centroids)<-c("breeding_lat","winter_lat","resident")

for(i in 1:nrow(centroids)){
centroids[i,1][is.na(centroids[i,1])]<-centroids[i,3]
}
	
for(i in 1:nrow(centroids)){
centroids[i,2][is.na(centroids[i,2])]<-centroids[i,3]
}	
	
lats<-centroids[,1:2]	
head(l
write.csv(lats,"~/Dropbox/Warbler.Molt.Migration/latitudes.csv")



####plot examples of raster layers


wrld = readShapeSpatial("~/Documents/Glenn/Furnariidae/Furn_Maps/TM_WORLD_BORDERS-0.2/TM_WORLD_BORDERS-0.2.shp")

pal<-wes_palette("Zissou",100,type="continuous")

pdf('~/Dropbox/Warbler.Molt.Migration/solar_radiation.pdf')
par(mfrow=c(4,3),mar=c(2,2,2,2))
plot(nwsolar[['Jan']],main='Jan',col=pal,legend=FALSE);plot(nwmap,add=T,lwd=.1,border='black')
plot(nwsolar[['Feb']],main='Feb',col=pal,legend=FALSE);plot(nwmap,add=T,lwd=.1,border='black')
plot(nwsolar[['Mar']],main='Mar',col=pal,legend=FALSE);plot(nwmap,add=T,lwd=.1,border='black')
plot(nwsolar[['Apr']],main='Apr',col=pal,legend=FALSE);plot(nwmap,add=T,lwd=.1,border='black')
plot(nwsolar[['May']],main='May',col=pal,legend=FALSE);plot(nwmap,add=T,lwd=.1,border='black')
plot(nwsolar[['Jun']],main='Jun',col=pal,legend=FALSE);plot(nwmap,add=T,lwd=.1,border='black')
plot(nwsolar[['Jul']],main='Jul',col=pal,legend=FALSE);plot(nwmap,add=T,lwd=.1,border='black')
plot(nwsolar[['Aug']],main='Aug',col=pal,legend=FALSE);plot(nwmap,add=T,lwd=.1,border='black')
plot(nwsolar[['Sep']],main='Sep',col=pal,legend=FALSE);plot(nwmap,add=T,lwd=.1,border='black')
plot(nwsolar[['Oct']],main='Oct',col=pal,legend=FALSE);plot(nwmap,add=T,lwd=.1,border='black')
plot(nwsolar[['Nov']],main='Nov',col=pal,legend=FALSE);plot(nwmap,add=T,lwd=.1,border='black')
plot(nwsolar[['Dec']],main='Dec',col=pal,legend=FALSE);plot(nwmap,add=T,lwd=.1,border='black')
dev.off()


