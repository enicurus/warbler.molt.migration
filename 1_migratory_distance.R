#obtain migratory distance
library(GISTools)
library(maptools)
library(raster)
library(geosphere)

map.dir = '~/Dropbox/Warbler.Molt.Migration/Parulidae_shapefiles/'
files = list.files(map.dir,pattern='*.shp')

taxa = as.character(sapply(files,function(x) paste(strsplit(x,'_')[[1]][1],'_',strsplit(x,'_')[[1]][2],sep='')))

colnames = c('Shapefile','strategy','maxlat_breed','maxlat_winter','maxlat_resid','minlat_breed','minlat_winter','minlat_resid','midlat_breed','midlat_winter','midlat_resid','centroid_lat_breed','centroid_long_breed','centroid_lat_winter','centroid_long_winter','centroid_lat_resid','centroid_long_resid','migdist1','migdist2','migdist3','migdist4','migdist5','migdist6')

data = data.frame(matrix(nrow=length(taxa),ncol=length(colnames)))
colnames(data) = colnames
rownames(data) = taxa
data$Shapefile = files



# m 	= migratory - only breeding and winter ranges (2,3)
# nm 	= non-migratory - only resident ranges (1)
# mix 	= mixed - resident, breeding, and wintering ranges (1,2,3)
# par 	= partial - resident and breeding ranges (1,2) 

i = 48
for(i in 1:nrow(data)){
	print(rownames(data)[i])
	path = paste(map.dir,data$Shapefile[i],sep='')
	shp = readShapeSpatial(path)
	projection(shp) = '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'
	
	if(all(1 %in% shp@data$SEASONAL & !(c(2,3) %in% shp@data$SEASONAL))){
		
		data[i,'strategy'] = 'nm'
		
		res = shp[shp@data$SEASONAL == 1,]
		data[i,'maxlat_resid'] = extent(res)@ymax
		data[i,'minlat_resid'] = extent(res)@ymin
		data[i,'midlat_resid'] = extent(res)@ymax - ((extent(res)@ymax - extent(res)@ymin)/2)
		data[i,c('centroid_long_resid','centroid_lat_resid')] = gCentroid(res,byid=F)@coords
		
	}else if(all(c(2,3) %in% shp@data$SEASONAL & !(1 %in% shp@data$SEASONAL))){
		
		data[i,'strategy'] = 'm'

		breed = shp[shp@data$SEASONAL == 2,]
		data[i,'maxlat_breed'] = extent(breed)@ymax
		data[i,'minlat_breed'] = extent(breed)@ymin
		data[i,'midlat_breed'] = extent(breed)@ymax - ((extent(breed)@ymax - extent(breed)@ymin)/2)
		data[i,c('centroid_long_breed','centroid_lat_breed')] = gCentroid(breed,byid=F)@coords		
		
		winter = shp[shp@data$SEASONAL == 3,]
		data[i,'maxlat_winter'] = extent(winter)@ymax
		data[i,'minlat_winter'] = extent(winter)@ymin
		data[i,'midlat_winter'] = extent(winter)@ymax - ((extent(winter)@ymax - extent(winter)@ymin)/2)
		data[i,c('centroid_long_winter','centroid_lat_winter')] = gCentroid(winter,byid=F)@coords	
	
	}else if(all(c(1,2,3) %in% shp@data$SEASONAL)){
		
		data[i,'strategy'] = 'mix'
		
		res = shp[shp@data$SEASONAL == 1,]
		data[i,'maxlat_resid'] = extent(res)@ymax
		data[i,'minlat_resid'] = extent(res)@ymin
		data[i,'midlat_resid'] = extent(res)@ymax - ((extent(res)@ymax - extent(res)@ymin)/2)
		data[i,c('centroid_long_resid','centroid_lat_resid')] = gCentroid(res,byid=F)@coords
		
		breed = shp[shp@data$SEASONAL == 2,]
		data[i,'maxlat_breed'] = extent(breed)@ymax
		data[i,'minlat_breed'] = extent(breed)@ymin
		data[i,'midlat_breed'] = extent(breed)@ymax - ((extent(breed)@ymax - extent(breed)@ymin)/2)
		data[i,c('centroid_long_breed','centroid_lat_breed')] = gCentroid(breed,byid=F)@coords	
		
		winter = shp[shp@data$SEASONAL == 3,]
		data[i,'maxlat_winter'] = extent(winter)@ymax
		data[i,'minlat_winter'] = extent(winter)@ymin
		data[i,'midlat_winter'] = extent(winter)@ymax - ((extent(winter)@ymax - extent(winter)@ymin)/2)
		data[i,c('centroid_long_winter','centroid_lat_winter')] = gCentroid(winter,byid=F)@coords
		
		
	}else if(all(c(1,2) %in% shp@data$SEASONAL)){
		
		data[i,'strategy'] = 'par'
		
		res = shp[shp@data$SEASONAL == 1,]
		data[i,'maxlat_resid'] = extent(res)@ymax
		data[i,'minlat_resid'] = extent(res)@ymin
		data[i,'midlat_resid'] = extent(res)@ymax - ((extent(res)@ymax - extent(res)@ymin)/2)
		data[i,c('centroid_long_resid','centroid_lat_resid')] = gCentroid(res,byid=F)@coords
		
		breed = shp[shp@data$SEASONAL == 2,]
		data[i,'maxlat_breed'] = extent(breed)@ymax
		data[i,'minlat_breed'] = extent(breed)@ymin
		data[i,'midlat_breed'] = extent(breed)@ymax - ((extent(breed)@ymax - extent(breed)@ymin)/2)
		data[i,c('centroid_long_breed','centroid_lat_breed')] = gCentroid(breed,byid=F)@coords
	}	#end of if statement	
} #i loop

data

#calculate migratory distance

# migdist1 = midlat distance = distance between midlat_breed and midlat_winter (units = degrees latitude)
# migdist2 = maxlat distance = distance between maxlat_breed and maxlat_winter (units = degrees latitude)
# migdist3 = minlat distance = distance between minlat_breed and minlat_winter (units = degrees latitude)
# migdist4 = max distance = distance between maxlat_breed and minlat_winter (units = degrees latitude)
# migdist5 = min distance = distance between minlat_breed and maxlat_winter (units = degrees latitude)
# migdist6 = great circle distance between centroids (units = kilometers)



for(i in 1:nrow(data)){
	
	if(data[i,'strategy'] == 'nm'){
		
data[i,c('migdist1','migdist2','migdist3','migdist4','migdist5','migdist6')] = 0
	
	}else if(data[i,'strategy'] == 'm'){
		
		data[i,'migdist1'] = data[i,'midlat_breed'] - data[i,'midlat_winter']
		data[i,'migdist2'] = data[i,'maxlat_breed'] - data[i,'maxlat_winter']
		data[i,'migdist3'] = data[i,'minlat_breed'] - data[i,'minlat_winter']
		data[i,'migdist4'] = data[i,'maxlat_breed'] - data[i,'minlat_winter']
		data[i,'migdist5'] = data[i,'minlat_breed'] - data[i,'maxlat_winter']
		data[i,'migdist6'] = distHaversine(data[i,c('centroid_long_breed','centroid_lat_breed')], data[i,c('centroid_long_winter','centroid_lat_winter')],r=6378)
		
	}else if(data[i,'strategy'] == 'mix'){
		
		data[i,'migdist1'] = data[i,'midlat_breed'] - data[i,'midlat_winter']
		data[i,'migdist2'] = data[i,'maxlat_breed'] - data[i,'maxlat_winter']
		data[i,'migdist3'] = data[i,'minlat_breed'] - data[i,'minlat_winter']
		data[i,'migdist4'] = data[i,'maxlat_breed'] - data[i,'minlat_winter']
		data[i,'migdist5'] = data[i,'minlat_breed'] - data[i,'maxlat_winter']
		data[i,'migdist6'] = distHaversine(data[i,c('centroid_long_breed','centroid_lat_breed')], data[i,c('centroid_long_winter','centroid_lat_winter')],r=6378)

	}else if(data[i,'strategy'] == 'par'){
		
		data[i,'migdist1'] = data[i,'midlat_breed'] - data[i,'midlat_resid']
		data[i,'migdist2'] = data[i,'maxlat_breed'] - data[i,'maxlat_resid']
		data[i,'migdist3'] = data[i,'minlat_breed'] - data[i,'minlat_resid']
		data[i,'migdist4'] = data[i,'maxlat_breed'] - data[i,'minlat_resid']
		data[i,'migdist5'] = data[i,'minlat_breed'] - data[i,'maxlat_resid']
		data[i,'migdist6'] = distHaversine(data[i,c('centroid_long_breed','centroid_lat_breed')],data[i,c('centroid_long_resid','centroid_lat_resid')],r=6378)

	}	#end of if statement	

}#i loop

pairs(data[data$strategy %in% c('m','nm','par','mix'),c('migdist1','migdist2','migdist3','migdist4','migdist5','migdist6')])

pairs(data[data$strategy %in% c('m','nm','par'),c('migdist1','migdist2','migdist3','migdist4','migdist5','migdist6')])

pairs(data[data$strategy %in% c('m','nm'),c('migdist1','migdist2','migdist3','migdist4','migdist5','migdist6')])

pairs(data[data$strategy == 'm',c('migdist1','migdist2','migdist3','migdist4','migdist5','migdist6','PC1','PC2')])



pairs(data[data$strategy == 'mix' ,c('migdist1','migdist2','migdist3','migdist4','migdist5','migdist6','PC1','PC2')])



my.prc = prcomp(data[,c('migdist1','migdist2','migdist3','migdist4','migdist5','migdist6')],center=T,scale=T,na.omit=T)
pc.values = predict(my.prc)

importance = data.frame(summary(my.prc)$importance[2,])
colnames(importance) = 'proportion.of.variance'

str(my.prc)
round(my.prc$rotation,2)

fit = lm(data$PC1 ~ data$migdist1)
summary(fit)

data = cbind(data,pc.values[,1:2])





write.table(data,file = '~/Documents/Glenn/Molt/migratory_distance.txt',sep='\t',row.names=T,col.names=T,quote=F)



