library(raster)
library(maps)
library(maptools)
library(sp)
library(GISTools)
library(geosphere)
##################################################################
#Import solar radiation data

nasa = read.table(file='~/Dropbox/Warbler.Molt.Migration/global_radiation.txt', skip=13, header=TRUE)
coordinates(nasa) = ~ Lon + Lat
proj4string(nasa) = CRS('+proj=longlat +ellps=WGS84')
summary(nasa)
extent = extent(nasa)
length(-180:179)
length(-90:89)
rast = raster(extent,ncol=360,nrow=180,crs='+proj=longlat +ellps=WGS84')
grid = rasterize(nasa,rast,fun='last')
grid = grid[[2:13]]


##################################################################
#Plot the solar radiation data

wrld = readShapeSpatial("~/Documents/Glenn/Furnariidae/Furn_Maps/TM_WORLD_BORDERS-0.2/TM_WORLD_BORDERS-0.2.shp")

ex = c(-170,-30,-60,89) #new world
nwsolar = crop(grid,extent(ex))
nwmap = crop(wrld,extent(ex))

pdf('~/Dropbox/Warbler.Molt.Migration/solar_radiation.pdf')
par(mfrow=c(2,2),mar=c(2,2,2,2))
plot(nwsolar[['Jan']],main='Jan')
plot(nwmap,add=T,lwd=.1,border='black')
plot(nwsolar[['Apr']],main='Apr')
plot(nwmap,add=T,lwd=.1,border='black')
plot(nwsolar[['Jul']],main='Jul')
plot(nwmap,add=T,lwd=.1,border='black')
plot(nwsolar[['Oct']],main='Oct')
plot(nwmap,add=T,lwd=.1,border='black')
dev.off()

par(mfrow=c(1,1),mar=c(2,2,2,2))

##################################################################
#Import average monthly daylight hours
day = read.table(file='~/Dropbox/Warbler.Molt.Migration/daylight.txt', skip=8, header=TRUE)


##################################################################
##################################################################

map.dir = '~/Dropbox/Warbler.Molt.Migration/Parulidae_shapefiles/'
files = list.files(map.dir,pattern='*.shp')

taxa = as.character(sapply(files,function(x) paste(strsplit(x,'_')[[1]][1],'_',strsplit(x,'_')[[1]][2],sep='')))

colnames = c('Shapefile','strategy','solar','daylight')

data = data.frame(matrix(nrow=length(taxa),ncol=length(colnames)))
colnames(data) = colnames
rownames(data) = taxa
data$Shapefile = files

# m 	= migratory - only wintering and breeding ranges (2,3)
# nm 	= non-migratory - only resident ranges (1)
# mix 	= mixed - resident, wintering and breeding ranges (1,2,3)
# par 	= partial - resident and breeding ranges (1,2) 


breeding = c('May','Jun','Jul','Aug')
nonbreed = c('Jan','Feb','Mar','Apr','Sep','Oct','Nov','Dec') 

head(data)
i = 29
i = 48
for(i in 1:nrow(data)){
	print(rownames(data)[i])
	path = paste(map.dir,data$Shapefile[i],sep='')
	shp = readShapeSpatial(path)
	projection(shp) = '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'
	
	sol.tmp = matrix(nrow=3,ncol=13)
	colnames(sol.tmp) = c('Area','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')
	
	day.tmp = matrix(nrow=3,ncol=13)
	colnames(day.tmp) = c('Area','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec') 

	#non-migratory - only resident ranges (1)
	if(all(1 %in% shp@data$SEASONAL & !(c(2,3) %in% shp@data$SEASONAL))){
		
		data[i,'strategy'] = 'nm'
		
		res = shp[shp@data$SEASONAL == 1, ]
		months = c(breeding,nonbreed)
		samp = spsample(res,1000,type='random')
		area = mean(sapply(slot(res,'polygons'),slot,'area'))
		sol.tmp[1,'Area'] = area
		day.tmp[1,'Area'] = area
		sol.resid = extract(grid[[months]],samp) 
		sol.tmp[1,months] = apply(sol.resid,2,FUN=mean)
		#for each latitude coordinate, round to nearest integer then extract row in data.frame 'day'
		#that it corresponds to, then average across all rows and columns
		day.resid=t(sapply(round(samp@coords[,2],0),function(x){
				as.numeric(day[day$Lat == x, months])	
				}))
		day.tmp[1,months] = apply(day.resid,2,FUN=mean)

		
		
		data[i,'solar'] = mean(apply(sol.tmp[,-1],2,FUN=weighted.mean,w=sol.tmp[,1],na.rm=T))
		data[i,'daylight'] = mean(apply(day.tmp[,-1],2,FUN=weighted.mean,w=day.tmp[,1],na.rm=T))
						
	#migratory - only breeding and wintering ranges (2,3)	
	}else if(all(c(2,3) %in% shp@data$SEASONAL & !(1 %in% shp@data$SEASONAL))){
		
		data[i,'strategy'] = 'm'

		breed = shp[shp@data$SEASONAL == 2,]
		months = c(breeding)
		samp = spsample(breed,1000,type='random')
		area = mean(sapply(slot(breed,'polygons'),slot,'area'))
		sol.tmp[2,'Area'] = area
		day.tmp[2,'Area'] = area
		sol.breed = extract(grid[[months]],samp) 
		sol.tmp[2,months] = apply(sol.breed,2,FUN=mean)
		day.breed=t(sapply(round(samp@coords[,2],0),function(x){
				as.numeric(day[day$Lat == x, months])	
				}))
		day.tmp[2,months] = apply(day.breed,2,FUN=mean)


		winter = shp[shp@data$SEASONAL == 3,]
		months = c(nonbreed)
		samp = spsample(winter,1000,type='random')
		area = mean(sapply(slot(winter,'polygons'),slot,'area'))
		sol.tmp[3,'Area'] = area
		day.tmp[3,'Area'] = area
		sol.winter = extract(grid[[months]],samp) 
		sol.tmp[3,months] = apply(sol.winter,2,FUN=mean)
		day.winter=t(sapply(round(samp@coords[,2],0),function(x){
				as.numeric(day[day$Lat == x, months])	
				}))
		day.tmp[3,months] = apply(day.winter,2,FUN=mean)


		
		data[i,'solar'] = mean(apply(sol.tmp[,-1],2,FUN=weighted.mean,w=sol.tmp[,1],na.rm=T))
		data[i,'daylight'] = mean(apply(day.tmp[,-1],2,FUN=weighted.mean,w=day.tmp[,1],na.rm=T))

	#mixed - resident, breeding and wintering ranges (1,2,3)
	}else if(all(c(1,2,3) %in% shp@data$SEASONAL)){
		
		data[i,'strategy'] = 'mix'

		res = shp[shp@data$SEASONAL == 1, ]
		months = c(breeding,nonbreed)
		samp = spsample(res,1000,type='random')
		area = mean(sapply(slot(res,'polygons'),slot,'area'))
		sol.tmp[1,'Area'] = area
		day.tmp[1,'Area'] = area
		sol.resid = extract(grid[[months]],samp) 
		sol.tmp[1,months] = apply(sol.resid,2,FUN=mean)
		day.resid=t(sapply(round(samp@coords[,2],0),function(x){
				as.numeric(day[day$Lat == x, months])	
				}))
		day.tmp[1,months] = apply(day.resid,2,FUN=mean)


		breed = shp[shp@data$SEASONAL == 2,]
		months = c(breeding)
		samp = spsample(breed,1000,type='random')
		area = mean(sapply(slot(breed,'polygons'),slot,'area'))
		sol.tmp[2,'Area'] = area
		day.tmp[2,'Area'] = area
		sol.breed = extract(grid[[months]],samp) 
		sol.tmp[2,months] = apply(sol.breed,2,FUN=mean)
		day.breed=t(sapply(round(samp@coords[,2],0),function(x){
				as.numeric(day[day$Lat == x, months])	
				}))
		day.tmp[2,months] = apply(day.breed,2,FUN=mean)


		winter = shp[shp@data$SEASONAL == 3,]
		months = c(nonbreed)
		samp = spsample(winter,1000,type='random')
		area = mean(sapply(slot(winter,'polygons'),slot,'area'))
		sol.tmp[3,'Area'] = area
		day.tmp[3,'Area'] = area
		sol.winter = extract(grid[[months]],samp) 
		sol.tmp[3,months] = apply(sol.winter,2,FUN=mean)
		day.winter=t(sapply(round(samp@coords[,2],0),function(x){
				as.numeric(day[day$Lat == x, months])	
				}))
		day.tmp[3,months] = apply(day.winter,2,FUN=mean)


		data[i,'solar'] = mean(apply(sol.tmp[,-1],2,FUN=weighted.mean,w=sol.tmp[,1],na.rm=T))
		data[i,'daylight'] = mean(apply(day.tmp[,-1],2,FUN=weighted.mean,w=day.tmp[,1],na.rm=T))

	#partial - resident and breeding ranges (1,2)
	}else if(all(c(1,2) %in% shp@data$SEASONAL)){
		
		data[i,'strategy'] = 'par'
		
		res = shp[shp@data$SEASONAL == 1, ]
		months = c(breeding,nonbreed)
		samp = spsample(res,1000,type='random')
		area = mean(sapply(slot(res,'polygons'),slot,'area'))
		sol.tmp[1,'Area'] = area
		day.tmp[1,'Area'] = area
		sol.resid = extract(grid[[months]],samp) 
		sol.tmp[1,months] = apply(sol.resid,2,FUN=mean)
		day.resid=t(sapply(round(samp@coords[,2],0),function(x){
				as.numeric(day[day$Lat == x, months])	
				}))
		day.tmp[1,months] = apply(day.resid,2,FUN=mean)

		
		breed = shp[shp@data$SEASONAL == 2,]
		months = c(breeding)
		samp = spsample(breed,1000,type='random')
		area = mean(sapply(slot(breed,'polygons'),slot,'area'))
		sol.tmp[2,'Area'] = area
		day.tmp[2,'Area'] = area
		sol.breed = extract(grid[[months]],samp) 
		sol.tmp[2,months] = apply(sol.breed,2,FUN=mean)
		day.breed=t(sapply(round(samp@coords[,2],0),function(x){
				as.numeric(day[day$Lat == x, months])	
				}))
		day.tmp[2,months] = apply(day.breed,2,FUN=mean)

		data[i,'solar'] = mean(apply(sol.tmp[,-1],2,FUN=weighted.mean,w=sol.tmp[,1],na.rm=T))
		data[i,'daylight'] = mean(apply(day.tmp[,-1],2,FUN=weighted.mean,w=day.tmp[,1],na.rm=T))
		
	}	#end of if statement	
} #i loop

data

write.table(data,'~/Dropbox/Warbler.Molt.Migration/solar_daylight.txt',quote=F,row.names=F,col.names=T,sep='\t')

plot(data$solar,data$daylight)

cer = day.tmp
tri = day.tmp
blk = day.tmp

mean(cer[2,breeding])
mean(tri[1,breeding])
mean(blk[2,breeding])

mean(cer[3,nonbreed])
mean(tri[1,nonbreed])
mean(blk[3,nonbreed])



