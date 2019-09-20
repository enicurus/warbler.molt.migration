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

wrld = readShapeSpatial("~/Dropbox/Warbler.Molt.Migration/nwmap.shp ")

ex = c(-170,-30,-60,89) #new world
nwsolar = crop(grid,extent(ex))
nwmap = crop(wrld,extent(ex))

#pdf('~/Dropbox/Warbler.Molt.Migration/solar_radiation_quad.pdf')
#par(mfrow=c(2,2),mar=c(2,2,2,2))
#plot(nwsolar[['Jan']],main='Jan')
#plot(nwmap,add=T,lwd=.1,border='black')
#plot(nwsolar[['Apr']],main='Apr')
#plot(nwmap,add=T,lwd=.1,border='black')
#plot(nwsolar[['Jul']],main='Jul')
#plot(nwmap,add=T,lwd=.1,border='black')
#plot(nwsolar[['Oct']],main='Oct')
#plot(nwmap,add=T,lwd=.1,border='black')
#dev.off()

##################################################################
#Import average monthly daylight hours
day = read.table(file='~/Dropbox/Warbler.Molt.Migration/daylight.txt', skip=8, header=TRUE)


##################################################################
##################################################################

map.dir = '~/Dropbox/Warbler.Molt.Migration/Parulidae_shapefiles/'
files = list.files(map.dir,pattern='*.shp')

taxa = as.character(sapply(files,function(x) paste(strsplit(x,'_')[[1]][1],'_',strsplit(x,'_')[[1]][2],sep='')))

colnames = c('Shapefile','strategy','solar','daylength')

data = data.frame(matrix(nrow=length(taxa),ncol=length(colnames)))
colnames(data) = colnames
rownames(data) = taxa
data$Shapefile = files

# m 	= migratory - only breeding and wintering ranges (2,3)
# nm 	= non-migratory - only resident ranges (1)
# mix 	= mixed - resident, breeding and wintering ranges (1,2,3)
# par 	= partial - resident and breeding ranges (1,2) 


breeding = c('May','Jun','Jul','Aug')
nonbreed = c('Jan','Feb','Nov','Dec') 


for(i in 1:nrow(data)){
	print(rownames(data)[i])
	path = paste(map.dir,data$Shapefile[i],sep='')
	shp = readShapeSpatial(path)
	projection(shp) = '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'
	
	polys = sort(unique(shp@data$SEASONAL))
	polys = as.numeric(polys[polys <= 3])
	
	if(rownames(data)[i] == 'Dendroica_petechia'){
		polys= c(2,3)
	}
	
	if(identical(c(1),polys)){
		data[i,'strategy'] = 'nm'
	}else if(identical(c(1,2),polys)){
		data[i,'strategy'] = 'par'	
	}else if(identical(c(1,2,3),polys)){
		data[i,'strategy'] = 'mix'	
	}else if(identical(c(2,3),polys)){
		data[i,'strategy'] = 'm'	
	}
	
	sol.tmp = matrix(nrow=3,ncol=13)
	colnames(sol.tmp) = c('Area','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')
	
	day.tmp = matrix(nrow=3,ncol=13)
	colnames(day.tmp) = c('Area','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec') 
		
	for(p in polys){
		poly = shp[shp@data$SEASONAL == p, ]
		samp = spsample(poly,1000,type='random')
		area = mean(sapply(slot(poly,'polygons'),slot,'area'))
		sol.tmp[p,'Area'] = day.tmp[p,'Area'] = area
		
		sol.poly = extract(grid,samp) 
		sol.tmp[p,2:13] = apply(sol.poly,2,FUN=mean)
		
		day.poly=t(sapply(round(samp@coords[,2],0),function(x){
				as.numeric(day[day$Lat == x, ]) 	
				}))
		day.tmp[p,2:13] = apply(day.poly[,2:13],2,FUN=mean)
	} # p loop
	
	
	day.tmp[2,nonbreed] = sol.tmp[2,nonbreed] = NA
	day.tmp[3,breeding] = sol.tmp[3,breeding] = NA
	
	data[i,'solar']		= mean(apply(sol.tmp[,-1],2,FUN=weighted.mean,w=sol.tmp[,1],na.rm=T))
	data[i,'daylength'] = mean(apply(day.tmp[,-1],2,FUN=weighted.mean,w=day.tmp[,1],na.rm=T))
		
}


write.table(data,'~/Dropbox/Warbler.Molt.Migration/solar_daylight.txt',quote=F,row.names=F,col.names=T,sep='\t')



