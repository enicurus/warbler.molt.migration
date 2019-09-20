
#######################################################

TAKEN FROM

http://www.r-bloggers.com/maps-of-solar-radiation/

#######################################################

library(raster)
library(maps)
library(maptools)

wrld = readShapeSpatial("~/Documents/Glenn/Furnariidae/Furn_Maps/TM_WORLD_BORDERS-0.2/TM_WORLD_BORDERS-0.2.shp")

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

ex = c(-170,-30,-60,89) #new world
nwsolar = crop(grid,extent(ex))
nwmap = crop(wrld,extent(ex))

writeRaster(nwsolar,"~/Dropbox/Warbler.Molt.Migration/solarRaster")

library(wesanderson)
pal<-wes_palette("Zissou",100,type="continuous")

pdf('~/Dropbox/Warbler.Molt.Migration/solar_radiation.pdf')
par(mfrow=c(2,2),mar=c(2,2,2,2))
plot(nwsolar[['Jan']],main='Jan',col=pal)
plot(nwmap,add=T,lwd=.1,border='black')
plot(nwsolar[['Apr']],main='Apr',col=pal)
plot(nwmap,add=T,lwd=.1,border='black')
plot(nwsolar[['Jul']],main='Jul',col=pal)
plot(nwmap,add=T,lwd=.1,border='black')
plot(nwsolar[['Oct']],main='Oct',col=pal)
plot(nwmap,add=T,lwd=.1,border='black')
dev.off()

par(mfrow=c(1,1),mar=c(2,2,2,2))


x = shp
plot(x)
plot(rasters[[1]],add=T)
plot(x,add=T)

res1 = res
res2 = res

shp=res1
shp=res2
rasters=grid[[c(breeding,nonbreed)]]
#function to obtain average raster value for all raster cells in shapefile
mean.solar.for.shpfile.1 = function(shp,rasters){
	#extract weighted mean solar for each polygon for each month
	fff = extract(rasters,shp,weights=T,fun=mean)
	round(t(fff),7) == round(tmp2,7)
	tmp1 = extract(rasters[[1]],shp,weights=T,fun=mean)
	
	tmp1 = extract(rasters,shp,weights=T)
	if(tmp1)
	
	sum(tmp1[[2]][,'weight'],tmp1[[1]][,'weight'],tmp1[[3]][,'weight'])
	
	tmp2 = sapply(tmp1, function(x) if (!is.null(x)) {
				x = tmp1[[1]]
				apply(x[,-ncol(x)],2,function(y){ 
					y = x[,1]
					weighted.mean(y,x[,'weight']) 
					})
				} else NA  )
	
	#calculate mean solar radiation/month for each polygon in shapefile
	tmp2 = apply(tmp1,1,FUN=mean)
	#calculate mean solar radiation/month for shapefile weighted by area of each polygon
	areas = sapply(slot(shp,'polygons'),slot,'area')
	return(weighted.mean(tmp2,areas))	
}

mean.solar.for.shpfile.2 = function(shp,rasters){
	
	#extract weighted mean solar for each polygon for each month
	xx = c()
	for(f in 1:dim(rasters)[3]){
		tmp1 =	extract(rasters[[f]],shp,weights=T,fun=mean)
		#calculate mean solar radiation/month for each polygon in shapefile
		xx = rbind(xx,tmp1)			
	}
	
	tmp2 = apply(xx,2,mean)
	
	#calculate mean solar radiation/month for shapefile weighted by area of each polygon
	areas = sapply(slot(shp,'polygons'),slot,'area')
	return(weighted.mean(tmp2,areas))	
}



#obtain migratory distance
library(GISTools)
library(maptools)
library(raster)
library(geosphere)

map.dir = '~/Documents/Glenn/Molt/Parulidae_maps/'
files = list.files(map.dir,pattern='*.shp')

taxa = as.character(sapply(files,function(x) paste(strsplit(x,'_')[[1]][1],'_',strsplit(x,'_')[[1]][2],sep='')))

colnames = c('Shapefile','strategy','solar.1','solar.2')

data = data.frame(matrix(nrow=length(taxa),ncol=length(colnames)))
colnames(data) = colnames
rownames(data) = taxa
data$Shapefile = files

# m 	= migratory - only wintering and breeding ranges (2,3)
# nm 	= non-migratory - only resident ranges (1)
# mix 	= mixed - resident, wintering and breeding ranges (1,2,3)
# par 	= partial - resident and breeding ranges (1,2) 


breeding = c('May','Jun','Jul','Aug')
nonbreed = c('Jan','Feb','Nov','Dec') 
i = 19
for(i in 2:10){
	print(rownames(data)[i])
	path = paste(map.dir,data$Shapefile[i],sep='')
	shp = readShapeSpatial(path)
	projection(shp) = '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'
	
	
	#non-migratory - only resident ranges (1)
	if(all(1 %in% shp@data$SEASONAL & !(c(2,3) %in% shp@data$SEASONAL))){
		
		data[i,'strategy'] = 'nm'
			
		res = shp[shp@data$SEASONAL == 1, ]
		sol.resid.1 = mean.solar.for.shpfile.1(shp=res,rasters=grid[[c(breeding,nonbreed)]])
		sol.resid.2 = mean.solar.for.shpfile.2(shp=res,rasters=grid[[c(breeding,nonbreed)]])

		data[i,'solar.1'] = sol.resid.1
		data[i,'solar.2'] = sol.resid.2

		
	#migratory - only breeding and wintering ranges (2,3)	
	}else if(all(c(2,3) %in% shp@data$SEASONAL & !(1 %in% shp@data$SEASONAL))){
		
		data[i,'strategy'] = 'm'

		breed = shp[shp@data$SEASONAL == 2,]
		sol.breed.1 = mean.solar.for.shpfile.1(shp=breed,rasters=grid[[breeding]])
		sol.breed.2 = mean.solar.for.shpfile.2(shp=breed,rasters=grid[[breeding]])

		winter = shp[shp@data$SEASONAL == 3,]
		sol.winter.1 = mean.solar.for.shpfile.1(shp=winter,rasters=grid[[nonbreed]])
		sol.winter.2 = mean.solar.for.shpfile.2(shp=winter,rasters=grid[[nonbreed]])

		data[i,'solar.1'] = mean(sol.breed.1,sol.winter.1)
		data[i,'solar.2'] = mean(sol.breed.2,sol.winter.2)	
	
	#mixed - resident, breeding and wintering ranges (1,2,3)
	}else if(all(c(1,2,3) %in% shp@data$SEASONAL)){
		
		data[i,'strategy'] = 'mix'
		
		resid = shp[shp@data$SEASONAL == 1,]
		sol.resid.1 = mean.solar.for.shpfile.1(shp=res,rasters=grid[[c(breeding,nonbreed)]])
		sol.resid.2 = mean.solar.for.shpfile.2(shp=res,rasters=grid[[c(breeding,nonbreed)]])
		
		breed = shp[shp@data$SEASONAL == 2,]
		sol.breed.1 = mean.solar.for.shpfile.1(shp=breed,rasters=grid[[breeding]])
		sol.breed.2 = mean.solar.for.shpfile.2(shp=breed,rasters=grid[[breeding]])

		winter = shp[shp@data$SEASONAL == 3,]
		sol.winter.1 = mean.solar.for.shpfile.1(shp=winter,rasters=grid[[nonbreed]])
		sol.winter.2 = mean.solar.for.shpfile.2(shp=winter,rasters=grid[[nonbreed]])

		data[i,'solar.1'] = mean(sol.resid.1,sol.breed.1,sol.winter.1)
		data[i,'solar.2'] = mean(sol.resid.2,sol.breed.2,sol.winter.2)

		
	}else if(all(c(1,2) %in% shp@data$SEASONAL)){
		
		data[i,'strategy'] = 'par'
		
		resid = shp[shp@data$SEASONAL == 1,]
		sol.resid.1 = mean.solar.for.shpfile.1(shp=res,rasters=grid[[c(breeding,nonbreed)]])
		sol.resid.2 = mean.solar.for.shpfile.2(shp=res,rasters=grid[[c(breeding,nonbreed)]])

		breed = shp[shp@data$SEASONAL == 2,]
		sol.breed.1 = mean.solar.for.shpfile.1(shp=breed,rasters=grid[[breeding]])
		sol.breed.2 = mean.solar.for.shpfile.2(shp=breed,rasters=grid[[breeding]])

		data[i,'solar.1'] = mean(sol.resid.1,sol.breed.1)
		data[i,'solar.2'] = mean(sol.resid.2,sol.breed.2)

		
	}	#end of if statement	
} #i loop

data[i,]



