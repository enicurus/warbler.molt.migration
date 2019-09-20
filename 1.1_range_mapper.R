#obtain migratory distance
library(GISTools)
library(maptools)
library(raster)
library(geosphere)

shp.dir = '~/Dropbox/Warbler.Molt.Migration/Parulidae_shapefiles/'
map.dir = '~/Dropbox/Warbler.Molt.Migration/Parulidae_maps/'

files = list.files(shp.dir,pattern='*.shp')

taxa = as.character(sapply(files,function(x) paste(strsplit(x,'_')[[1]][1],'_',strsplit(x,'_')[[1]][2],sep='')))

wrld = readShapeSpatial("~/Documents/Glenn/Furnariidae/Furn_Maps/TM_WORLD_BORDERS-0.2/TM_WORLD_BORDERS-0.2.shp")
ex = c(-170,-30,-60,89) #new world
nwmap = crop(wrld,extent(ex))

i = 43
for(i in 1:length(files)){
	print(taxa[i])
	path = paste(shp.dir,files[i],sep='')
	shp = readShapeSpatial(path)
	projection(shp) = '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'
	
	pdf(paste(map.dir,taxa[i],'.pdf',sep=''))
	par(mar=c(2,2,2,2))
	plot(nwmap,lwd=0.1,border='white')

	for(f in unique(shp$SEASONAL)){
		if(f == 1){
			plot(shp[shp$SEASONAL == 1,],col='green',add=T,lwd=0.01)
		}else if(f == 2){
			plot(shp[shp$SEASONAL == 2,],col='orange',add=T,lwd=0.01)
		}else if(f == 3){
			plot(shp[shp$SEASONAL == 3,],col='blue',add=T,lwd=0.01)
		}else if(f == 4){
			#plot(shp[shp$SEASONAL == 4,],col='grey',add=T)
		}	
	}
	plot(nwmap,lwd=0.1,border='black',add=T)
	dev.off()	

}





