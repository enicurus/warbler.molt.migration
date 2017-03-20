#Combining all the data into one data frame for analysis

lats<-read.csv("/Users/ryanterrill/Dropbox/Warbler.Molt.Migration/Latitudes.csv")
clim<-read.csv("/Users/ryanterrill/Dropbox/Warbler.Molt.Migration/climate.csv")
hab<-read.csv("/Users/ryanterrill/Dropbox/Warbler.Molt.Migration/habitat_nest_stratum.csv")
mass<-read.csv("/Users/ryanterrill/Dropbox/Warbler.Molt.Migration/parulidMass.csv")
PA<-read.csv("/Users/ryanterrill/Dropbox/Warbler.Molt.Migration/Molt\ extent/Parulidae_PA.csv")
PF<-read.csv("/Users/ryanterrill/Dropbox/Warbler.Molt.Migration/Molt\ extent/Parulidae_PF.csv")
tempDimorph<-read.csv("/Users/ryanterrill/Dropbox/Warbler.Molt.Migration/Parulidae_temporal_dimorphism.csv")
sexDimorph<-read.csv("/Users/ryanterrill/Dropbox/Warbler.Molt.Migration/Parulidae_sexual_dimorphism.csv")
migDist<-read.table("/Users/ryanterrill/Dropbox/Warbler.Molt.Migration/migratory_distance.txt",header=TRUE)

#solar<-read.table("/Users/ryanterrill/Dropbox/Warbler.Molt.Migration/solar_daylight.txt",header=TRUE)##This is messed up, use daylength from Data file or redo


migration<-data.frame(migDist$migdist1)
rownames(migration)<-rownames(migDist)
names(migration)<-("migDist")

mass$mass<-as.numeric(gsub("-999",NA,mass$mass))
Mass<-mass$mass;names(Mass)<-mass$species

###have to split up the shapefile name

files<-as.character(solar$Shapefile)
taxa = as.character(sapply(files,function(x) paste(strsplit(x,'_')[[1]][1],'_',strsplit(x,'_')[[1]][2],sep='')))

sol<-solar$solar;names(sol)<-taxa
daylength<-solar$daylength;names(daylength)<-taxa
solar<-data.frame(sol,daylength)

rownames(sexDimorph)<-sexDimorph$tree_taxon
sexDimorph<-sexDimorph[2:13]
sexDimorph$extent<-rowSums(sexDimorph)
sdExtent<-sexDimorph$extent;names(sdExtent)<-rownames(sexDimorph)

rownames(tempDimorph)<-tempDimorph$tree_taxon
tempDimorph<-tempDimorph[2:13]
tempDimorph$extent<-rowSums(tempDimorph)
tdExtent<-tempDimorph$extent;names(tdExtent)<-rownames(tempDimorph)

rownames(PF)<-PF$tree_taxon
PF<-PF[,2:13]
PF$PFextent<-rowSums(PF)
pfextent<-PF$PFextent;names(pfextent)<-rownames(PF)


rownames(PA)<-PA$tree_taxon
PA<-PA[,2:13]
PA$PAextent<-rowSums(PA)
paextent<-PA$PAextent;names(paextent)<-rownames(PA)



habitat<-hab[,2:7];rownames(habitat)<-hab$species

climate<-clim[,2:11];rownames(climate)<-clim[,1]

latitude<-lats[,2:3];rownames(latitude)<-lats[,1]

habitat<-hab[,2:7];rownames(habitat)<-hab$species

data1<-data.frame(latitude,climate,solar,migration)###daylight got fucked - redo or get original values
data2<-data.frame(paextent,pfextent,Mass,sdExtent,tdExtent44)
data3<-merge(data1,data2,by=0,all=TRUE)




Data<-data3[,2:ncol(data3)];rownames(Data)<-data3[,1]

Data$absBreedLat<-abs(Data$breeding_lat)
Data$absWinLat<-abs(Data$winter_lat)
Data$SunExpos<-(Data$sol*Data$daylength)


pdf("~/Dropbox/Warbler.Molt.Migration/DataPairsplot.pdf")
pairs(Data,cex=.1,cex.labels=.1)
dev.off()


pdf("~/Dropbox/Warbler.Molt.Migration/corrPlot.pdf")
library(corrplot)
corMatrix<-cor(Data,use="pairwise.complete.obs")
corrplot(corMatrix,method="ellipse",sig.level=0.05)
dev.off()

Data<-Data[c("paextent","pfextent","tdExtent","sdExtent","daylength","migDist","absBreedLat","breeding_lat","absWinLat","sol","breedSolar","winSolar","meanSolar","winter_lat","precip_Breed","precip_Winter","minTemp_Breed","minTemp_Winter","maxTemp_Breed","maxTemp_Winter","avgTemp_Breed","avgTemp_Winter","elevation_Breed","elevation_Winter","Mass")]

write.csv(Data,"~/Dropbox/Warbler.Molt.Migration/Data.csv")

