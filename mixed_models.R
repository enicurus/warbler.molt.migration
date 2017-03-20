####

Mixed models

###PGLS, AIC weights for all likely models for the four predictor variables
###Mixed Models built from examining PGLS response of individual variables


library(MuMIn)
library(MCMCglmm)

##Maybe just use PGLS and design models by hand
## Use correlation matrix to build models

##################################################################
###PA extent models
##################################################################


source("/Users/ryanterrill/Dropbox/Warbler.Molt.Migration/model_selection.R")

pat<-treedata(warTree,na.omit((PA[,1,drop=FALSE])))
PAtree<-pat$phy


Data<-read.csv("~/Dropbox/Warbler.Molt.Migration/Data_treeNames.csv",row.names=1)
Data$names<-rownames(Data)

###BM is the supported model for PA
PA.dat<-comparative.data(warTree,Data,names,vcv=TRUE,vcv.dim=3,warn.dropped=TRUE)


## PA extent predictor mixed models
PAmods<-list()
PAmods[[1]]<-pgls(data=PA.dat,paextent~daylength)	
PAmods[[2]]<-pgls(data=PA.dat,paextent~migDist)	
PAmods[[3]]<-pgls(data=PA.dat,paextent~breeding_forage_stratum)	
PAmods[[4]]<-pgls(data=PA.dat,paextent~winter_strat)
PAmods[[5]]<-pgls(data=PA.dat,paextent~breedSolar)	
PAmods[[6]]<-pgls(data=PA.dat,paextent~precip_Breed)
PAmods[[7]]<-pgls(data=PA.dat,paextent~minTemp_Breed)
PAmods[[8]]<-pgls(data=PA.dat,paextent~avgTemp_Breed)
PAmods[[9]]<-pgls(data=PA.dat,paextent~Mass)
PAmods[[10]]<-pgls(data=PA.dat,paextent~daylength*migDist)	
PAmods[[11]]<-pgls(data=PA.dat,paextent~daylength*breedSolar)	
PAmods[[12]]<-pgls(data=PA.dat,paextent~daylength*breeding_forage_stratum)	
PAmods[[13]]<-pgls(data=PA.dat,paextent~daylength*winter_strat)	
PAmods[[14]]<-pgls(data=PA.dat,paextent~daylength*precip_Breed)	
PAmods[[15]]<-pgls(data=PA.dat,paextent~daylength*minTemp_Breed)	
PAmods[[16]]<-pgls(data=PA.dat,paextent~daylength*breeding_forage_stratum*winter_strat)	
PAmods[[17]]<-pgls(data=PA.dat,paextent~daylength*breeding_forage_stratum*winter_strat*breedSolar)	
PAmods[[18]]<-pgls(data=PA.dat,paextent~breeding_forage_stratum*winter_strat*breedSolar)	
PAmods[[19]]<-pgls(data=PA.dat,paextent~daylength*breeding_forage_stratum*winter_strat*Mass)	
PAmods[[20]]<-pgls(data=PA.dat,paextent~daylength*breedSolar*sol)
PAmods[[21]]<-pgls(data=PA.dat,paextent~daylength*sol)
PAmods[[22]]<-pgls(data=PA.dat,paextent~daylength*winSolar)
PAmods[[23]]<-pgls(data=PA.dat,paextent~migDist*breedSolar)
PAmods[[24]]<-pgls(data=PA.dat,paextent~migDist*sol)
PAmods[[25]]<-pgls(data=PA.dat,paextent~migDist*winSolar)
PAmods[[26]]<-pgls(data=PA.dat,paextent~migDist*daylength*breedSolar)

PAmodsAIC<-matrix(nrow=26,ncol=4)
for (i in 1:26){
	PAmodsAIC[i,2]<-PAmods[[i]]$aicc
	PAmodsAIC[i,3]<-summary(PAmods[[i]])$adj.r.squared
	PAmodsAIC[i,4]<-pf(summary(PAmods[[i]])$fstatistic[1],summary(PAmods[[i]])$fstatistic[2],summary(PAmods[[i]])$fstatistic[3],lower.tail=F)
	PAmodsAIC[i,1]<-paste("PA_Model",i,sep="_")
		}
		
PAmodsAIC<-data.frame(PAmodsAIC)
PAmodsAIC[,2]<-as.numeric(as.character(PAmodsAIC[,2]))
PAmodsAIC[,3]<-as.numeric(as.character(PAmodsAIC[,3]))


names(PAmodsAIC)<-c("model","AICc","adj_r_sq","p-value")

paM<-PAmodsAIC[order(PAmodsAIC$AICc),]
write.csv(paM,"~/Dropbox/Warbler.Molt.Migration/PA_models.csv")
###this all works and makes sense. Consider rerunning some of the top models by redion or modular unit - ie does solar radiation affect top of bird more than bottom?

##################################################################
### Preformative molt extent models
##################################################################


PFmods<-list()
PFmods[[1]]<-pgls(data=PA.dat,pfextent~sdExtent)
PFmods[[2]]<-pgls(data=PA.dat,pfextent~daylength)
PFmods[[3]]<-pgls(data=PA.dat,pfextent~absWinLat)
PFmods[[4]]<-pgls(data=PA.dat,pfextent~breeding_forage_stratum)
PFmods[[5]]<-pgls(data=PA.dat,pfextent~nest_site_stratum)
PFmods[[6]]<-pgls(data=PA.dat,pfextent~winter_strat)
PFmods[[7]]<-pgls(data=PA.dat,pfextent~breedSolar)
PFmods[[8]]<-pgls(data=PA.dat,pfextent~SunExpos)
PFmods[[9]]<-pgls(data=PA.dat,pfextent~precip_Breed)
PFmods[[10]]<-pgls(data=PA.dat,pfextent~minTemp_Breed)
PFmods[[11]]<-pgls(data=PA.dat,pfextent~avgTemp_Breed)
PFmods[[12]]<-pgls(data=PA.dat,pfextent~maxTemp_Breed)
PFmods[[13]]<-pgls(data=PA.dat,pfextent~sdExtent*daylength)
PFmods[[14]]<-pgls(data=PA.dat,pfextent~daylength*absWinLat)
PFmods[[15]]<-pgls(data=PA.dat,pfextent~daylength*breeding_forage_stratum)
PFmods[[16]]<-pgls(data=PA.dat,pfextent~daylength*nest_site_stratum)
PFmods[[17]]<-pgls(data=PA.dat,pfextent~daylength*breedSolar)
PFmods[[18]]<-pgls(data=PA.dat,pfextent~daylength*sol)
PFmods[[19]]<-pgls(data=PA.dat,pfextent~breeding_forage_stratum*breedSolar)
PFmods[[20]]<-pgls(data=PA.dat,pfextent~daylength*nest_site_stratum)


PFmodsAIC<-matrix(nrow=20,ncol=4)
for (i in 1:20){
	PFmodsAIC[i,2]<-PFmods[[i]]$aicc
	PFmodsAIC[i,3]<-summary(PFmods[[i]])$adj.r.squared
	PFmodsAIC[i,4]<-pf(summary(PFmods[[i]])$fstatistic[1],summary(PFmods[[i]])$fstatistic[2],summary(PFmods[[i]])$fstatistic[3],lower.tail=F)
	PFmodsAIC[i,1]<-paste("PF_Model",i,sep="_")
		}
		
PFmodsAIC<-data.frame(PFmodsAIC)
PFmodsAIC[,2]<-as.numeric(as.character(PFmodsAIC[,2]))
PFmodsAIC[,3]<-as.numeric(as.character(PFmodsAIC[,3]))


names(PFmodsAIC)<-c("model","AICc","adj_r_sq","p-value")

PFM<-PFmodsAIC[order(PFmodsAIC$AICc),]
write.csv(PFM,"~/Dropbox/Warbler.Molt.Migration/PF_models.csv")

##################################################################
##Sexual Dichromatism extent models
##################################################################


SDmods<-list()

SDmods[[1]]<-pgls(data=PA.dat,sdExtent~daylength)
SDmods[[2]]<-pgls(data=PA.dat,sdExtent~migDist)
SDmods[[3]]<-pgls(data=PA.dat,sdExtent~absBreedLat)
SDmods[[4]]<-pgls(data=PA.dat,sdExtent~absWinLat)
SDmods[[5]]<-pgls(data=PA.dat,sdExtent~breeding_forage_stratum)
SDmods[[6]]<-pgls(data=PA.dat,sdExtent~nest_site_stratum)
SDmods[[7]]<-pgls(data=PA.dat,sdExtent~winter_hab)
SDmods[[8]]<-pgls(data=PA.dat,sdExtent~sol)
SDmods[[9]]<-pgls(data=PA.dat,sdExtent~breedSolar)
SDmods[[10]]<-pgls(data=PA.dat,sdExtent~SunExpos)
SDmods[[11]]<-pgls(data=PA.dat,sdExtent~minTemp_Breed)
SDmods[[12]]<-pgls(data=PA.dat,sdExtent~maxTemp_Breed)
SDmods[[13]]<-pgls(data=PA.dat,sdExtent~avgTemp_Breed)
SDmods[[14]]<-pgls(data=PA.dat,sdExtent~Mass)
SDmods[[15]]<-pgls(data=PA.dat,sdExtent~daylength*breeding_forage_stratum)
SDmods[[16]]<-pgls(data=PA.dat,sdExtent~daylength*nest_site_stratum)
SDmods[[17]]<-pgls(data=PA.dat,sdExtent~daylength*breeding_forage_stratum)
SDmods[[18]]<-pgls(data=PA.dat,sdExtent~absBreedLat*breeding_forage_stratum)
SDmods[[19]]<-pgls(data=PA.dat,sdExtent~absBreedLat*nest_site_stratum)
SDmods[[20]]<-pgls(data=PA.dat,sdExtent~winter_hab*nest_site_stratum)
SDmods[[21]]<-pgls(data=PA.dat,sdExtent~absBreedLat*Mass)
SDmods[[22]]<-pgls(data=PA.dat,sdExtent~absWinLat*breeding_forage_stratum)


SDmodsAIC<-matrix(nrow=22,ncol=4)
for (i in 1:22){
	SDmodsAIC[i,2]<-SDmods[[i]]$aicc
	SDmodsAIC[i,3]<-summary(SDmods[[i]])$adj.r.squared
	SDmodsAIC[i,4]<-pf(summary(SDmods[[i]])$fstatistic[1],summary(SDmods[[i]])$fstatistic[2],summary(SDmods[[i]])$fstatistic[3],lower.tail=F)
	SDmodsAIC[i,1]<-paste("SD_Model",i,sep="_")
		}
		
SDmodsAIC<-data.frame(SDmodsAIC)
SDmodsAIC[,2]<-as.numeric(as.character(SDmodsAIC[,2]))
SDmodsAIC[,3]<-as.numeric(as.character(SDmodsAIC[,3]))


names(SDmodsAIC)<-c("model","AICc","adj_r_sq","p-value")

SDmod<-SDmodsAIC[order(SDmodsAIC$AICc),]
write.csv(SDmod,"~/Dropbox/Warbler.Molt.Migration/SD_models.csv")



##################################################################
##Seasonal Dichromatism extent models
##################################################################


TDmods<-list()

TDmods[[1]]<-pgls(data=PA.dat,tdExtent~daylength)
TDmods[[2]]<-pgls(data=PA.dat,tdExtent~migDist)
TDmods[[3]]<-pgls(data=PA.dat,tdExtent~breeding_forage_stratum)
TDmods[[4]]<-pgls(data=PA.dat,tdExtent~winter_strat)
TDmods[[5]]<-pgls(data=PA.dat,tdExtent~absBreedLat)
TDmods[[6]]<-pgls(data=PA.dat,tdExtent~breedSolar)
TDmods[[7]]<-pgls(data=PA.dat,tdExtent~SunExpos)
TDmods[[8]]<-pgls(data=PA.dat,tdExtent~minTemp_Breed)
TDmods[[9]]<-pgls(data=PA.dat,tdExtent~maxTemp_Breed)
TDmods[[10]]<-pgls(data=PA.dat,tdExtent~avgTemp_Breed)
TDmods[[11]]<-pgls(data=PA.dat,tdExtent~elevation_Winter)
TDmods[[12]]<-pgls(data=PA.dat,tdExtent~elevation_Winter)
TDmods[[13]]<-pgls(data=PA.dat,tdExtent~paextent)
TDmods[[14]]<-pgls(data=PA.dat,tdExtent~paextent*winter_strat)
TDmods[[15]]<-pgls(data=PA.dat,tdExtent~paextent*breeding_forage_stratum)
TDmods[[16]]<-pgls(data=PA.dat,tdExtent~paextent*winter_strat*breeding_forage_stratum)
TDmods[[17]]<-pgls(data=PA.dat,tdExtent~paextent*winter_strat*daylength)
TDmods[[18]]<-pgls(data=PA.dat,tdExtent~winter_strat*breeding_forage_stratum)
TDmods[[19]]<-pgls(data=PA.dat,tdExtent~winter_strat*avgTemp_Breed)
TDmods[[20]]<-pgls(data=PA.dat,tdExtent~migDist*winter_strat)
TDmods[[21]]<-pgls(data=PA.dat,tdExtent~migDist*breeding_forage_stratum)
TDmods[[22]]<-pgls(data=PA.dat,tdExtent~migDist*breeding_forage_stratum*winter_strat)






TDmodsAIC<-matrix(nrow=22,ncol=4)
for (i in 1:22){
	TDmodsAIC[i,2]<-TDmods[[i]]$aicc
	TDmodsAIC[i,3]<-summary(TDmods[[i]])$adj.r.squared
	TDmodsAIC[i,4]<-pf(summary(TDmods[[i]])$fstatistic[1],summary(TDmods[[i]])$fstatistic[2],summary(TDmods[[i]])$fstatistic[3],lower.tail=F)
	TDmodsAIC[i,1]<-paste("TD_Model",i,sep="_")
		}
		
TDmodsAIC<-data.frame(TDmodsAIC)
TDmodsAIC[,2]<-as.numeric(as.character(TDmodsAIC[,2]))
TDmodsAIC[,3]<-as.numeric(as.character(TDmodsAIC[,3]))


names(TDmodsAIC)<-c("model","AICc","adj_r_sq","p-value")

TDmod<-TDmodsAIC[order(TDmodsAIC$AICc),]
write.csv(TDmod,"~/Dropbox/Warbler.Molt.Migration/TD_models.csv")


pdf("~/Dropbox/Warbler.Molt.Migration/TD_extent~winter_strat.pdf")


ggplot(PA.dat$data,aes(x=tdExtent,y=paextent*winter_strat*breeding_forage_stratum))+geom_jitter(width=.1,height=.1,size=4,alpha=.5)+theme_bw()+geom_abline()

dev.off()




d<-data.frame(PA.dat$data$tdExtent,PA.dat$data$paextent,PA.dat$data$winter_strat,PA.dat$data$breeding_forage_stratum,PA.dat$data$migDist)
names(d)<-c("tdextent","paextent","winter_strat","breeding_forage_stratum","Migratory_distance")
td16<-melt(d,id="tdextent")

pdf("~/Dropbox/Warbler.Molt.Migration/tdextenttopVars.pdf")
ggplot(td16,aes(y=tdextent,x=value))+geom_jitter()+theme_bw()+geom_smooth(method="lm",alpha=.2,colour="gray")+facet_wrap(~variable,scales="free")
dev.off()



#### Great, this works and makes sense too. Follow up with some mixed models on feather regions or mod units###




###in path analysis, bring together all the standardized regression coefficients for figure

