####

Mixed models

###PGLS, AIC weights for all likely models for the four predictor variables
###Mixed Models built from examining PGLS response of individual variables

library(caper)
library(MuMIn)
library(MCMCglmm)

##Maybe just use PGLS and design models by hand
## Use correlation matrix to build models

##################################################################
###PA extent models
##################################################################


source("~/Dropbox/Warbler.Molt.Migration/model_selection.R")

pat<-treedata(warTree,na.omit((PA[,1,drop=FALSE])))
PAtree<-pat$phy


Data<-read.csv("~/Dropbox/Warbler.Molt.Migration/Data_treeNames.csv",row.names=1)
Data$names<-rownames(Data)

###BM is the supported model for PA
PA.dat<-comparative.data(warTree,Data,names,vcv=TRUE,vcv.dim=3,warn.dropped=TRUE)


## PA extent predictor mixed models
PAmods<-list()
PAmods[[1]]<-pgls(data=PA.dat,paextent~daylength,lambda="ML")	
PAmods[[2]]<-pgls(data=PA.dat,paextent~migDist,lambda="ML")
PAmods[[3]]<-pgls(data=PA.dat,paextent~breeding_forage_stratum,lambda="ML")
PAmods[[4]]<-pgls(data=PA.dat,paextent~winter_strat,lambda="ML")
PAmods[[5]]<-pgls(data=PA.dat,paextent~breedSolar,lambda="ML")
PAmods[[6]]<-pgls(data=PA.dat,paextent~precip_Breed,lambda="ML")
PAmods[[7]]<-pgls(data=PA.dat,paextent~minTemp_Breed,lambda="ML")
PAmods[[8]]<-pgls(data=PA.dat,paextent~avgTemp_Breed,lambda="ML")
PAmods[[9]]<-pgls(data=PA.dat,paextent~Mass,lambda="ML")
PAmods[[10]]<-pgls(data=PA.dat,paextent~daylength*migDist,lambda="ML")
PAmods[[11]]<-pgls(data=PA.dat,paextent~daylength*breedSolar,lambda="ML")
PAmods[[12]]<-pgls(data=PA.dat,paextent~daylength*breeding_forage_stratum,lambda="ML")
PAmods[[13]]<-pgls(data=PA.dat,paextent~daylength*winter_strat,lambda="ML")	
PAmods[[14]]<-pgls(data=PA.dat,paextent~daylength*precip_Breed,lambda="ML")	
PAmods[[15]]<-pgls(data=PA.dat,paextent~daylength*minTemp_Breed,lambda="ML")	
PAmods[[16]]<-pgls(data=PA.dat,paextent~daylength*breeding_forage_stratum*winter_strat,lambda="ML")	
PAmods[[17]]<-pgls(data=PA.dat,paextent~daylength*breeding_forage_stratum*winter_strat*breedSolar,lambda="ML")	
PAmods[[18]]<-pgls(data=PA.dat,paextent~breeding_forage_stratum*winter_strat*breedSolar,lambda="ML")	
PAmods[[19]]<-pgls(data=PA.dat,paextent~daylength*breeding_forage_stratum*winter_strat*Mass,lambda="ML")	
PAmods[[20]]<-pgls(data=PA.dat,paextent~daylength*breedSolar*sol,lambda="ML")
PAmods[[21]]<-pgls(data=PA.dat,paextent~daylength*sol,lambda="ML")
PAmods[[22]]<-pgls(data=PA.dat,paextent~daylength*winSolar,lambda="ML")
PAmods[[23]]<-pgls(data=PA.dat,paextent~migDist*breedSolar,lambda="ML")
PAmods[[24]]<-pgls(data=PA.dat,paextent~migDist*sol,lambda="ML")
PAmods[[25]]<-pgls(data=PA.dat,paextent~migDist*winSolar,lambda="ML")
PAmods[[26]]<-pgls(data=PA.dat,paextent~migDist*daylength*breedSolar,lambda="ML")

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
PFmods[[1]]<-pgls(data=PA.dat,pfextent~sdExtent,lambda="ML")
PFmods[[2]]<-pgls(data=PA.dat,pfextent~daylength,lambda="ML")
PFmods[[3]]<-pgls(data=PA.dat,pfextent~absWinLat,lambda="ML")
PFmods[[4]]<-pgls(data=PA.dat,pfextent~breeding_forage_stratum,lambda="ML")
PFmods[[5]]<-pgls(data=PA.dat,pfextent~nest_site_stratum,lambda="ML")
PFmods[[6]]<-pgls(data=PA.dat,pfextent~winter_strat,lambda="ML")
PFmods[[7]]<-pgls(data=PA.dat,pfextent~breedSolar,lambda="ML")
PFmods[[8]]<-pgls(data=PA.dat,pfextent~SunExpos,lambda="ML")
PFmods[[9]]<-pgls(data=PA.dat,pfextent~precip_Breed,lambda="ML")
PFmods[[10]]<-pgls(data=PA.dat,pfextent~minTemp_Breed,lambda="ML")
PFmods[[11]]<-pgls(data=PA.dat,pfextent~avgTemp_Breed,lambda="ML")
PFmods[[12]]<-pgls(data=PA.dat,pfextent~maxTemp_Breed,lambda="ML")
PFmods[[13]]<-pgls(data=PA.dat,pfextent~sdExtent*daylength,lambda="ML")
PFmods[[14]]<-pgls(data=PA.dat,pfextent~daylength*absWinLat,lambda="ML")
PFmods[[15]]<-pgls(data=PA.dat,pfextent~daylength*breeding_forage_stratum,lambda="ML")
PFmods[[16]]<-pgls(data=PA.dat,pfextent~daylength*nest_site_stratum,lambda="ML")
PFmods[[17]]<-pgls(data=PA.dat,pfextent~daylength*breedSolar,lambda="ML")
PFmods[[18]]<-pgls(data=PA.dat,pfextent~daylength*sol,lambda="ML")
PFmods[[19]]<-pgls(data=PA.dat,pfextent~breeding_forage_stratum*breedSolar,lambda="ML")
PFmods[[20]]<-pgls(data=PA.dat,pfextent~daylength*nest_site_stratum,lambda="ML")


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

SDmods[[1]]<-pgls(data=PA.dat,sdExtent~daylength,lambda="ML")
SDmods[[2]]<-pgls(data=PA.dat,sdExtent~migDist,lambda="ML")
SDmods[[3]]<-pgls(data=PA.dat,sdExtent~absBreedLat,lambda="ML")
SDmods[[4]]<-pgls(data=PA.dat,sdExtent~absWinLat,lambda="ML")
SDmods[[5]]<-pgls(data=PA.dat,sdExtent~breeding_forage_stratum,lambda="ML")
SDmods[[6]]<-pgls(data=PA.dat,sdExtent~nest_site_stratum,lambda="ML")
SDmods[[7]]<-pgls(data=PA.dat,sdExtent~winter_hab,lambda="ML")
SDmods[[8]]<-pgls(data=PA.dat,sdExtent~sol,lambda="ML")
SDmods[[9]]<-pgls(data=PA.dat,sdExtent~breedSolar,lambda="ML")
SDmods[[10]]<-pgls(data=PA.dat,sdExtent~SunExpos,lambda="ML")
SDmods[[11]]<-pgls(data=PA.dat,sdExtent~minTemp_Breed,lambda="ML")
SDmods[[12]]<-pgls(data=PA.dat,sdExtent~maxTemp_Breed,lambda="ML")
SDmods[[13]]<-pgls(data=PA.dat,sdExtent~avgTemp_Breed,lambda="ML")
SDmods[[14]]<-pgls(data=PA.dat,sdExtent~Mass,lambda="ML")
SDmods[[15]]<-pgls(data=PA.dat,sdExtent~daylength*breeding_forage_stratum,lambda="ML")
SDmods[[16]]<-pgls(data=PA.dat,sdExtent~daylength*nest_site_stratum,lambda="ML")
SDmods[[17]]<-pgls(data=PA.dat,sdExtent~daylength*breeding_forage_stratum,lambda="ML")
SDmods[[18]]<-pgls(data=PA.dat,sdExtent~absBreedLat*breeding_forage_stratum,lambda="ML")
SDmods[[19]]<-pgls(data=PA.dat,sdExtent~absBreedLat*nest_site_stratum,lambda="ML")
SDmods[[20]]<-pgls(data=PA.dat,sdExtent~winter_hab*nest_site_stratum,lambda="ML")
SDmods[[21]]<-pgls(data=PA.dat,sdExtent~absBreedLat*Mass,lambda="ML")
SDmods[[22]]<-pgls(data=PA.dat,sdExtent~absWinLat*breeding_forage_stratum,lambda="ML")


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

TDmods[[1]]<-pgls(data=PA.dat,tdExtent~daylength,lambda="ML")
TDmods[[2]]<-pgls(data=PA.dat,tdExtent~migDist,lambda="ML")
TDmods[[3]]<-pgls(data=PA.dat,tdExtent~breeding_forage_stratum,lambda="ML")
TDmods[[4]]<-pgls(data=PA.dat,tdExtent~winter_strat,lambda="ML")
TDmods[[5]]<-pgls(data=PA.dat,tdExtent~absBreedLat,lambda="ML")
TDmods[[6]]<-pgls(data=PA.dat,tdExtent~breedSolar,lambda="ML")
TDmods[[7]]<-pgls(data=PA.dat,tdExtent~SunExpos,lambda="ML")
TDmods[[8]]<-pgls(data=PA.dat,tdExtent~minTemp_Breed,lambda="ML")
TDmods[[9]]<-pgls(data=PA.dat,tdExtent~maxTemp_Breed,lambda="ML")
TDmods[[10]]<-pgls(data=PA.dat,tdExtent~avgTemp_Breed,lambda="ML")
TDmods[[11]]<-pgls(data=PA.dat,tdExtent~elevation_Winter,lambda="ML")
TDmods[[12]]<-pgls(data=PA.dat,tdExtent~elevation_Winter,lambda="ML")
TDmods[[13]]<-pgls(data=PA.dat,tdExtent~paextent,lambda="ML")
TDmods[[14]]<-pgls(data=PA.dat,tdExtent~paextent*winter_strat,lambda="ML")
TDmods[[15]]<-pgls(data=PA.dat,tdExtent~paextent*breeding_forage_stratum,lambda="ML")
TDmods[[16]]<-pgls(data=PA.dat,tdExtent~paextent*winter_strat*breeding_forage_stratum,lambda="ML")
TDmods[[17]]<-pgls(data=PA.dat,tdExtent~paextent*winter_strat*daylength,lambda="ML")
TDmods[[18]]<-pgls(data=PA.dat,tdExtent~winter_strat*breeding_forage_stratum,lambda="ML")
TDmods[[19]]<-pgls(data=PA.dat,tdExtent~winter_strat*avgTemp_Breed,lambda="ML")
TDmods[[20]]<-pgls(data=PA.dat,tdExtent~migDist*winter_strat,lambda="ML")
TDmods[[21]]<-pgls(data=PA.dat,tdExtent~migDist*breeding_forage_stratum,lambda="ML")
TDmods[[22]]<-pgls(data=PA.dat,tdExtent~migDist*breeding_forage_stratum*winter_strat,lambda="ML")






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
