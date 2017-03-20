#phylogenetic correlation matrix

library(caper)
library(phytools)

#####Save all correlation and p-value matrices to external files here

###this will be a function that constructs a phylogenetic correlation matrix
### We need PGLS so that we can run it under Brownian and OU models


####This makes a correlation matrix using PGLS - convert it to the right class so its plottable and rerun all the correlation matrix work
#### Can be altered to fit OU and BM models
pgls.matrix<-function(matrix,phy,names){
	comdat<-comparative.data(phy,matrix,names,vcv=TRUE,vcv.dim=3,warn.dropped=FALSE)
	r<-nrow(matrix)
	c<-ncol(comdat$data)
	out<-list()
	out$r2.mat<-out$p.mat<-matrix(NA,nrow=c,ncol=c)
	rownames(out$r2.mat)<-colnames(out$r2.mat)<-rownames(out$p.mat)<-colnames(out$p.mat)<-colnames(comdat$data)
			for(j in 1:c){
				for(k in 1:c){
			tmp<-pgls(comdat$data[,j]~comdat$data[,k],comdat)
			out$r2.mat[j,k]<-summary(tmp)$adj.r.squared
			out$p.mat[j,k]<-pf(summary(tmp)$fstatistic[1],summary(tmp)$fstatistic[2],summary(tmp)$fstatistic[3],lower.tail=F)
			cat("\r","row",j,"column",k)
			}
		}
		return(out)
		
	}
	

	
	
pgls.mat<-pgls.matrix(Data,warTree,Data$names)	

library(ape)
library(geiger)
library(nlme)
library(phytools)
library(corrplot)
library(wesanderson)
library(reshape)
library(Hmisc)

pal<-wes_palette("Zissou",100,type="continuous")


Data<-read.csv("~/Dropbox/Warbler.Molt.Migration/Data_treeNames.csv",row.names=1)

Data$names<-rownames(Data)


warTree<-read.nexus("/Users/ryanterrill/Dropbox/Warbler.Molt.Migration/Parulidae_phylogeny/Trees/r8stree.tre")
warTree<-multi2di(warTree)





pdf("~/Dropbox/Warbler.Molt.Migration/PGLSmultiCorr_signif.pdf")
corrplot(pgls.mat$r2.mat,type="upper",tl.pos="tp",tl.col="black",method="square",p.mat=pgls.mat$p.mat, sig.level=0.05,insig="blank",col=pal)
corrplot(pgls.mat$r2.mat,add=TRUE, type="lower", method="number",
diag=FALSE,tl.pos="n", cl.pos="n",tl.col="black",tl.cex=.3,number.cex=.35,p.mat=pgls.mat$p.mat, sig.level=0.05,insig="blank",col=pal)

dev.off()

pdf("~/Dropbox/Warbler.Molt.Migration/PGLSmultiCorr_all.pdf")
corrplot(pgls.mat$r2.mat,type="upper",tl.pos="tp",tl.col="black",method="square",col=pal)
corrplot(pgls.mat$r2.mat,add=TRUE, type="lower", method="number",
diag=FALSE,tl.pos="n", cl.pos="n",tl.col="black",tl.cex=.3,number.cex=.35,col=pal)

dev.off()


#### non-phylo controlled ###
resVars<-c("paextent","pfextent","tdExtent","sdExtent")

res<-Data[,resVars]
pred<-Data[,!(names(Data)%in%resVars)]

plotData<-Data
plotData$species<-rownames(Data)

paMelt<-melt(plotData,id=c("species","paextent"))


pdf("~/Dropbox/Warbler.Molt.Migration/PAcorrs.pdf")
ggplot(paMelt,aes(value,paextent))+geom_point(col="black",cex=.4)+facet_wrap(~variable,scales="free")+theme_void()+theme(strip.text.x = element_text(size = 7),axis.text.x = element_text(size=5),strip.background=element_rect(fill="white"))+geom_smooth(method="lm",color="gray",lwd=.5)+labs(title="Correlates of extent of prealternate molt")+theme(plot.title = element_text(hjust = 0.5))+theme(axis.line=element_line())
dev.off()


pfMelt<-melt(plotData,id=c("species","pfextent"))


pdf("~/Dropbox/Warbler.Molt.Migration/PFcorrs.pdf")
ggplot(pfMelt,aes(value,pfextent))+geom_point(col="black",cex=.4)+facet_wrap(~variable,scales="free")+theme_void()+theme(strip.text.x = element_text(size = 7),axis.text.x = element_text(size=5),strip.background=element_rect(fill="white"))+geom_smooth(method="lm",color="gray",lwd=.5)+labs(title="Correlates of extent of preformative molt")+theme(plot.title = element_text(hjust = 0.5))+theme(axis.line=element_line())
dev.off()


tdMelt<-melt(plotData,id=c("species","tdExtent"))


pdf("~/Dropbox/Warbler.Molt.Migration/tdcorrs.pdf")
ggplot(tdMelt,aes(value,tdExtent))+geom_point(col="black",cex=.4)+facet_wrap(~variable,scales="free")+theme_void()+theme(strip.text.x = element_text(size = 7),axis.text.x = element_text(size=5),strip.background=element_rect(fill="white"))+geom_smooth(method="lm",color="gray",lwd=.5)+labs(title="Correlates of extent of seasonal dimorphism")+theme(plot.title = element_text(hjust = 0.5))+theme(axis.line=element_line())
dev.off()



sdMelt<-melt(plotData,id=c("species","sdExtent"))


pdf("~/Dropbox/Warbler.Molt.Migration/sdcorrs.pdf")
ggplot(sdMelt,aes(value,sdExtent))+geom_point(col="black",cex=.4)+facet_wrap(~variable,scales="free")+theme_void()+theme(strip.text.x = element_text(size = 7),axis.text.x = element_text(size=5),strip.background=element_rect(fill="white"))+geom_smooth(method="lm",color="gray",lwd=.5)+labs(title="Correlates of extent of sexual dimorphism")+theme(plot.title = element_text(hjust = 0.5))+theme(axis.line=element_line())
dev.off()


pdf("~/Dropbox/Warbler.Molt.Migration/picDataPairsplot.pdf")
pairs(picData,cex=.1,cex.labels=.1)
dev.off()



##############################################################################################################

#plots above are not phylogenetically independent - plots below are phylogenetic independent contrasts


###Make sure everything is PGLS instead of PIC 

##############################################################################################################



##pgls.mat$r2.mat
##pgls.mat$p.mat





################################################################################################################
###Correlations of individual body pearts within traits, between response variables, and with all variables
###############################################################################################################

PA_tr<-read.csv("/Users/ryanterrill/Dropbox/Warbler.Molt.Migration/Molt\ extent/Parulidae_PA.csv",row.names=1)
PF_tr<-read.csv("/Users/ryanterrill/Dropbox/Warbler.Molt.Migration/Molt\ extent/Parulidae_PF.csv",row.names=1)
SD_tr<-read.csv("~/Dropbox/Warbler.Molt.Migration/Parulidae_sexual_dimorphism.csv",row.names=1)
TD_tr<-read.csv("~/Dropbox/Warbler.Molt.Migration/Parulidae_temporal_dimorphism.csv",row.names=1)




taxa.key=read.delim('~/Dropbox/Warbler.Molt.Migration/Parulidae_taxonomies.txt',stringsAsFactors=F)


######################
###Need to line these up with Tree names###

zeroCols<-c("PA_primaries","PA_primary_coverts","PA_alula","treeName")
keeps<-setdiff(colnames(PA_tr),zeroCols)
PA_tr_var<-PA_tr[keeps]
PA_tr_var$names<-rownames(PA_tr_var)
PA_tr_c<-PA_tr_var[complete.cases(PA_tr_var),]


pglsPA<-pgls.matrix(PA_tr_c,warTree,PA_tr_c$names)



pdf("~/Dropbox/Warbler.Molt.Migration/PAphyloModularity.pdf")

corrplot(pglsPA$r2.mat,order="hclust",col=pal,method="square",p.mat=pglsPA$p.mat,insig="blank",sig.level=.05,title="Prealternate molt",mar=c(0,0,1,0),tl.col="black")
dev.off()	
#####################


######################
colSums(PF_tr,na.rm=TRUE)
zeroPF<-c("treeName")
keeps<-setdiff(colnames(PF_tr),zeroCols)
PF_tr_var<-PF_tr[keeps]
PF_tr_var$names<-rownames(PF_tr_var)

PF_tr_c<-PF_tr_var[complete.cases(PF_tr_var),]



pglsPF<-pgls.matrix(PF_tr_c,warTree,PF_tr_c$names)



pdf("~/Dropbox/Warbler.Molt.Migration/PFphyloModularity.pdf")

corrplot(pglsPF$r2.mat,order="hclust",col=pal,method="square",p.mat=pglsPF$p.mat,insig="blank",sig.level=.05,title="Preformative Molt",mar=c(0,0,1,0),tl.col="black")
dev.off()	
#####################
######################

colSums(SD_tr,na.rm=TRUE)

SD_tr$names<-rownames(SD_tr)

pglsSD<-pgls.matrix(SD_tr,warTree,SDData[complete.cases(SDData),]$names)

pdf("~/Dropbox/Warbler.Molt.Migration/SDphyloModularity.pdf")

corrplot(pglsSD$r2.mat,order="hclust",col=pal,method="square",p.mat=pglsSD$p.mat,insig="blank",sig.level=.05,title="Sexual Dichromatism",mar=c(0,0,1,0),tl.col="black")
dev.off()	
#####################
######################

colSums(TD_tr,na.rm=TRUE)


zeroCols<-c("TD_primary_coverts","TD_alula","TD_rectrices","TD_primaries","TD_secondaries")
keeps<-setdiff(colnames(TD_tr),zeroCols)
TD_tr_var<-TD_tr[keeps]
TD_tr_var$names<-rownames(TD_tr_var)


pglsTD<-pgls.matrix(TD_tr_var,warTree,TD_tr_var$names)

pdf("~/Dropbox/Warbler.Molt.Migration/TDphyloModularity.pdf")

corrplot(pglsTD$r2.mat,order="hclust",col=pal,method="square",p.mat=pglsTD$p.mat,insig="blank",sig.level=.05,title="Seasonal Dichromatism",mar=c(0,0,1,0),tl.col="black")
dev.off()
#####################

########################################################################################################################################################################

###Now phylo correlation between each tract and predictor variables




###########################################


###rewritten correlation matrix code

unev.cor.mtest <- function(mat, conf.level = 0.95){
  mat <- as.matrix(mat)
  	r<-nrow(mat)
    n <- ncol(mat)
    p.mat <- lowCI.mat <- uppCI.mat <- matrix(NA, r, n)
    for(i in 1:r){
        for(j in 1:n){
            tmp <- cor.test(mat[,i], mat[,j], conf.level = conf.level)
            p.mat[i,j] <- tmp$p.value
        }
    }
    return(list(p.mat))
}


############################################



paDat<-merge(PA_tr,Data,by=0)

rownames(paDat)<-paDat$Row.names

paDat<-paDat[,-1]

zeroCols<-c("PA_primaries","PA_primary_coverts","PA_alula","treeName")
keeps<-setdiff(colnames(paDat),zeroCols)
paDat<-paDat[keeps]

paDat<-paDat[complete.cases(paDat),]

PAcorr<-pgls.matrix(paDat,warTree,paDat$names)



pdf("~/Dropbox/Warbler.Molt.Migration/PA_tract_predictors.pdf")

corrplot(PAcorr$r2.mat,col=pal,method="square",p.mat=PAcorr$p.mat,insig="blank",sig.level=.05,title="Prealternate Molt PGLS predictors",mar=c(0.2,0.2,1,0.2), tl.col="black",tl.cex=.5)

dev.off()



############################################

############################################



PFDat<-merge(PF_tr,Data,by=0)

rownames(PFDat)<-PFDat$Row.names

PFDat<-PFDat[,-1]

zeroCols<-c("PF_primaries","PF_primary_coverts","PF_alula","treeName")
keeps<-setdiff(colnames(PFDat),zeroCols)
PFDat<-PFDat[keeps]

PFDat<-PFDat[complete.cases(PFDat),]

PFcorr<-pgls.matrix(PFDat,warTree,PFDat$names)



pdf("~/Dropbox/Warbler.Molt.Migration/PF_tract_predictors.pdf")

corrplot(PFcorr$r2.mat,col=pal,method="square",p.mat=PFcorr$p.mat,insig="blank",sig.level=.05,title="Preformative Molt PGLS predictors",mar=c(0.2,0.2,1,0.2), tl.col="black",tl.cex=.5)

dev.off()



############################################	


############################################



SDDat<-merge(SD_tr,Data,by=0)

rownames(SDDat)<-SDDat$Row.names

SDDat<-SDDat[,-1]

zeroCols<-c("SD_primaries","SD_primary_coverts","SD_alula","treeName","names.x","names.y")
keeps<-setdiff(colnames(SDDat),zeroCols)

SDDat<-SDDat[keeps]
SDDat$names<-rownames(SDDat)

SDDat<-SDDat[complete.cases(SDDat),]



SDcorr<-pgls.matrix(SDDat,warTree,SDDat$names)



pdf("~/Dropbox/Warbler.Molt.Migration/SD_tract_predictors.pdf")

corrplot(SDcorr$r2.mat,col=pal,method="square",p.mat=SDcorr$p.mat,insig="blank",sig.level=.05,title="Sexual Dichromatism PGLS predictors",mar=c(0.2,0.2,1,0.2), tl.col="black",tl.cex=.5)

dev.off()



############################################


############################################



TDDat<-merge(TD_tr,Data,by=0)

rownames(TDDat)<-TDDat$Row.names

TDDat<-TDDat[,-1]

zeroCols<-c("TD_primaries","TD_primary_coverts","TD_rectrices","TD_secondaries","TD_alula")
keeps<-setdiff(colnames(TDDat),zeroCols)
TDDat<-TDDat[keeps]

TDDat<-TDDat[complete.cases(TDDat),]

TDcorr<-pgls.matrix(TDDat,warTree,TDDat$names)



pdf("~/Dropbox/Warbler.Molt.Migration/TD_tract_predictors.pdf")

corrplot(TDcorr$r2.mat,col=pal,method="square",p.mat=TDcorr$p.mat,insig="blank",sig.level=.05,title="Seasonal Dichromatism PGLS predictors",mar=c(0.2,0.2,1,0.2), tl.col="black",tl.cex=.5)

dev.off()



############################################

pdf("~/Dropbox/Warbler.Molt.Migration/TD_extent~winter_strat.pdf")


ggplot(TDDat,aes(x=tdExtent,y=winter_strat))+geom_jitter(width=.1,height=.1,size=4,alpha=.5)+theme_bw()

dev.off()


pdf("~/Dropbox/Warbler.Molt.Migration/TD_extent~breeding_strat.pdf")

ggplot(TDDat,aes(x=tdExtent,y=breeding_forage_stratum))+geom_jitter(width=.1,height=.1,size=4,alpha=.5)+theme_bw()

dev.off()
