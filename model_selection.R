#Model selection for molts and dichromatism
###Write this to work with the warbler tree - but see about running it over the PP trees



library(gieger)
library(reshape)
library(ggplot2)

warTree<-read.nexus("/Users/ryanterrill/Dropbox/Warbler.Molt.Migration/Parulidae_phylogeny/Trees/r8stree.tre")
warTree<-multi2di(warTree)
warTree<-chronos(warTree)
class(warTree)<-"phylo"

taxa.key=read.delim('~/Dropbox/Warbler.Molt.Migration/Parulidae_taxonomies.txt',stringsAsFactors=F)


PA_mods<-read.csv("~/Dropbox/Warbler.Molt.Migration/PA_modVals.csv",row.names=1)

PA_mods$Birdlife<-rownames(PA_mods)
setdiff(taxa.key$Birdlife,rownames(PA_mods))
p1<-merge(PA_mods,taxa.key,by="Birdlife")
rownames(p1)<-p1$AOU
PA_mods<-p1[!names(p1)%in%c("AOU","Old","Birdlife")]

PF_mods<-read.csv("~/Dropbox/Warbler.Molt.Migration/PF_modVals.csv",row.names=1)

PF_mods$Birdlife<-rownames(PF_mods)
setdiff(taxa.key$Birdlife,rownames(PF_mods))
p1<-merge(PF_mods,taxa.key,by="Birdlife")
rownames(p1)<-p1$AOU
PF_mods<-p1[!names(p1)%in%c("AOU","Old","Birdlife")]

SD_mods<-read.csv("~/Dropbox/Warbler.Molt.Migration/SD_modVals.csv",row.names=1)

SD_mods$Birdlife<-rownames(SD_mods)
setdiff(taxa.key$Birdlife,rownames(SD_mods))
p1<-merge(SD_mods,taxa.key,by="Birdlife")
rownames(p1)<-p1$AOU
SD_mods<-p1[!names(p1)%in%c("AOU","Old","Birdlife")]


TD_mods<-read.csv("~/Dropbox/Warbler.Molt.Migration/TD_modVals.csv",row.names=1)

TD_mods$Birdlife<-rownames(TD_mods)
setdiff(taxa.key$Birdlife,rownames(TD_mods))
p1<-merge(TD_mods,taxa.key,by="Birdlife")
rownames(p1)<-p1$AOU
TD_mods<-p1[!names(p1)%in%c("AOU","Old","Birdlife")]


PA_tr<-read.csv("/Users/ryanterrill/Dropbox/Warbler.Molt.Migration/Molt\ extent/Parulidae_PA.csv",row.names=1)

PA_tr$Birdlife<-rownames(PA_tr)
setdiff(taxa.key$Birdlife,rownames(PA_tr))
p1<-merge(PA_tr,taxa.key,by="Birdlife")
rownames(p1)<-p1$AOU
PA_tr<-p1[!names(p1)%in%c("AOU","Old","Birdlife")]


PF_tr<-read.csv("/Users/ryanterrill/Dropbox/Warbler.Molt.Migration/Molt\ extent/Parulidae_PF.csv",row.names=1)

PF_tr$Birdlife<-rownames(PF_tr)
setdiff(taxa.key$Birdlife,rownames(PF_tr))
p1<-merge(PF_tr,taxa.key,by="Birdlife")
rownames(p1)<-p1$AOU
PF_tr<-p1[!names(p1)%in%c("AOU","Old","Birdlife")]


SD_tr<-read.csv("~/Dropbox/Warbler.Molt.Migration/Parulidae_sexual_dimorphism.csv",row.names=1)

SD_tr$Birdlife<-rownames(SD_tr)
setdiff(taxa.key$Birdlife,rownames(SD_tr))
p1<-merge(SD_tr,taxa.key,by="Birdlife")
rownames(p1)<-p1$AOU
SD_tr<-p1[!names(p1)%in%c("AOU","Old","Birdlife")]

TD_tr<-read.csv("~/Dropbox/Warbler.Molt.Migration/Parulidae_temporal_dimorphism.csv",row.names=1)

TD_tr$Birdlife<-rownames(TD_tr)
setdiff(taxa.key$Birdlife,rownames(TD_tr))
p1<-merge(TD_tr,taxa.key,by="Birdlife")
rownames(p1)<-p1$AOU
TD_tr<-p1[!names(p1)%in%c("AOU","Old","Birdlife")]



PA_ex<-rowSums(PA_tr);write.csv(PA_ex,"~/Dropbox/Warbler.Molt.Migration/PA_extent.csv")
PF_ex<-rowSums(PF_tr);write.csv(PF_ex,"~/Dropbox/Warbler.Molt.Migration/PF_extent.csv")
SD_ex<-rowSums(SD_tr);write.csv(SD_ex,"~/Dropbox/Warbler.Molt.Migration/SD_extent.csv")
TD_ex<-rowSums(TD_tr);write.csv(TD_ex,"~/Dropbox/Warbler.Molt.Migration/TD_extent.csv")



PA_pres<-PA_ex;PA_pres[PA_pres>0]<-1;write.csv(PA_pres,"~/Dropbox/Warbler.Molt.Migration/PA_presence.csv")
PF_pres<-PF_ex;PF_pres[PF_pres>0]<-1;write.csv(PF_pres,"~/Dropbox/Warbler.Molt.Migration/PF_presence.csv")
SD_pres<-SD_ex;SD_pres[SD_pres>0]<-1;write.csv(SD_pres,"~/Dropbox/Warbler.Molt.Migration/SD_presence.csv")
TD_pres<-TD_ex;TD_pres[TD_pres>0]<-1;write.csv(TD_pres,"~/Dropbox/Warbler.Molt.Migration/TD_presence.csv")


colnames(PA_tr)<-sub("PA_","",names(PA_tr))
colnames(PF_tr)<-sub("PF_","",names(PF_tr))
colnames(SD_tr)<-sub("SD_","",names(SD_tr))
colnames(TD_tr)<-sub("TD_","",names(TD_tr))

PA<-data.frame(PA_pres[complete.cases(PA_pres)],PA_ex[complete.cases(PA_ex)],PA_mods,PA_tr[complete.cases(PA_tr),])
colnames(PA)<-c("presence","extent",colnames(PA_mods),colnames(PA_tr))

PF<-data.frame(PF_pres[complete.cases(PF_pres)],PF_ex[complete.cases(PF_ex)],PF_mods,PF_tr[complete.cases(PF_tr),])
colnames(PF)<-c("presence","extent",colnames(PF_mods),colnames(PF_tr))

SD<-data.frame(SD_pres[complete.cases(SD_pres)],SD_ex[complete.cases(SD_ex)],SD_mods,SD_tr[complete.cases(SD_tr),])
colnames(SD)<-c("presence","extent",colnames(SD_mods),colnames(SD_tr))


TD<-data.frame(TD_pres[complete.cases(TD_pres)],TD_ex[complete.cases(TD_ex)],TD_mods,TD_tr[complete.cases(TD_tr),])
colnames(TD)<-c("presence","extent",colnames(TD_mods),colnames(TD_tr))




##AICs for models of evolution


modelAICc<-function(matrix,phy){
	tree.dat<-list()
	AICc<-matrix(nrow=3,ncol=ncol(matrix))
	rownames(AICc)<-c("BM","OU","EB")
	colnames(AICc)<-colnames(matrix)
	for(i in 1:ncol(matrix)){
		tree.dat[[i]]<-treedata(phy,na.omit(matrix[,i,drop=FALSE]))
			AICc[1,i]<-fitContinuous(tree.dat[[i]]$phy,tree.dat[[i]]$dat,model="BM")$opt$aicc
			AICc[2,i]<-fitContinuous(tree.dat[[i]]$phy,tree.dat[[i]]$dat,model="OU")$opt$aicc
			AICc[3,i]<-fitContinuous(tree.dat[[i]]$phy,tree.dat[[i]]$dat,model="EB")$opt$aicc
			}
			return(AICc)
		}
		
		
AICc_weight_matrix<-function(matrix){
	out<-matrix(nrow=nrow(matrix),ncol=ncol(matrix))
	rownames(out)<-rownames(matrix)
	colnames(out)<-colnames(matrix)
	aiccD<-list()
	aicc<-list()
	aw<-list()
	aiccW<-list()
	for(i in 1: ncol(matrix)){
		aicc[[i]]<-matrix[,i]
		aiccD[[i]]<-aicc[[i]]-min(aicc[[i]])
		aw[[i]]<-exp(-0.5*aiccD[[i]])
		aiccW[[i]]<-aw[[i]]/sum(aw[[i]])
	out[,i]<-aiccW[[i]]
	}
	return(out)
}

####PRESENCE NEEDS TO USE FITDISCRETE


##might as well put all groups together in four matrices		
PA_AICc<-modelAICc(PA,warTree)
wtAIC_PA<-AICc_weight_matrix(PA_AICc)

write.csv(PA_AICc,"~/Dropbox/Warbler.Molt.Migration/PA_AICc.csv")
write.csv(wtAIC_PA,"~/Dropbox/Warbler.Molt.Migration/PA_weights.csv")

PF_AICc<-modelAICc(PF,warTree)
wtAIC_PF<-AICc_weight_matrix(PF_AICc)

write.csv(PF_AICc,"~/Dropbox/Warbler.Molt.Migration/PF_AICc.csv")
write.csv(wtAIC_PF,"~/Dropbox/Warbler.Molt.Migration/PF_weights.csv")

SD_AICc<-modelAICc(SD,warTree)
wtAIC_SD<-AICc_weight_matrix(SD_AICc)

write.csv(SD_AICc,"~/Dropbox/Warbler.Molt.Migration/SD_AICc.csv")
write.csv(wtAIC_SD,"~/Dropbox/Warbler.Molt.Migration/SD_weights.csv")

TD_AICc<-modelAICc(TD,warTree)
wtAIC_TD<-AICc_weight_matrix(TD_AICc)

write.csv(TD_AICc,"~/Dropbox/Warbler.Molt.Migration/TD_AICc.csv")
write.csv(wtAIC_TD,"~/Dropbox/Warbler.Molt.Migration/TD_weights.csv")




### combine these and plot by region and parameter, which model is favored overall? 
### are different models more favored for different parameters?


PA_melt<-melt(wtAIC_PA);PA_melt$par<-"PA"
PF_melt<-melt(wtAIC_PF);PF_melt$par<-"PF"
SD_melt<-melt(wtAIC_SD);SD_melt$par<-"SD"
TD_melt<-melt(wtAIC_TD);TD_melt$par<-"TD"

AIC_wts<-rbind(PA_melt,PF_melt,SD_melt,TD_melt)
colnames(AIC_wts)<-c("Model","variable","value","par")

pdf("~/Dropbox/Warbler.Molt.Migration/presence_extent_AICc_weights.pdf")
ggplot(subset(AIC_wts,variable%in%c("extent")),aes(x=Model,y=value,col=Model))+geom_point(aes(size=value,alpha=value))+facet_grid(par~.)+theme_bw()+theme(panel.grid=element_blank(),axis.text.x = element_text(angle=45, hjust=1))+ggtitle("Extent Model AIC weights")+theme(plot.title = element_text(hjust = 0.5))
dev.off()

lvs<-levels(AIC_wts$variable);lvs<-lvs[lvs!=c("presence","extent")]

tracts<-c("alula","back","belly","breast","gr_cov","head","med_cov","primaries","primary_coverts","rectrices","secondaries","tertials")

s<-subset(AIC_wts,variable%in%lvs)

tr<-subset(AIC_wts,variable%in%tracts)

mods<-c("mod_1","mod_2","mod_3","mod_4","mod_5","mod_6")

md<-subset(AIC_wts,variable%in%mods)

pdf("~/Dropbox/Warbler.Molt.Migration/tracts_AICc_weights.pdf")
ggplot(tr,aes(x=variable,y=value,col=AICc_weight,fill=AICc_weight))+geom_point(aes(size=value,alpha=value))+facet_grid(par~.)+theme_bw()+theme(panel.grid=element_blank(),axis.text.x = element_text(angle=45, hjust=1))+geom_line(aes(group=AICc_weight),alpha=.2)
dev.off()

pdf("~/Dropbox/Warbler.Molt.Migration/mods_AICc_weights.pdf")
ggplot(md,aes(x=variable,y=value,col=AICc_weight,fill=AICc_weight))+geom_point(aes(size=value,alpha=value))+facet_grid(par~.)+theme_bw()+theme(panel.grid=element_blank(),axis.text.x = element_text(angle=45, hjust=1))+geom_line(aes(group=AICc_weight))
dev.off()



pdf("~/Dropbox/Warbler.Molt.Migration/sum_Units_AICc_weights.pdf")
ggplot(s,aes(x=par,y=value,col=AICc_weight,fill=AICc_weight))+geom_boxplot(alpha=.5)+theme_bw()+theme(panel.grid=element_blank(),axis.text.x = element_text(angle=45, hjust=1))+geom_point(aes(fill=AICc_weight),position = position_jitterdodge(),alpha=.8)
dev.off()





#### Calculate rateparameters for BM, OU, EB
#### plot model parameters along with AICc weights for just BM and OU


rateParam<-function(matrix,phy){
	tree.dat<-list()
	sigsq<-matrix(nrow=3,ncol=ncol(matrix))
	rownames(sigsq)<-c("BM","OU","EB")
	colnames(sigsq)<-colnames(matrix)
	for(i in 1:ncol(matrix)){
		tree.dat[[i]]<-treedata(phy,na.omit(matrix[,i,drop=FALSE]))
			sigsq[1,i]<-fitContinuous(tree.dat[[i]]$phy,tree.dat[[i]]$dat,model="BM")$opt$sigsq
			sigsq[2,i]<-fitContinuous(tree.dat[[i]]$phy,tree.dat[[i]]$dat,model="OU")$opt$sigsq
			sigsq[3,i]<-fitContinuous(tree.dat[[i]]$phy,tree.dat[[i]]$dat,model="EB")$opt$sigsq
			}
			return(sigsq)
		}
		
PA_sigsq<-rateParam(PA,warTree)
write.csv(PA_AICc,"~/Dropbox/Warbler.Molt.Migration/PA_rateParams.csv")

PF_sigsq<-rateParam(PF,warTree)
write.csv(PF_AICc,"~/Dropbox/Warbler.Molt.Migration/PF_rateParams.csv")

SD_sigsq<-rateParam(SD,warTree)
write.csv(SD_AICc,"~/Dropbox/Warbler.Molt.Migration/SD_rateParams.csv")

TD_sigsq<-rateParam(TD,warTree)
write.csv(TD_AICc,"~/Dropbox/Warbler.Molt.Migration/TD_rateParams.csv")


####Plotting parameter rates


PA_sigsq_melt<-melt(PA_sigsq);PA_sigsq_melt$par<-"PA"
PF_sigsq_melt<-melt(PF_sigsq);PF_sigsq_melt$par<-"PF"
SD_sigsq_melt<-melt(SD_sigsq);SD_sigsq_melt$par<-"SD"
TD_sigsq_melt<-melt(TD_sigsq);TD_sigsq_melt$par<-"TD"

sigsq_all<-rbind(PA_sigsq_melt,PF_sigsq_melt,SD_sigsq_melt,TD_sigsq_melt)


colnames(sigsq_all)<-c("Model","variable","value","par")

pdf("~/Dropbox/Warbler.Molt.Migration/presence_extent_Model_parameters.pdf")
ggplot(subset(sigsq_all,variable%in%c("extent")),aes(x=Model,y=value,col=Model))+geom_point(cex=6,alpha=.8)+facet_grid(par~.)+theme_bw()+theme(panel.grid=element_blank(),axis.text.x = element_text(angle=45, hjust=1))+ggtitle("Extent Model parameter rates")+theme(plot.title = element_text(hjust = 0.5))
dev.off()

lvs<-levels(sigsq_all$variable);lvs<-lvs[lvs!=c("presence","extent")]

tracts<-c("alula","back","belly","breast","gr_cov","head","med_cov","primaries","primary_coverts","rectrices","secondaries","tertials")

s_sig<-subset(sigsq_all,variable%in%lvs)

tr_sig<-subset(sigsq_all,variable%in%tracts)

mods<-c("mod_1","mod_2","mod_3","mod_4","mod_5","mod_6")

md_sig<-subset(sigsq_all,variable%in%mods)

pdf("~/Dropbox/Warbler.Molt.Migration/tracts_parameter_Values.pdf")
ggplot(tr_sig,aes(x=variable,y=value,col=Model,fill=Model))+geom_point(cex=2)+facet_grid(par~.,scales="free_y")+theme_bw()+theme(panel.grid=element_blank(),axis.text.x = element_text(angle=45, hjust=1))+geom_line(aes(group=Model))
dev.off()

pdf("~/Dropbox/Warbler.Molt.Migration/mods_parameters_Values.pdf")
ggplot(md_sig,aes(x=variable,y=value,col=Model,fill=Model))+geom_point(cex=2)+facet_grid(par~.,scales="free_y")+theme_bw()+theme(panel.grid=element_blank(),axis.text.x = element_text(angle=45, hjust=1))+geom_line(aes(group=Model))
dev.off()



pdf("~/Dropbox/Warbler.Molt.Migration/sum_Units_parameter_values.pdf")
ggplot(s_sig,aes(x=par,y=value,col=Model,fill=Model))+geom_boxplot(alpha=.5,outlier.shape=NA)+theme_bw()+theme(panel.grid=element_blank(),axis.text.x = element_text(angle=45, hjust=1))+scale_y_continuous(limits = c(0,.75))
dev.off()







PA_all<-data.frame(PA_melt,PA_sigsq_melt);PA_all<-PA_all[,c(1:4,7)];colnames(PA_all)<-c("Model","variable","AIC","par","rate")
PF_all<-data.frame(PF_melt,PF_sigsq_melt);PF_all<-PF_all[,c(1:4,7)];colnames(PF_all)<-c("Model","variable","AIC","par","rate")
SD_all<-data.frame(SD_melt,SD_sigsq_melt);SD_all<-SD_all[,c(1:4,7)];colnames(SD_all)<-c("Model","variable","AIC","par","rate")
TD_all<-data.frame(TD_melt,TD_sigsq_melt);TD_all<-TD_all[,c(1:4,7)];colnames(TD_all)<-c("Model","variable","AIC","par","rate")

all_mods<-rbind(PA_all,PF_all,SD_all,TD_all)

tr_all<-subset(all_mods,variable%in%tracts)

pdf("~/Dropbox/Warbler.Molt.Migration/weight_par_Values.pdf")
ggplot(tr_all,aes(x=variable,y=rate,col=Model,fill=Model))+geom_point(aes(size=AIC,alpha=AIC))+facet_grid(par~.,scales="free_y")+theme_bw()+theme(panel.grid=element_blank(),axis.text.x = element_text(angle=45, hjust=1))+geom_line(aes(group=Model),alpha=.2)+ggtitle("AIC weighted parameter rates")+theme(plot.title = element_text(hjust = 0.5))
dev.off()




pdf("~/Dropbox/Warbler.Molt.Migration/presence_extent_Model__weighted_parameters.pdf")
ggplot(subset(all_mods,variable%in%c("extent")),aes(x=Model,y=rate,col=Model))+geom_point(aes(size=AIC,alpha=AIC))+facet_grid(par~.)+theme_bw()+theme(panel.grid=element_blank(),axis.text.x = element_text(angle=45, hjust=1))+ggtitle("Extent Model parameter rates")+theme(plot.title = element_text(hjust = 0.5))
dev.off()

lvs<-levels(sigsq_all$variable);lvs<-lvs[lvs!=c("presence","extent")]

tracts<-c("alula","back","belly","breast","gr_cov","head","med_cov","primaries","primary_coverts","rectrices","secondaries","tertials")

s_all<-subset(sigsq_all,variable%in%lvs)

mods<-c("mod_1","mod_2","mod_3","mod_4","mod_5","mod_6")

md_all<-subset(all_mods,variable%in%mods)



pdf("~/Dropbox/Warbler.Molt.Migration/mods_parameters_Values.pdf")
ggplot(md_all,aes(x=variable,y=rate,col=Model,fill=Model))+geom_point(aes(size=AIC,alpha=AIC))+facet_grid(par~.,scales="free_y")+theme_bw()+theme(panel.grid=element_blank(),axis.text.x = element_text(angle=45, hjust=1))+geom_line(aes(group=Model))+ggtitle("AIC weighted parameter rates")+theme(plot.title = element_text(hjust = 0.5))
dev.off()



pdf("~/Dropbox/Warbler.Molt.Migration/sum_Units_parameter_values.pdf")
ggplot(s_sig,aes(x=par,y=value,col=Model,fill=Model))+geom_boxplot(alpha=.5,outlier.shape=NA)+theme_bw()+theme(panel.grid=element_blank(),axis.text.x = element_text(angle=45, hjust=1))+scale_y_continuous(limits = c(0,.75))
dev.off()


####calculate  AIC WEIGHTED PARAMETERS: ###http://treethinkers.org/tutorials/morphological-evolution-in-r/  ####

PABM_wt_sum<-subset(all_mods,par=="PA");PABM_wt_sum<-subset(PABM_wt_sum,Model=="BM");PABM_wt_sum$wtRate<-PABM_wt_sum$AIC*PABM_wt_sum$rate;PABM_wt_sum<-sum(PABM_wt_sum$wtRate)
PAOU_wt_sum<-subset(all_mods,par=="PA");PAOU_wt_sum<-subset(PAOU_wt_sum,Model=="OU")
PAEB_wt_sum<-subset(all_mods,par=="PA");PAEB_wt_sum<-subset(PAEB_wt_sum,Model=="EB")



weightedSum<-function(matrix,PAR,MODEL){
	out<-matrix(nrow=1,ncol=3)
	temp<-subset(matrix,par==PAR)
	temp<-subset(temp,Model==MODEL)
	wtrate<-temp$AIC*temp$rate
	wtSum<-sum(wtrate)
	out[1,1]<-PAR
	out[1,2]<-MODEL
	out[1,3]<-wtSum
	out<-data.frame(out)
	colnames(out)<-c("par","Model","value")
	return(out)
	}

PA_BM_ws<-weightedSum(all_mods,"PA","BM")
PA_OU_ws<-weightedSum(all_mods,"PA","OU")
PA_EB_ws<-weightedSum(all_mods,"PA","EB")

PF_BM_ws<-weightedSum(all_mods,"PF","BM")
PF_OU_ws<-weightedSum(all_mods,"PF","OU")
PF_EB_ws<-weightedSum(all_mods,"PF","EB")

SD_BM_ws<-weightedSum(all_mods,"SD","BM")
SD_OU_ws<-weightedSum(all_mods,"SD","OU")
SD_EB_ws<-weightedSum(all_mods,"SD","EB")

TD_BM_ws<-weightedSum(all_mods,"TD","BM")
TD_OU_ws<-weightedSum(all_mods,"TD","OU")
TD_EB_ws<-weightedSum(all_mods,"TD","EB")

wtSums<-rbind(PA_BM_ws,PA_OU_ws,PA_EB_ws,PF_BM_ws,PF_OU_ws,PF_EB_ws,SD_BM_ws,SD_OU_ws,SD_EB_ws,TD_BM_ws,TD_OU_ws,TD_EB_ws)
wtSums$value<-as.numeric(as.character(wtSums$value))

pdf("~/Dropbox/Warbler.Molt.Migration/AICweighted_rateParameterSums.pdf")
ggplot(wtSums,aes(x=par,y=log(value),col=Model))+geom_point(aes(size=log(value),alpha=log(value)))+theme_bw()+ggtitle("AIC weighted parameter rates sums")+theme(plot.title = element_text(hjust = 0.5))
dev.off()


##Also maybe plot parameter rates with size of dot weighted by AICc weight

####Phylogenetic signal

pSigMatrix<-function(matrix,phy){
	tree.dat<-list()
	psig<-matrix(nrow=3,ncol=ncol(matrix))
	rownames(psig)<-c("Lambda","logL","P")
	colnames(psig)<-colnames(matrix)
	for(i in 1:ncol(matrix)){
		tree.dat[[i]]<-treedata(phy,na.omit(matrix[,i,drop=FALSE]))
			psig[1,i]<-phylosig(tree.dat[[i]]$phy,tree.dat[[i]]$data,test=TRUE,method="lambda")$lambda
			psig[2,i]<-phylosig(tree.dat[[i]]$phy,tree.dat[[i]]$data,test=TRUE,method="lambda")$logL
			psig[3,i]<-phylosig(tree.dat[[i]]$phy,tree.dat[[i]]$data,test=TRUE,method="lambda")$P
			}
			return(psig)
		}


PA_lambda<-pSigMatrix(PA,warTree)
write.csv(PA_lambda,"~/Dropbox/Warbler.Molt.Migration/PA_lambdas.csv")

PF_lambda<-pSigMatrix(PF,warTree)
write.csv(PF_lambda,"~/Dropbox/Warbler.Molt.Migration/PF_lambdas.csv")

SD_lambda<-pSigMatrix(SD,warTree)
write.csv(SD_lambda,"~/Dropbox/Warbler.Molt.Migration/SD_lambdas.csv")

TD_lambda<-pSigMatrix(TD,warTree)
write.csv(TD_lambda,"~/Dropbox/Warbler.Molt.Migration/TD_lambdas.csv")



PA_lam<-melt(PA_lambda);PA_lam$par<-"PA"
PF_lam<-melt(PF_lambda);PF_lam$par<-"PF"
SD_lam<-melt(SD_lambda);SD_lam$par<-"SD"
TD_lam<-melt(TD_lambda);TD_lam$par<-"TD"


all_lambdas<-rbind(PA_lam,PF_lam,SD_lam,TD_lam);colnames(all_lambdas)<-c("opt","variable","value","par")

pdf("~/Dropbox/Warbler.Molt.Migration/extent_Lambda.pdf")
ggplot(subset(all_lambdas,variable%in%c("extent")),aes(x=variable,y=value,col=opt,fill=opt))+geom_point(size=4)+facet_grid(opt~par,scales="free_y")+theme_bw()+theme(panel.grid=element_blank(),axis.text.x = element_text(angle=45, hjust=1))+ggtitle("Phylogenetic Signal")+theme(plot.title = element_text(hjust = 0.5))
dev.off()



pdf("~/Dropbox/Warbler.Molt.Migration/tract_Lambda.pdf")
ggplot(subset(all_lambdas,variable%in%tracts),aes(x=variable,y=value,col=opt,fill=opt))+geom_point(size=4)+facet_grid(opt~par,scales="free_y")+theme_bw()+theme(panel.grid=element_blank(),axis.text.x = element_text(angle=45, hjust=1,size=6))+ggtitle("Phylogenetic Signal")+theme(plot.title = element_text(hjust = 0.5))
dev.off()

pdf("~/Dropbox/Warbler.Molt.Migration/mods_Lambda.pdf")
ggplot(subset(all_lambdas,variable%in%mods),aes(x=variable,y=value,col=opt,fill=opt))+geom_point(size=4)+facet_grid(opt~par,scales="free_y")+theme_bw()+theme(panel.grid=element_blank(),axis.text.x = element_text(angle=45, hjust=1,size=6))+ggtitle("Phylogenetic Signal")+theme(plot.title = element_text(hjust = 0.5))
dev.off()



####compare ER and ARD models
###Function saves matrix of log likelihoods for ER, and ARD models
### calculate p value for ARD with likelihood ratio test 1-pchisq(2*abs(ERloglik-ARDloglik,1)


###function that calculates probabilities of ARD model


pARD<-function(matrix,phy,TYPE){
	#TYPE = "discrete" or "continuous"
	tree.dat<-list()
	logLiks<-matrix(nrow=3,ncol=ncol(matrix))
	rownames(logLiks)<-c("ER","ARD","pARD")
	colnames(logLiks)<-colnames(matrix)
	for(i in 1:ncol(matrix)){
		tree.dat[[i]]<-treedata(phy,na.omit(matrix[,i,drop=FALSE]))
			logLiks[1,i]<-ace(x=tree.dat[[i]]$data,phy=tree.dat[[i]]$phy,type=TYPE,model="ER")$loglik
			logLiks[2,i]<-ace(x=tree.dat[[i]]$data,phy=tree.dat[[i]]$phy,type=TYPE,model="ARD")$loglik
			logLiks[3,i]<-(1-pchisq(2*abs(logLiks[1,i]-logLiks[2,i]),1))

			}
			return(logLiks)
		}
		
##function that returns transition paramaters under best model####		

PApresPARD<-pARD(PA[,1,drop=FALSE],warTree,"discrete")
PFpresPARD<-pARD(PF[,1,drop=FALSE],warTree,"discrete")
SDpresPARD<-pARD(SD[,1,drop=FALSE],warTree,"discrete")
TDpresPARD<-pARD(TD[,1,drop=FALSE],warTree,"discrete")

ARDprob<-matrix(nrow=4,ncol=1)
rownames(ARDprob)<-c("PA","PF","SD","TD")
colnames(ARDprob)<-"p(ARD)"
ARDprob[1,1]<-PApresPARD[3,1]
ARDprob[2,1]<-PFpresPARD[3,1]
ARDprob[3,1]<-SDpresPARD[3,1]
ARDprob[4,1]<-TDpresPARD[3,1]


