###threshold modelling
###From Liam Revell's tutorial here: http://lukejharmon.github.io/ilhabela/instruction/2015/07/05/threshold-model/
#### Modelling exogenous correlates of PA and seasonal dichromatism under threshold model

PA_tr<-read.csv("/Users/ryanterrill/Dropbox/Warbler.Molt.Migration/Molt\ extent/Parulidae_PA.csv",row.names=1)
PF_tr<-read.csv("/Users/ryanterrill/Dropbox/Warbler.Molt.Migration/Molt\ extent/Parulidae_PF.csv",row.names=1)
SD_tr<-read.csv("~/Dropbox/Warbler.Molt.Migration/Parulidae_sexual_dimorphism.csv",row.names=1)
TD_tr<-read.csv("~/Dropbox/Warbler.Molt.Migration/Parulidae_temporal_dimorphism.csv",row.names=1)
Data<-read.csv("~/Dropbox/Warbler.Molt.Migration/Data_treeNames.csv",row.names=1)
Data$names<-rownames(Data)

library(phytools)
require(gridExtra)


warTree<-read.nexus("/Users/ryanterrill/Dropbox/Warbler.Molt.Migration/Parulidae_phylogeny/Trees/r8stree.tre")
warTree<-multi2di(warTree)
warTree<-chronos(warTree)
class(warTree)<-"phylo"

warTree

##Transform PA tracts to binary discrete character
PA_tr[PA_tr==.5]<-1
TD_tr[TD_tr==.5]<-1




###tree and vars from tutorial
tree<-pbtree(n=50,scale=1)
r<-0.75 # simulate using a high correlation
V<-matrix(c(1,r,r,1),2,2)
X<-sim.corrs(tree,V)

colnames(PA_tr)<-sub("PA_","",names(PA_tr))
pa.dat<-merge(PA_tr,Data,by=0)
rownames(pa.dat)<-pa.dat$Row.names
pa.dat<-pa.dat[,-1]

paTree<-treedata(warTree,pa.dat)


####Plot presence of PA by tract with migration


paMig<-data.frame(pa.dat[,1:18])
paMig<-paMig[,-(13:17)]
paMig<-paMig[complete.cases(paMig),]
paMig.melt<-melt(paMig,id="migDist")

pdf("~/Dropbox/Warbler.Molt.Migration/Migratory_distance_PATracts.pdf")
ggplot(paMig.melt,aes(y=migDist,x=factor(value),fill=factor(value)))+geom_boxplot(alpha=.3)+geom_jitter(aes(colour=factor(value)))+facet_wrap(~variable)+theme_bw()+labs(x="Presence of Prealternate Molt",y="Migratory Distance")
dev.off()


colnames(TD_tr)<-sub("TD_","",names(TD_tr))
td.dat<-merge(TD_tr,Data,by=0)
rownames(td.dat)<-td.dat$Row.names
td.dat<-td.dat[,-1]



####Plot presence of td by tract with migration


tdMig<-data.frame(td.dat[,1:23])
tdMig<-tdMig[,-(13:22)]
tdMig<-tdMig[complete.cases(tdMig),]
tdMig.melt<-melt(tdMig,id="breeding_forage_stratum")

pdf("~/Stratum_tdTracts.pdf")
ggplot(tdMig.melt,aes(y=breeding_forage_stratum,x=factor(value),fill=factor(value)))+geom_boxplot(alpha=.3)+geom_jitter(aes(colour=factor(value)))+facet_wrap(~variable)+theme_bw()+labs(x="Presence of Seasonal Dichromatism",y="Foraging Stratum")
dev.off()

paMig.melt$variable <- factor(paMig.melt$variable, levels = c("head","breast","back","belly","med_cov","gr_cov","tertials","secondaries","rectrices","primaries","primary_coverts","alula"))
tdMig.melt$variable <- factor(tdMig.melt$variable, levels = c("head","breast","back","belly","med_cov","gr_cov","tertials","secondaries","rectrices","primaries","primary_coverts","alula"))



pm<-ggplot(paMig.melt,aes(y=migDist,x=factor(value),color=factor(value)))+geom_boxplot(outlier.shape=NA,color="black",alpha=.3,lwd=.1,aes(fill=factor(value)))+geom_jitter(shape=21,width=.1,height=.1,color="black",aes(fill=factor(value)),size=.9)+facet_grid(.~variable)+theme_bw()+labs(x="Presence of Prealternate Molt",y="Migratory Distance")+scale_fill_grey(start=1,end=0)+theme(legend.position="none")+theme(strip.text.x = element_text(size = 4),axis.text=element_text(size=5))
th<-ggplot(tdMig.melt,aes(y=breeding_forage_stratum,x=factor(value),color=factor(value)))+geom_boxplot(outlier.shape=NA,color="black",alpha=.3,lwd=.1,aes(fill=factor(value)))+geom_jitter(shape=21,width=.1,height=.1,color="black",aes(fill=factor(value)),size=.9)+facet_grid(.~variable)+theme_bw()+labs(x="Presence of Seasonal Dichromatism",y="Foraging Stratum")+scale_fill_grey(start=1,end=0)+theme(legend.position="none")+theme(strip.text.x = element_text(size = 4),axis.text=element_text(size=5))

pdf("~/Dropbox/Warbler.Molt.Migration/PA_TD_tractPredictors.pdf")
grid.arrange(pm,th, nrow=3)
dev.off()


######
This looks great - run phylogenetic ANOVA by tract for p-values

####
ANOVA for tracts PA- Migdist and DayLength - TD and stratum (plus any other significant non autocorrellated models)
###test anova

paHead<-paMig$head;names(paHead)<-rownames(paMig)
paHead<-factor(paHead,levels=c(0,1))
paMigDist<-paMig$migDist;names(paMigDist)<-rownames(paMig)

anovaTree<-drop.tip(warTree,setdiff(warTree$tip.label,rownames(paMig)))

paHeadanova<-phylANOVA(anovaTree,x=paHead,y=paMigDist,nsim=1000,posthoc=T)

#migratory distance vs. PA in feather regions ANOVAs

anovaMatrix<-function(responseMatrix,predictor,phy){
		responseMatrix<-responseMatrix[complete.cases(responseMatrix),]
		responseMatrix<-responseMatrix[,colSums(responseMatrix)!=0]
		tree<-drop.tip(phy,setdiff(phy$tip.label,rownames(responseMatrix)))
		out<-list()
		out$params<-list()
		out$p<-matrix(nrow=ncol(responseMatrix),ncol=1)
		out$f<--matrix(nrow=ncol(responseMatrix),ncol=1)
		rM<-list()
		for(i in 1:ncol(responseMatrix)){
			rM[[i]]<-responseMatrix[,i]
				names(rM)[i]<-colnames(responseMatrix)[i]
				names(rM[[i]])<-rownames(responseMatrix)
			out$params[[i]]<-phylANOVA(tree,x=factor(rM[[i]],levels=c(0,1)),y=predictor,posthoc=TRUE)
			out$p[i,]<-out$params[[i]]$Pf
			out$f[i,]<-out$params[[i]]$F
			}
		names(out$params)<-colnames(responseMatrix)
		rownames(out$p)<-colnames(responseMatrix)
		return(out)
		}
			
	##put p values in output matrix by response and predictor
	
	
paMiganova<-anovaMatrix(PA_tr,paMigDist,warTree)
write.csv(paMiganova$p,"~/Dropbox/Warbler.Molt.Migration/PAMigANOVA.csv")


tdStrat<-tdMig[,13];names(tdStrat)<-rownames(tdMig)

tdStratAnova<-anovaMatrix(tdMig[,-13],tdStrat,warTree)
write.csv(tdStratAnova$p,"~/Dropbox/Warbler.Molt.Migration/TDstratANOVA.csv")

###condust ANOVAs for tracts for PA and TD for top pgls models, including mixed models



#####td tracts thresholds######

pa.tr.c<-PA_tr[complete.cases(PA_tr),]
pa.tr.tree<-drop.tip(warTree,setdiff(warTree$tip.label,rownames(pa.tr.c)))

mcmcPA.tr<-threshBayes(pa.tr.tree,pa.tr.c,ngen=ngen,control=list(sample=sample))


##test run - just checking R value of correlation using MCMC


sample<-500
ngen<-200000 ## in 'real' studies this should be larger
burnin<-0.2*ngen



#####PAHead ~ migDist######
pa.mig.head<-as.matrix(data.frame(pa.dat$PA_head,pa.dat$migDist))
rownames(pa.mig.head)<-rownames(pa.dat)
pa.mig.head<-pa.mig.head[complete.cases(pa.mig.head),]
pa.mig.tree<-drop.tip(warTree,setdiff(warTree$tip.label,rownames(pa.mig.head)))



test<-data.frame(colMeans(mcmcPAHead$liab[,2:52]),colMeans(mcmcPAHead$liab[,53:103]))
test1<-data.frame(colMeans(mcmcPAHead$liab[,2:52]),pa.mig.head[,1])

##test run - just checking R value of correlation using MCMC


sample<-500
ngen<-200000 ## in 'real' studies this should be larger
burnin<-0.2*ngen



mcmcPAHead<-threshBayes(pa.mig.tree,pa.mig.head,types=c("disc","cont"),ngen=ngen,control=list(sample=sample))
mean(mcmcPAHead$par[(burnin/sample+1):nrow(mcmcCC$par),"r"])



##plot pesterior density for correlation

pdf("~/Dropbox/Warbler.Molt.Migration/PA-head~migDist_corrPost.pdf")

plot(density(mcmcPAHead$par[(burnin/sample+1):nrow(mcmcCC$par),
    "r"],bw=0.1),xlab="r",main="posterior density for r - head PA vs migratory distance")
lines(c(r,r),c(0,1000),lty="dashed")
dev.off()


## plot our likelihood profile
pdf("~/Dropbox/Warbler.Molt.Migration/PA-head~migDist_threshold_PosteriorProbability.pdf")
plot(density(mcmcPAHead$par[,4]),main="Migratory distance threshold for PA - head",xlab="degrees latitude travelled",ylab="Transition Density 0->1")
lines(c(mean(mcmcPAHead$par[,4]),.5),c(0,1000),lty="dashed")

plot(density(mcmcPAHead$par[,4]*mcmcPAHead$par[,5]))

dev.off()

###plot threshold location
###Great, this looks right - calculate threshold for all tracts for PA

#####PA_breast ~ migDist######
pa.mig.breast<-as.matrix(data.frame(pa.dat$PA_breast,pa.dat$migDist))
rownames(pa.mig.breast)<-rownames(pa.dat)
pa.mig.breast<-pa.mig.breast[complete.cases(pa.mig.breast),]
pa.mig.tree<-drop.tip(warTree,setdiff(warTree$tip.label,rownames(pa.mig.breast)))
mcmcPAbreast<-threshBayes(pa.mig.tree,pa.mig.breast,types=c("disc","cont"),ngen=ngen,control=list(sample=sample))

plot(density(mcmcPAbreast$par[,5]),main="Migratory distance threshold for PA - breast",xlab="degrees latitude travelled",ylab="Transition Density 0->1")
lines(c(mean(mcmcPAbreast$par[,5]),.5),c(0,1000),lty="dashed")

#####PA_belly ~ migDist######
pa.mig.belly<-as.matrix(data.frame(pa.dat$PA_belly,pa.dat$migDist))
rownames(pa.mig.belly)<-rownames(pa.dat)
pa.mig.belly<-pa.mig.belly[complete.cases(pa.mig.belly),]
pa.mig.tree<-drop.tip(warTree,setdiff(warTree$tip.label,rownames(pa.mig.belly)))
mcmcPAbelly<-threshBayes(pa.mig.tree,pa.mig.belly,types=c("disc","cont"),ngen=ngen,control=list(sample=sample))


plot(density(mcmcPAbelly$par[,5]),main="Migratory distance threshold for PA - belly",xlab="degrees latitude travelled",ylab="Transition Density 0->1")
lines(c(mean(mcmcPAbelly$par[,5]),.5),c(0,1000),lty="dashed")


#####PA_back ~ migDist######
pa.mig.back<-as.matrix(data.frame(pa.dat$PA_back,pa.dat$migDist))
rownames(pa.mig.back)<-rownames(pa.dat)
pa.mig.back<-pa.mig.back[complete.cases(pa.mig.back),]
pa.mig.tree<-drop.tip(warTree,setdiff(warTree$tip.label,rownames(pa.mig.back)))
mcmcPAback<-threshBayes(pa.mig.tree,pa.mig.back,types=c("disc","cont"),ngen=ngen,control=list(sample=sample))



plot(density(mcmcPAback$par[,6]),main="PA in Back ~ Migratory distance",xlab="correlation coefficient",ylab="Density",xlim=c(-.5,1))
lines(c(mean(mcmcPAback$par[,6]),.5),c(0,1000),lty="dashed")


#####PA_rectrices ~ migDist######
pa.mig.rectrices<-as.matrix(data.frame(pa.dat$PA_rectrices,pa.dat$migDist))
rownames(pa.mig.rectrices)<-rownames(pa.dat)
pa.mig.rectrices<-pa.mig.rectrices[complete.cases(pa.mig.rectrices),]
pa.mig.tree<-drop.tip(warTree,setdiff(warTree$tip.label,rownames(pa.mig.rectrices)))
mcmcPArectrices<-threshBayes(pa.mig.tree,pa.mig.rectrices,types=c("disc","cont"),ngen=ngen,control=list(sample=sample))


plot(density(mcmcPArectrices$par[,6]),main="PA in Rectrices ~ Migratory distance",xlab="correlation coefficient",ylab="Density",xlim=c(-.5,1))
lines(c(mean(mcmcPArectrices$par[,6]),.5),c(0,1000),lty="dashed")




###Threshold Ancestral state and transition estimation


####so make matrix for R valus for thresholds for each tract by top few models
#### Then look at location of thresholds for each


