###Ancestral state reconstructions and stochastic character maps

library(phytools)
library(geiger)
library(reshape)
library(ggplot2)
library(ggthemes)



#############################################################################
#######    Same set-up as model_selection.R      ############################
#####copied here so scripts can be run independently of each other ##########
#############################################################################

oneMat<-matrix(nrow=51,ncol=2)
oneMat[,1]<-rep(0,51)
oneMat[,2]<-rep(1,51)
colnames(oneMat)<-c("0","1")

zeroMat<-matrix(nrow=51,ncol=2)
zeroMat[,1]<-rep(1,51)
zeroMat[,2]<-rep(0,51)
colnames(zeroMat)<-c("0","1")


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

###################################################################################
########### Ancestral State Reconstruction and Stochastic Character Mapping #######
###################################################################################



###### Discrete ACEs - presence of molts and dichromatism, and presence by body part and modular unit###
###### Convert all integers >1 to 1 ######


###get likelihoods of root state and rate paramatersfor each; for best model -- ARD for PA and TD, ER for SD

###plot each variable and by tract

pres<-Data[,1:4,drop=FALSE]
pres<-pres[complete.cases(pres),]

pres[pres>0]<-"P"
pres[pres==0]<-"A"
presTree<-drop.tip(warTree,setdiff(warTree$tip.label,rownames(pres)))


###Order all the data up here######
##################################
presnum<-data.frame(1:51);rownames(presnum)<-presTree$tip.label

presnum<-merge(pres,presnum,by=0)
presnum<-presnum[order(presnum$X1.51),]
pres<-presnum[,2:5];rownames(pres)<-presnum$Row.names

cols<-c("white","black");names(cols)<-c("0","1")

###PA Presence#####

##########
###Add rates of gains and losses to these plots - plot separately and combine in affinity###
##########
PA_ex<-pres$paextent;names(PA_ex)<-rownames(pres)

PApresER<-fitDiscrete(presTree,PA_ex,model="ER")
PApresARD<-fitDiscrete(presTree,PA_ex,model="ARD")
PA_ex_lik<-1-pchisq(2*(PApresARD$opt$lnL-2*PApresER$opt$lnL),1)### export these liklihood ratio test P values to a table
PApres<-ace(PA_ex,presTree,model="ARD",type="discrete") ## export rates of all of these to table, or plot with std err

pdf("~/Dropbox/Warbler.Molt.Migration/PApresence.pdf")
plot(presTree,cex=.5,main = "Presence of Prealternate Molt",label.offset=.01)
cols<-c("white","black");names(cols)<-c(0,1)
nodelabels(node=1:tree$Nnode+Ntip(tree),
    pie=PApres$lik.anc,piecol=cols,cex=0.5)
tiplabels(pie=to.matrix(PA_ex,sort(unique(PA_ex))),piecol=cols,cex=0.3,offset=-.1)
dev.off()

###SD Presence#####
SD_ex<-pres$sdExtent;names(SD_ex)<-rownames(pres)
SDpresER<-fitDiscrete(presTree,SD_ex,model="ER")
SDpresARD<-fitDiscrete(presTree,SD_ex,model="ARD")
SD_ex_lik<-1-pchisq(2*(SDpresARD$opt$lnL-2*SDpresER$opt$lnL),1)### export these liklihood ratio test P values to a table
SDpres<-ace(SD_ex,presTree,model="ARD",type="discrete")

pdf("~/Dropbox/Warbler.Molt.Migration/SDpresence.pdf")
plot(presTree,cex=.5,main = "Presence of Sexual Dichromatism",label.offset=.01)
nodelabels(node=1:tree$Nnode+Ntip(tree),
    pie=SDpres$lik.anc,piecol=cols,cex=0.5)
tiplabels(pie=to.matrix(SD_ex,sort(unique(SD_ex))),piecol=cols,cex=0.3,offset=-.1)
dev.off()


###TD Presence#####

TD_ex<-pres$tdExtent;names(TD_ex)<-rownames(pres)
TDpresER<-fitDiscrete(presTree,TD_ex,model="ER")
TDpresARD<-fitDiscrete(presTree,TD_ex,model="ARD")
TD_ex_lik<-1-pchisq(2*(TDpresARD$opt$lnL-2*TDpresER$opt$lnL),1)### export these liklihood ratio test P values to a table
TD_ex_lik
TDpres<-ace(TD_ex,presTree,model="ARD",type="discrete")

TDnum<-data.frame(1:51);rownames(TDnum)<-presTree$tip.label


pdf("~/Dropbox/Warbler.Molt.Migration/TDpresence.pdf")
plot(presTree,cex=.5,main = "Presence of Seasonal Dichromatism",label.offset=.01, show.node.label=TRUE)
nodelabels(node=1:presTree$Nnode+Ntip(presTree),
    pie=TDpres$lik.anc,piecol=cols,cex=0.5)
tiplabels(pie = to.matrix(TD_ex, sort(unique(TD_ex))), piecol = cols, 
    cex = 0.3)
dev.off()


#################################################
###### Figure for PA by feather region ############
#################################################


PA_tr<-read.csv("/Users/ryanterrill/Dropbox/Warbler.Molt.Migration/Molt\ extent/Parulidae_PA.csv",row.names=1)

PA_or<-data.frame(1:51);rownames(PA_or)<-presTree$tip.label

PA_tr<-merge(PA_or,PA_tr,by=0)
PA_tr<-PA_tr[order(PA_tr$X1.51),]
rownames(PA_tr)<-PA_tr$Row.names;PA_tr<-PA_tr[,3:14]

PA_tr[PA_tr==.5]<-1
PA_tr[PA_tr==.6]<-1

hPA<-PA_tr[,1];names(hPA)<-rownames(PA_tr)
headPA<-ace(hPA,presTree,model="ER",type="discrete")

brPA<-PA_tr[,2];names(brPA)<-rownames(PA_tr)
breastPA<-ace(brPA,presTree,model="ER",type="discrete")

bePA<-PA_tr[,3];names(bePA)<-rownames(PA_tr)
bellyPA<-ace(bePA,presTree,model="ER",type="discrete")

baPA<-PA_tr[,4];names(baPA)<-rownames(PA_tr)
backPA<-ace(baPA,presTree,model="ER",type="discrete")

rePA<-PA_tr[,5];names(rePA)<-rownames(PA_tr)
rectPA<-ace(rePA,presTree,model="ER",type="discrete")


###primaries == all absent
prPA<-PA_tr[,6];names(prPA)<-rownames(PA_tr)


sePA<-PA_tr[,7];names(sePA)<-rownames(PA_tr)
secPA<-ace(sePA,presTree,model="ER",type="discrete")

trPA<-PA_tr[,8];names(trPA)<-rownames(PA_tr)
tertPA<-ace(trPA,presTree,model="ER",type="discrete")

mcPA<-PA_tr[,9];names(mcPA)<-rownames(PA_tr)
medcovPA<-ace(mcPA,presTree,model="ER",type="discrete")

gcPA<-PA_tr[,10];names(gcPA)<-rownames(PA_tr)
grcovPA<-ace(gcPA,presTree,model="ER",type="discrete")

###p covs == all absent
pcPA<-PA_tr[,11];names(pcPA)<-rownames(PA_tr)

###alula = all absent
alPA<-PA_tr[,12];names(alPA)<-rownames(PA_tr)

#####Plot by feather region
pdf("~/Dropbox/Warbler.Molt.Migration/PA_ACE_byTract.pdf")

par(mar=c(1,4,3,10))

tcex=.3

plot(presTree,cex=.5,main = "Prealternate Molt by Feather Region",label.offset=.32, show.node.label=TRUE)
nodelabels(pie=headPA$lik.anc,piecol=cols,adj=c(0.48,0.5),cex=0.2)
nodelabels(pie=breastPA$lik.anc,piecol=cols,adj=c(0.49,0.5),cex=0.2)
nodelabels(pie=bellyPA$lik.anc,piecol=cols,adj=c(0.5,0.5),cex=0.2)
nodelabels(pie=backPA$lik.anc,piecol=cols,adj=c(0.51,0.5),cex=0.2)
nodelabels(pie=rectPA$lik.anc,piecol=cols,adj=c(0.52,0.5),cex=0.2)
nodelabels(pie = to.matrix(prPA, sort(unique(hPA))),piecol=cols,adj=c(0.53,0.5),cex=0.2)
nodelabels(pie=secPA$lik.anc,piecol=cols,adj=c(0.54,0.5),cex=0.2)
nodelabels(pie=tertPA$lik.anc,piecol=cols,adj=c(0.55,0.5),cex=0.2)
nodelabels(pie=medcovPA$lik.anc,piecol=cols,adj=c(0.56,0.5),cex=0.2)
nodelabels(pie=grcovPA$lik.anc,piecol=cols,adj=c(0.57,0.5),cex=0.2)
nodelabels(pie = to.matrix(pcPA, sort(unique(hPA))),piecol=cols,adj=c(0.58,0.5),cex=0.2)
nodelabels(pie = to.matrix(alPA, sort(unique(hPA))),piecol=cols,adj=c(0.59,0.5),cex=0.2)



tiplabels(pie = to.matrix(hPA, sort(unique(hPA))), piecol = cols,adj=c(.58,0.5), cex = 0.3)
tiplabels(pie = to.matrix(brPA, sort(unique(hPA))), piecol = cols,adj=c(.60,0.5), cex = 0.3)
tiplabels(pie = to.matrix(bePA, sort(unique(bePA))), piecol = cols,adj=c(.62,0.5), cex = 0.3)
tiplabels(pie = to.matrix(baPA, sort(unique(baPA))), piecol = cols,adj=c(.64,0.5), cex = 0.3)
tiplabels(pie = to.matrix(rePA, sort(unique(rePA))), piecol = cols,adj=c(.66,0.5), cex = 0.3)
tiplabels(pie = to.matrix(prPA, sort(unique(prPA))), piecol = cols,adj=c(.68,0.5), cex = 0.3)
tiplabels(pie = to.matrix(sePA, sort(unique(sePA))), piecol = cols,adj=c(.70,0.5), cex = 0.3)
tiplabels(pie = to.matrix(trPA, sort(unique(trPA))), piecol = cols,adj=c(.72,0.5), cex = 0.3)
tiplabels(pie = to.matrix(mcPA, sort(unique(mcPA))), piecol = cols,adj=c(.74,0.5), cex = 0.3)
tiplabels(pie = to.matrix(gcPA, sort(unique(gcPA))), piecol = cols,adj=c(.76,0.5), cex = 0.3)
tiplabels(pie = to.matrix(pcPA, sort(unique(pcPA))), piecol = cols,adj=c(.78,0.5), cex = 0.3)
tiplabels(pie = to.matrix(alPA, sort(unique(alPA))),piecol=cols,adj=c(0.80,0.5),cex=0.3)

dev.off()


###### 
#plot mean and SD for transition rates for each region

###### 
#plot mean and SD for transition rates for each region

rateMat<-matrix(nrow=12,ncol=3)
rownames(rateMat)<-colnames(PA_tr)
colnames(rateMat)<-c("rate","se","reg")


##set up matrix with rates and SEs
rateMat[1,1]<-headPA$rates
rateMat[1,2]<-headPA$se

rateMat[2,1]<-breastPA$rates[1]
rateMat[2,2]<-breastPA$se[1]

rateMat[3,1]<-bellyPA$rates[1]
rateMat[3,2]<-bellyPA$se[1]


rateMat[4,1]<-backPA$rates[1]
rateMat[4,2]<-backPA$se[1]


rateMat[5,1]<-rectPA$rates[1]
rateMat[5,2]<-rectPA$se[1]


rateMat[6,1]<-0
rateMat[6,2]<-NaN

rateMat[7,1]<-secPA$rates[1]
rateMat[7,2]<-secPA$se[1]


rateMat[8,1]<-tertPA$rates[1]
rateMat[8,2]<-tertPA$se[1]

rateMat[9,1]<-medcovPA$rates
rateMat[9,2]<-medcovPA$se


rateMat[10,1]<-grcovPA$rates[1]
rateMat[10,2]<-grcovPA$se[1]

rateMat[11,1]<-0
rateMat[11,2]<-NaN

rateMat[12,1]<-0
rateMat[12,2]<-NaN



rateMat[,3]<-rownames(rateMat)

rateMat<-data.frame(rateMat)

####need to change matrix so mean and sd are columns


rateMat$reg<-gsub("PA_","",rateMat$reg)
rateMat$reg<-gsub("rectrices","rr",rateMat$reg)
rateMat$reg<-gsub("primaries","pp",rateMat$reg)
rateMat$reg<-gsub("secondaries","ss",rateMat$reg)
rateMat$reg<-gsub("med_cov","m_cov",rateMat$reg)
rateMat$reg<-gsub("primary_coverts","p_cov",rateMat$reg)
rateMat[is.na(rateMat)]<-0

rateMat$reg<-factor(rateMat$reg, level = c("head","breast","belly","back","rr","pp","ss","tertials","m_cov","gr_cov","p_cov","alula"))

limits <- aes(ymax = as.numeric(as.character(rate)) + as.numeric(as.character(se)), ymin=as.numeric(as.character(rate)) - as.numeric(as.character(se)))

write.csv(t(rateMat),"~/Dropbox/Warbler.Molt.Migration/PA_tractRates.csv")

pdf("~/Dropbox/Warbler.Molt.Migration/PA_tracttransitions.pdf")

ggplot(rateMat,aes(y=as.numeric(as.character(rate)),x=reg))+geom_errorbar(limits)+geom_point(size=6)+theme_tufte()+theme(panel.border=element_rect(fill="transparent"),axis.text.x = element_text(size=12,angle=45,vjust=.5),plot.margin=unit(c(1,1,12,1),"cm"))+scale_shape_manual(values=c(16,1))+xlab("rate of evolution")+ylab("")

dev.off()



lamMat<-matrix(nrow=14,ncol=2)
rownames(lamMat)<-colnames(PA_tr)
colnames(lamMat)<-c("lambda","reg")

lamMat[1,1]<-phylosig(presTree,hPA)
lamMat[2,1]<-phylosig(presTree,brPA)
lamMat[3,1]<-phylosig(presTree,bePA)
lamMat[4,1]<-phylosig(presTree,baPA)
lamMat[5,1]<-phylosig(presTree,rePA)
lamMat[6,1]<-phylosig(presTree,prPA)
lamMat[7,1]<-phylosig(presTree,sePA)
lamMat[8,1]<-phylosig(presTree,trPA)
lamMat[9,1]<-phylosig(presTree,mcPA)
lamMat[10,1]<-phylosig(presTree,gcPA)
lamMat[11,1]<-phylosig(presTree,pcPA)
lamMat[12,1]<-phylosig(presTree,alPA)
lamMat[13,1]<-phylosig(presTree,PA_ex)
lamMat[14,1]<-phylosig(warTree,PA_pr)




lamMat[,2]<-rownames(lamMat)
lamMat[6,1]<-1
lamMat[11,1]<-1
lamMat[12,1]<-1

lamMat<-data.frame(lamMat)

lamMat$reg<-gsub("PA_","",lamMat$reg)
lamMat$reg<-gsub("rectrices","rr",lamMat$reg)
lamMat$reg<-gsub("primaries","pp",lamMat$reg)
lamMat$reg<-gsub("secondaries","ss",lamMat$reg)
lamMat$reg<-gsub("med_cov","m_cov",lamMat$reg)
lamMat$reg<-gsub("primary_coverts","p_cov",lamMat$reg)

lamMat$reg<- c("head","breast","belly","back","rr","pp","ss","tertials","m_cov","gr_cov","p_cov","alula","extent","presence")

write.csv(lamMat,"~/Dropbox/Warbler.Molt.Migration/PAlambda.csv")

pdf("~/Dropbox/Warbler.Molt.Migration/phylsignal.pdf")
ggplot(lamMat,aes(y=as.numeric(as.character(lambda)),x=reg))+geom_point(size=4)+theme_tufte()+theme(panel.border=element_rect(fill="transparent"),axis.text.x = element_text(size=12,angle=45,vjust=.5),plot.margin=unit(c(1,1,12,1),"cm"))+scale_shape_manual(values=c(16,1))+xlab("phylogenetic signal")+ylab("")
dev.off()




PA_ex<-Data[,1];names(PA_ex)<-rownames(Data)
PA_ex<-PA_ex[complete.cases(PA_ex)]


pa.cm<-contMap(presTree,PA_ex,outline=TRUE,lwd=6)




pa.cm<-setMap(pa.cm,colors=c("white","black"))


pdf("~/Dropbox/Warbler.Molt.Migration/PAContMap.pdf")
plot(pa.cm,fsize=0)
dev.off()

pdf("~/Dropbox/Warbler.Molt.Migration/PAextentTree.pdf")
plotTree.wBars(pa.cm$tree,PA_ex,method="plotSimmap",colors=pa.cm$cols,scale=.05)
dev.off()

#################################################
###### Figure for PF by feather region ############
#################################################


PF_tr<-read.csv("/Users/ryanterrill/Dropbox/Warbler.Molt.Migration/Molt\ extent/parulidae_PF.csv",row.names=1)

PF_or<-data.frame(1:51);rownames(PF_or)<-presTree$tip.label

PF_tr<-merge(PF_or,PF_tr,by=0)
PF_tr<-PF_tr[order(PF_tr$X1.51),]
rownames(PF_tr)<-PF_tr$Row.names;PF_tr<-PF_tr[,3:14]

PF_tr[PF_tr==.5]<-1
PF_tr[PF_tr==.6]<-1

hPF<-PF_tr[,1];names(hPF)<-rownames(PF_tr)
headPF<-ace(hPF,presTree,model="ER",type="discrete")

brPF<-PF_tr[,2];names(brPF)<-rownames(PF_tr)
breastPF<-ace(brPF,presTree,model="ER",type="discrete")

bePF<-PF_tr[,3];names(bePF)<-rownames(PF_tr)
bellyPF<-ace(bePF,presTree,model="ER",type="discrete")

baPF<-PF_tr[,4];names(baPF)<-rownames(PF_tr)
backPF<-ace(baPF,presTree,model="ER",type="discrete")

rePF<-PF_tr[,5];names(rePF)<-rownames(PF_tr)
rectPF<-ace(rePF,presTree,model="ER",type="discrete")


prPF<-PF_tr[,6];names(prPF)<-rownames(PF_tr)
primPF<-ace(prPF,presTree,model="ER",type="discrete")


sePF<-PF_tr[,7];names(sePF)<-rownames(PF_tr)
secPF<-ace(sePF,presTree,model="ER",type="discrete")

trPF<-PF_tr[,8];names(trPF)<-rownames(PF_tr)
tertPF<-ace(trPF,presTree,model="ER",type="discrete")

med_covs - all present
mcPF<-PF_tr[,9];names(mcPF)<-rownames(PF_tr)


gcPF<-PF_tr[,10];names(gcPF)<-rownames(PF_tr)
grcovPF<-ace(gcPF,presTree,model="ER",type="discrete")


pcPF<-PF_tr[,11];names(pcPF)<-rownames(PF_tr)
prcovPF<-ace(pcPF,presTree,model="ER",type="discrete")


alPF<-PF_tr[,12];names(alPF)<-rownames(PF_tr)
alulPF<-ace(alPF,presTree,model="ER",type="discrete")

#####Plot by feather region
pdf("~/Dropbox/Warbler.Molt.Migration/PF_ACE_byTract.pdf")

par(mar=c(1,4,3,10))


plot(presTree,cex=.5,main = "Preformative Molt by Feather Region",label.offset=.32, show.node.label=TRUE)
nodelabels(pie=headPF$lik.anc,piecol=cols,adj=c(0.48,0.5),cex=0.2)
nodelabels(pie=breastPF$lik.anc,piecol=cols,adj=c(0.49,0.5),cex=0.2)
nodelabels(pie=bellyPF$lik.anc,piecol=cols,adj=c(0.50,0.5),cex=0.2)
nodelabels(pie=backPF$lik.anc,piecol=cols,adj=c(0.51,0.5),cex=0.2)
nodelabels(pie=rectPF$lik.anc,piecol=cols,adj=c(0.52,0.5),cex=0.2)
nodelabels(pie=primPF$lik.anc,piecol=cols,adj=c(0.53,0.5),cex=0.2)
nodelabels(pie=secPF$lik.anc,piecol=cols,adj=c(0.54,0.5),cex=0.2)
nodelabels(pie=tertPF$lik.anc,piecol=cols,adj=c(0.55,0.5),cex=0.2)
nodelabels(pie = to.matrix(mcPF, sort(unique(hPF))),piecol=cols,adj=c(0.56,0.5),cex=0.2)
nodelabels(pie=grcovPF$lik.anc,piecol=cols,adj=c(0.57,0.5),cex=0.2)
nodelabels(pie=prcovPF$lik.anc,piecol=cols,adj=c(0.58,0.5),cex=0.2)
nodelabels(pie=alulPF$lik.anc,piecol=cols,adj=c(0.59,0.5),cex=0.2)


oneMat<-matrix(nrow=51,ncol=2)
oneMat[,1]<-rep(0,51)
oneMat[,2]<-rep(1,51)
colnames(oneMat)<-c("0","1")

tiplabels(pie = oneMat, piecol = cols,adj=c(.58,0.5), cex=0.3)
tiplabels(pie = oneMat, piecol = cols,adj=c(.60,0.5), cex=0.3)
tiplabels(pie = oneMat, piecol = cols,adj=c(.62,0.5), cex=0.3)
tiplabels(pie = oneMat, piecol = cols,adj=c(.64,0.5), cex=0.3)
tiplabels(pie = to.matrix(rePF, sort(unique(rePF))), piecol = cols,adj=c(.66,0.5), cex=0.3)
tiplabels(pie = to.matrix(prPF, sort(unique(prPF))), piecol = cols,adj=c(.68,0.5), cex=0.3)
tiplabels(pie = to.matrix(sePF, sort(unique(sePF))), piecol = cols,adj=c(.70,0.5), cex=0.3)
tiplabels(pie = to.matrix(trPF, sort(unique(trPF))), piecol = cols,adj=c(.72,0.5), cex=0.3)
tiplabels(pie = to.matrix(mcPF, sort(unique(mcPF))), piecol = cols,adj=c(.74,0.5), cex=0.3)
tiplabels(pie = to.matrix(gcPF, sort(unique(gcPF))), piecol = cols,adj=c(.76,0.5), cex=0.3)
tiplabels(pie = to.matrix(pcPF, sort(unique(pcPF))), piecol = cols,adj=c(.78,0.5), cex=0.3)
tiplabels(pie = to.matrix(alPF, sort(unique(alPF))),piecol=cols,adj=c(.80,0.5),cex=0.3)

dev.off()


###### 
plot mean and SD for transition rates for each region

rateMat<-matrix(nrow=12,ncol=3)
rownames(rateMat)<-colnames(PF_tr)
colnames(rateMat)<-c("rate","se","reg")


##set up matrix with rates and SEs
rateMat[1,1]<-headPF$rates
rateMat[1,2]<-headPF$se

rateMat[2,1]<-breastPF$rates[1]
rateMat[2,2]<-breastPF$se[1]

rateMat[3,1]<-bellyPF$rates[1]
rateMat[3,2]<-bellyPF$se[1]


rateMat[4,1]<-backPF$rates[1]
rateMat[4,2]<-backPF$se[1]


rateMat[5,1]<-rectPF$rates[1]
rateMat[5,2]<-rectPF$se[1]


rateMat[6,1]<-primPF$rates[1]
rateMat[6,2]<-primPF$se[1]

rateMat[7,1]<-secPF$rates[1]
rateMat[7,2]<-secPF$se[1]


rateMat[8,1]<-tertPF$rates[1]
rateMat[8,2]<-tertPF$se[1]

rateMat[9,]<-rep(as.numeric(0),3)


rateMat[10,1]<-grcovPF$rates[1]
rateMat[10,2]<-grcovPF$se[1]

rateMat[11,1]<-prcovPF$rates[1]
rateMat[11,2]<-prcovPF$se[1]

rateMat[12,1]<-alulPF$rates[1]
rateMat[12,2]<-alulPF$se[1]



rateMat[,3]<-rownames(rateMat)

rateMat<-data.frame(rateMat)

####need to change matrix so mean and sd are columns


rateMat$reg<-gsub("PF_","",rateMat$reg)
rateMat$reg<-gsub("rectrices","rr",rateMat$reg)
rateMat$reg<-gsub("primaries","pp",rateMat$reg)
rateMat$reg<-gsub("secondaries","ss",rateMat$reg)
rateMat$reg<-gsub("med_cov","m_cov",rateMat$reg)
rateMat$reg<-gsub("primary_coverts","p_cov",rateMat$reg)

rateMat$reg<-factor(rateMat$reg, level = c("head","breast","belly","back","rr","pp","ss","tertials","m_cov","gr_cov","p_cov","alula"))

limits <- aes(ymax = as.numeric(as.character(rate)) + as.numeric(as.character(se)), ymin=as.numeric(as.character(rate)) - as.numeric(as.character(se)))

write.csv(t(rateMat),"~/Dropbox/Warbler.Molt.Migration/PF_tractRates.csv")


pdf("~/Dropbox/Warbler.Molt.Migration/PF_tracttransitions.pdf")

ggplot(rateMat,aes(y=as.numeric(as.character(rate)),x=reg))+geom_errorbar(limits)+geom_point(size=6)+theme_tufte()+theme(panel.border=element_rect(fill="transparent"),axis.text.x = element_text(size=12,angle=45,vjust=.5),plot.margin=unit(c(1,1,12,1),"cm"))+scale_shape_manual(values=c(16,1))+xlab("rate of evolution")+ylab("")

dev.off()

PF_ex<-Data[,2];names(PF_ex)<-rownames(Data)
PF_ex<-PF_ex[complete.cases(PF_ex)]

PF.cm<-contMap(presTree,PF_ex,outline=TRUE,lwd=6)
PF.cm<-setMap(PF.cm,colors=c("white","black"))

pdf("~/Dropbox/Warbler.Molt.Migration/PFextentTree.pdf")
plotTree.wBars(PF.cm$tree,PF_ex,method="plotSimmap",colors=PF.cm$cols,scale=.05)
dev.off()

pdf("~/Dropbox/Warbler.Molt.Migration/PFextentplot.pdf")
plot(PF.cm,fsize=0)
dev.off()







lamMat<-matrix(nrow=14,ncol=2)
rownames(lamMat)<-colnames(PF_tr)
colnames(lamMat)<-c("lambda","reg")

lamMat[1,1]<-phylosig(presTree,hPF)
lamMat[2,1]<-phylosig(presTree,brPF)
lamMat[3,1]<-phylosig(presTree,bePF)
lamMat[4,1]<-phylosig(presTree,baPF)
lamMat[5,1]<-phylosig(presTree,rePF)
lamMat[6,1]<-phylosig(presTree,prPF)
lamMat[7,1]<-phylosig(presTree,sePF)
lamMat[8,1]<-phylosig(presTree,trPF)
lamMat[9,1]<-phylosig(presTree,mcPF)
lamMat[10,1]<-phylosig(presTree,gcPF)
lamMat[11,1]<-phylosig(presTree,pcPF)
lamMat[12,1]<-phylosig(presTree,alPF)
lamMat[13,1]<-phylosig(presTree,PF_ex)
lamMat[14,1]<-phylosig(presTree,PF_pr)




lamMat[,2]<-rownames(lamMat)
lamMat[6,1]<-1
lamMat[11,1]<-1
lamMat[12,1]<-1

lamMat<-data.frame(lamMat)

lamMat$reg<-gsub("PF_","",lamMat$reg)
lamMat$reg<-gsub("rectrices","rr",lamMat$reg)
lamMat$reg<-gsub("primaries","pp",lamMat$reg)
lamMat$reg<-gsub("secondaries","ss",lamMat$reg)
lamMat$reg<-gsub("med_cov","m_cov",lamMat$reg)
lamMat$reg<-gsub("primary_coverts","p_cov",lamMat$reg)

lamMat$reg<-c("head","breast","belly","back","rr","pp","ss","tertials","m_cov","gr_cov","p_cov","alula","extent","presence")

write.csv(lamMat,"~/Dropbox/Warbler.Molt.Migration/PFlambda.csv")



pdf("~/Dropbox/Warbler.Molt.Migration/PF_phylsignal.pdf")
ggplot(lamMat,aes(y=as.numeric(as.character(lambda)),x=reg))+geom_point(size=4)+theme_tufte()+theme(panel.border=element_rect(fill="transparent"),axis.text.x = element_text(size=12,angle=45,vjust=.5),plot.margin=unit(c(1,1,12,1),"cm"))+scale_shape_manual(values=c(16,1))+xlab("phylogenetic signal")+ylab("")
dev.off()



#################################################
###### Figure for TD by feather region ############
#################################################


TD_tr<-read.csv("/Users/ryanterrill/Dropbox/Warbler.Molt.Migration/Parulidae_temporal_dimorphism.csv",row.names=1)

TD_or<-data.frame(1:51);rownames(TD_or)<-presTree$tip.label

TD_tr<-merge(TD_or,TD_tr,by=0)
TD_tr<-TD_tr[order(TD_tr$X1.51),]
rownames(TD_tr)<-TD_tr$Row.names;TD_tr<-TD_tr[,3:14]

TD_tr[TD_tr==.5]<-1
TD_tr[TD_tr==.6]<-1

hTD<-TD_tr[,1];names(hTD)<-rownames(TD_tr)
headTD<-ace(hTD,presTree,model="ER",type="discrete")

brTD<-TD_tr[,2];names(brTD)<-rownames(TD_tr)
breastTD<-ace(brTD,presTree,model="ER",type="discrete")

beTD<-TD_tr[,3];names(beTD)<-rownames(TD_tr)
bellyTD<-ace(beTD,presTree,model="ER",type="discrete")

baTD<-TD_tr[,4];names(baTD)<-rownames(TD_tr)
backTD<-ace(baTD,presTree,model="ER",type="discrete")

reTD<-TD_tr[,5];names(reTD)<-rownames(TD_tr)
### all absent


prTD<-TD_tr[,6];names(prTD)<-rownames(TD_tr)
### all absent

seTD<-TD_tr[,7];names(seTD)<-rownames(TD_tr)
### all absent


trTD<-TD_tr[,8];names(trTD)<-rownames(TD_tr)
tertTD<-ace(trTD,presTree,model="ER",type="discrete")

mcTD<-TD_tr[,9];names(trTD)<-rownames(TD_tr)
medcovTD<-ace(mcTD,presTree,model="ER",type="discrete")


gcTD<-TD_tr[,10];names(gcTD)<-rownames(TD_tr)
grcovTD<-ace(gcTD,presTree,model="ER",type="discrete")


pcTD<-TD_tr[,11];names(pcTD)<-rownames(TD_tr)
### all absent


alTD<-TD_tr[,12];names(alTD)<-rownames(TD_tr)
## all absent




#####Plot by feather region
pdf("~/Dropbox/Warbler.Molt.Migration/TD_ACE_byTract.pdf")

par(mar=c(1,4,3,10))


plot(presTree,cex=.5,main = "Seasonal Dichromatism by Feather Region",label.offset=.32, show.node.label=TRUE)
nodelabels(pie=headTD$lik.anc,piecol=cols,adj=c(0.48,0.5),cex=0.2)
nodelabels(pie=breastTD$lik.anc,piecol=cols,adj=c(0.49,0.5),cex=0.2)
nodelabels(pie=bellyTD$lik.anc,piecol=cols,adj=c(0.50,0.5),cex=0.2)
nodelabels(pie=backTD$lik.anc,piecol=cols,adj=c(0.51,0.5),cex=0.2)
nodelabels(pie=zeroMat,piecol=cols,adj=c(0.52,0.5),cex=0.2)
nodelabels(pie=zeroMat,piecol=cols,adj=c(0.53,0.5),cex=0.2)
nodelabels(pie=zeroMat,piecol=cols,adj=c(0.54,0.5),cex=0.2)
nodelabels(pie=tertTD$lik.anc,piecol=cols,adj=c(0.55,0.5),cex=0.2)
nodelabels(pie=medcovTD$lik.anc,piecol=cols,adj=c(0.56,0.5),cex=0.2)
nodelabels(pie=grcovTD$lik.anc,piecol=cols,adj=c(0.57,0.5),cex=0.2)
nodelabels(pie=zeroMat,piecol=cols,adj=c(0.58,0.5),cex=0.2)
nodelabels(pie=zeroMat,piecol=cols,adj=c(0.59,0.5),cex=0.2)




tiplabels(pie = to.matrix(hTD, sort(unique(hTD))), piecol = cols,adj=c(.58,0.5), cex=0.3)
tiplabels(pie = to.matrix(brTD, sort(unique(brTD))), piecol = cols,adj=c(.60,0.5), cex=0.3)
tiplabels(pie = to.matrix(beTD, sort(unique(beTD))), piecol = cols,adj=c(.62,0.5), cex=0.3)
tiplabels(pie = to.matrix(baTD, sort(unique(baTD))), piecol = cols,adj=c(.64,0.5), cex=0.3)
tiplabels(pie = to.matrix(reTD, sort(unique(reTD))), piecol = cols,adj=c(.66,0.5), cex=0.3)
tiplabels(pie = to.matrix(prTD, sort(unique(prTD))), piecol = cols,adj=c(.68,0.5), cex=0.3)
tiplabels(pie = to.matrix(seTD, sort(unique(seTD))), piecol = cols,adj=c(.70,0.5), cex=0.3)
tiplabels(pie = to.matrix(trTD, sort(unique(trTD))), piecol = cols,adj=c(.72,0.5), cex=0.3)
tiplabels(pie = to.matrix(mcTD, sort(unique(mcTD))), piecol = cols,adj=c(.74,0.5), cex=0.3)
tiplabels(pie = to.matrix(gcTD, sort(unique(gcTD))), piecol = cols,adj=c(.76,0.5), cex=0.3)
tiplabels(pie = to.matrix(pcTD, sort(unique(pcTD))), piecol = cols,adj=c(.78,0.5), cex=0.3)
tiplabels(pie = to.matrix(alTD, sort(unique(alTD))),piecol=cols,adj=c(.80,0.5),cex=0.3)

dev.off()


###### 
#plot mean and TD for transition rates for each region

rateMat<-matrix(nrow=12,ncol=3)
rownames(rateMat)<-colnames(TD_tr)
colnames(rateMat)<-c("rate","se","reg")


##set up matrix with rates and SEs
rateMat[1,1]<-headTD$rates
rateMat[1,2]<-headTD$se

rateMat[2,1]<-breastTD$rates[1]
rateMat[2,2]<-breastTD$se[1]

rateMat[3,1]<-bellyTD$rates[1]
rateMat[3,2]<-bellyTD$se[1]


rateMat[4,1]<-backTD$rates[1]
rateMat[4,2]<-backTD$se[1]


rateMat[5,1]<-rectTD$rates[1]
rateMat[5,2]<-rectTD$se[1]


rateMat[6,1]<-primTD$rates[1]
rateMat[6,2]<-primTD$se[1]

rateMat[7,1]<-secTD$rates[1]
rateMat[7,2]<-secTD$se[1]


rateMat[8,1]<-tertTD$rates[1]
rateMat[8,2]<-tertTD$se[1]

rateMat[9,]<-rep(as.numeric(0),3)


rateMat[10,1]<-grcovTD$rates[1]
rateMat[10,2]<-grcovTD$se[1]

rateMat[11,1]<-prcovTD$rates[1]
rateMat[11,2]<-prcovTD$se[1]

rateMat[12,1]<-alulTD$rates[1]
rateMat[12,2]<-alulTD$se[1]



rateMat[,3]<-rownames(rateMat)

rateMat<-data.frame(rateMat)

####need to change matrix so mean and TD are columns


rateMat$reg<-gsub("TD_","",rateMat$reg)
rateMat$reg<-gsub("rectrices","rr",rateMat$reg)
rateMat$reg<-gsub("primaries","pp",rateMat$reg)
rateMat$reg<-gsub("secondaries","ss",rateMat$reg)
rateMat$reg<-gsub("med_cov","m_cov",rateMat$reg)
rateMat$reg<-gsub("primary_coverts","p_cov",rateMat$reg)

rateMat$reg<-factor(rateMat$reg, level = c("head","breast","belly","back","rr","pp","ss","tertials","m_cov","gr_cov","p_cov","alula"))

limits <- aes(ymax = as.numeric(as.character(rate)) + as.numeric(as.character(se)), ymin=as.numeric(as.character(rate)) - as.numeric(as.character(se)))

write.csv(t(rateMat),"~/Dropbox/Warbler.Molt.Migration/TD_tractRates.csv")




pdf("~/Dropbox/Warbler.Molt.Migration/TD_tracttransitions.pdf")

ggplot(rateMat,aes(y=as.numeric(as.character(rate)),x=reg))+geom_errorbar(limits)+geom_point(size=6)+theme_tufte()+theme(panel.border=element_rect(fill="transparent"),axis.text.x = element_text(size=12,angle=45,vjust=.5),plot.margin=unit(c(1,1,12,1),"cm"))+scale_shape_manual(values=c(16,1))+xlab("rate of evolution")+ylab("")

dev.off()

TD_ex<-Data[,3];names(TD_ex)<-rownames(Data)
TD_ex<-TD_ex[complete.cases(TD_ex)]

TD.cm<-contMap(presTree,TD_ex,outline=TRUE,lwd=6)
TD.cm<-setMap(TD.cm,colors=c("white","black"))

pdf("~/Dropbox/Warbler.Molt.Migration/TDextentTree.pdf")
plotTree.wBars(TD.cm$tree,TD_ex,method="plotSimmap",colors=TD.cm$cols,scale=.05)
dev.off()

pdf("~/Dropbox/Warbler.Molt.Migration/TDextentplot.pdf")
plot(TD.cm,fsize=0)
dev.off()







lamMat<-matrix(nrow=12,ncol=2)
rownames(lamMat)<-colnames(TD_tr)
colnames(lamMat)<-c("lambda","reg")

lamMat[1,1]<-phylosig(presTree,hTD)
lamMat[2,1]<-phylosig(presTree,brTD)
lamMat[3,1]<-phylosig(presTree,beTD)
lamMat[4,1]<-phylosig(presTree,baTD)
lamMat[5,1]<-phylosig(presTree,reTD)
lamMat[6,1]<-phylosig(presTree,prTD)
lamMat[7,1]<-phylosig(presTree,seTD)
lamMat[8,1]<-phylosig(presTree,trTD)
lamMat[9,1]<-phylosig(presTree,mcTD)
lamMat[10,1]<-phylosig(presTree,gcTD)
lamMat[11,1]<-phylosig(presTree,pcTD)
lamMat[12,1]<-phylosig(presTree,alTD)
lamMat[13,1]<-phylosig(presTree,TD_ex)






lamMat[,2]<-rownames(lamMat)
lamMat[6,1]<-1
lamMat[11,1]<-1
lamMat[12,1]<-1

lamMat<-data.frame(lamMat)

lamMat$reg<-gsub("TD_","",lamMat$reg)
lamMat$reg<-gsub("rectrices","rr",lamMat$reg)
lamMat$reg<-gsub("primaries","pp",lamMat$reg)
lamMat$reg<-gsub("secondaries","ss",lamMat$reg)
lamMat$reg<-gsub("med_cov","m_cov",lamMat$reg)
lamMat$reg<-gsub("primary_coverts","p_cov",lamMat$reg)

lamMat$reg<-factor(lamMat$reg, level = c("head","breast","belly","back","rr","pp","ss","tertials","m_cov","gr_cov","p_cov","alula"))



pdf("~/Dropbox/Warbler.Molt.Migration/TD_phylsignal.pdf")
ggplot(lamMat,aes(y=as.numeric(as.character(lambda)),x=reg))+geom_point(size=4)+theme_tufte()+theme(panel.border=element_rect(fill="transparent"),axis.text.x = element_text(size=12,angle=45,vjust=.5),plot.margin=unit(c(1,1,12,1),"cm"))+scale_shape_manual(values=c(16,1))+xlab("phylogenetic signal")+ylab("")
dev.off()


#################################################
###### Figure for SD by feather region ############
#################################################


SD_tr<-read.csv("/Users/ryanterrill/Dropbox/Warbler.Molt.Migration/Parulidae_sexual_dimorphism.csv",row.names=1)

SD_or<-data.frame(1:51);rownames(SD_or)<-presTree$tip.label

SD_tr<-merge(SD_or,SD_tr,by=0)
SD_tr<-SD_tr[order(SD_tr$X1.51),]
rownames(SD_tr)<-SD_tr$Row.names;SD_tr<-SD_tr[,3:14]

SD_tr[SD_tr==.5]<-1
SD_tr[SD_tr==.6]<-1

hSD<-SD_tr[,1];names(hSD)<-rownames(SD_tr)
headSD<-ace(hSD,presTree,model="ER",type="discrete")

brSD<-SD_tr[,2];names(brSD)<-rownames(SD_tr)
breastSD<-ace(brSD,presTree,model="ER",type="discrete")

beSD<-SD_tr[,3];names(beSD)<-rownames(SD_tr)
bellySD<-ace(beSD,presTree,model="ER",type="discrete")

baSD<-SD_tr[,4];names(baSD)<-rownames(SD_tr)
backSD<-ace(baSD,presTree,model="ER",type="discrete")

reSD<-SD_tr[,5];names(reSD)<-rownames(SD_tr)
### all absent


prSD<-SD_tr[,6];names(prSD)<-rownames(SD_tr)
### all absent

seSD<-SD_tr[,7];names(seSD)<-rownames(SD_tr)
### all absent


trSD<-SD_tr[,8];names(trSD)<-rownames(SD_tr)
tertSD<-ace(trSD,presTree,model="ER",type="discrete")

mcSD<-SD_tr[,9];names(trSD)<-rownames(SD_tr)
medcovSD<-ace(mcSD,presTree,model="ER",type="discrete")


gcSD<-SD_tr[,10];names(gcSD)<-rownames(SD_tr)
grcovSD<-ace(gcSD,presTree,model="ER",type="discrete")


pcSD<-SD_tr[,11];names(pcSD)<-rownames(SD_tr)
### all absent


alSD<-SD_tr[,12];names(alSD)<-rownames(SD_tr)
## all absent




#####Plot by feather region
pdf("~/Dropbox/Warbler.Molt.Migration/SD_ACE_byTract.pdf")

par(mar=c(1,4,3,10))


plot(presTree,cex=.5,main = "Seasonal Dichromatism by Feather Region",label.offset=.32, show.node.label=TRUE)
nodelabels(pie=headSD$lik.anc,piecol=cols,adj=c(0.48,0.5),cex=0.2)
nodelabels(pie=breastSD$lik.anc,piecol=cols,adj=c(0.49,0.5),cex=0.2)
nodelabels(pie=bellySD$lik.anc,piecol=cols,adj=c(0.50,0.5),cex=0.2)
nodelabels(pie=backSD$lik.anc,piecol=cols,adj=c(0.51,0.5),cex=0.2)
nodelabels(pie=zeroMat,piecol=cols,adj=c(0.52,0.5),cex=0.2)
nodelabels(pie=zeroMat,piecol=cols,adj=c(0.53,0.5),cex=0.2)
nodelabels(pie=zeroMat,piecol=cols,adj=c(0.54,0.5),cex=0.2)
nodelabels(pie=tertSD$lik.anc,piecol=cols,adj=c(0.55,0.5),cex=0.2)
nodelabels(pie=medcovSD$lik.anc,piecol=cols,adj=c(0.56,0.5),cex=0.2)
nodelabels(pie=grcovSD$lik.anc,piecol=cols,adj=c(0.57,0.5),cex=0.2)
nodelabels(pie=zeroMat,piecol=cols,adj=c(0.58,0.5),cex=0.2)
nodelabels(pie=zeroMat,piecol=cols,adj=c(0.59,0.5),cex=0.2)




tiplabels(pie = to.matrix(hSD, sort(unique(hSD))), piecol = cols,adj=c(.58,0.5), cex=0.3)
tiplabels(pie = to.matrix(brSD, sort(unique(brSD))), piecol = cols,adj=c(.60,0.5), cex=0.3)
tiplabels(pie = to.matrix(beSD, sort(unique(beSD))), piecol = cols,adj=c(.62,0.5), cex=0.3)
tiplabels(pie = to.matrix(baSD, sort(unique(baSD))), piecol = cols,adj=c(.64,0.5), cex=0.3)
tiplabels(pie = to.matrix(reSD, sort(unique(reSD))), piecol = cols,adj=c(.66,0.5), cex=0.3)
tiplabels(pie = to.matrix(prSD, sort(unique(prSD))), piecol = cols,adj=c(.68,0.5), cex=0.3)
tiplabels(pie = to.matrix(seSD, sort(unique(seSD))), piecol = cols,adj=c(.70,0.5), cex=0.3)
tiplabels(pie = to.matrix(trSD, sort(unique(trSD))), piecol = cols,adj=c(.72,0.5), cex=0.3)
tiplabels(pie = to.matrix(mcSD, sort(unique(mcSD))), piecol = cols,adj=c(.74,0.5), cex=0.3)
tiplabels(pie = to.matrix(gcSD, sort(unique(gcSD))), piecol = cols,adj=c(.76,0.5), cex=0.3)
tiplabels(pie = to.matrix(pcSD, sort(unique(pcSD))), piecol = cols,adj=c(.78,0.5), cex=0.3)
tiplabels(pie = to.matrix(alSD, sort(unique(alSD))),piecol=cols,adj=c(.80,0.5),cex=0.3)

dev.off()


###### 
#plot mean and SD for transition rates for each region

rateMat<-matrix(nrow=12,ncol=3)
rownames(rateMat)<-colnames(SD_tr)
colnames(rateMat)<-c("rate","se","reg")


##set up matrix with rates and SEs
rateMat[1,1]<-headSD$rates
rateMat[1,2]<-headSD$se

rateMat[2,1]<-breastSD$rates[1]
rateMat[2,2]<-breastSD$se[1]

rateMat[3,1]<-bellySD$rates[1]
rateMat[3,2]<-bellySD$se[1]


rateMat[4,1]<-backSD$rates[1]
rateMat[4,2]<-backSD$se[1]


rateMat[5,1]<-rectSD$rates[1]
rateMat[5,2]<-rectSD$se[1]


rateMat[6,1]<-primSD$rates[1]
rateMat[6,2]<-primSD$se[1]

rateMat[7,1]<-secSD$rates[1]
rateMat[7,2]<-secSD$se[1]


rateMat[8,1]<-tertSD$rates[1]
rateMat[8,2]<-tertSD$se[1]

rateMat[9,]<-rep(as.numeric(0),3)


rateMat[10,1]<-grcovSD$rates[1]
rateMat[10,2]<-grcovSD$se[1]

rateMat[11,1]<-prcovSD$rates[1]
rateMat[11,2]<-prcovSD$se[1]

rateMat[12,1]<-alulSD$rates[1]
rateMat[12,2]<-alulSD$se[1]



rateMat[,3]<-rownames(rateMat)

rateMat<-data.frame(rateMat)

####need to change matrix so mean and SD are columns


rateMat$reg<-gsub("SD_","",rateMat$reg)
rateMat$reg<-gsub("rectrices","rr",rateMat$reg)
rateMat$reg<-gsub("primaries","pp",rateMat$reg)
rateMat$reg<-gsub("secondaries","ss",rateMat$reg)
rateMat$reg<-gsub("med_cov","m_cov",rateMat$reg)
rateMat$reg<-gsub("primary_coverts","p_cov",rateMat$reg)

rateMat$reg<-factor(rateMat$reg, level = c("head","breast","belly","back","rr","pp","ss","tertials","m_cov","gr_cov","p_cov","alula"))


write.csv(t(rateMat),"~/Dropbox/Warbler.Molt.Migration/SD_tractRates.csv")








modelACE<-function(matrix,phy,MODEL){
	rpars<-list()
		tree.dat<-treedata(phy,na.omit(matrix))
		dat<-as.vector(tree.dat$data);names(dat)<-rownames(tree.dat$data)
			rpars$sig2<-anc.ML(tree.dat$phy,dat,model=MODEL)$sig2
			rpars$ace<-anc.ML(tree.dat$phy,dat,model=MODEL)$ace
			return(rpars)
		}

PArates<-modelACE(PA[,1,drop=FALSE],warTree,MODEL="OU")
SDrates<-modelACE(PF[,1,drop=FALSE],warTree,MODEL="EB")
SDrates<-modelACE(SD[,1,drop=FALSE],warTree,MODEL="OU")
TDrates<-modelACE(TD[,1,drop=FALSE],warTree,MODEL="BM")

### plot rates by model?

#estimate gains and losses with stochastic character maps
#presence of molt and dichr in whole body and by tract




simmap<-function(matrix,phy,MODEL){
		tree.dat<-treedata(phy,na.omit(matrix))
		dat<-as.vector(tree.dat$data);names(dat)<-rownames(tree.dat$data)
			sim<-make.simmap(tree.dat$phy,dat,model=MODEL,nsim=1000)
			return(sim)
		}

paSimExtent<-simmap(PA[,1,drop=FALSE],warTree,MODEL="ARD")
PAs<-describe.simmap(paSimExtent,plot=FALSE)
pfSimExtent<-simmap(PF[,1,drop=FALSE],warTree,MODEL="ARD")
PFs<-describe.simmap(pfSimExtent,plot=FALSE)
sdSimExtent<-simmap(SD[,1,drop=FALSE],warTree,MODEL="ER")
SDs<-describe.simmap(sdSimExtent,plot=FALSE)
tdSimExtent<-simmap(TD[,1,drop=FALSE],warTree,MODEL="ARD")
TDs<-describe.simmap(tdSimExtent,plot=FALSE)

##save the simmaps
##plot the simmaps


colors<-setNames(c("lightblue","red"),c("0","1"))
plot(PAs,colors=colors,fsize=.8)
plot(SDs,colors=colors,cex=.3,fsize=.6)
plot(TDs,colors=colors,cex=.3,fsize=.6)






TD_ex<-Data[,3];names(TD_ex)<-rownames(Data)
TD_ex<-TD_ex[complete.cases(TD_ex)]

TD.cm<-contMap(presTree,TD_ex,outline=TRUE,lwd=6)
TD.cm<-setMap(TD.cm,colors=c("white","black"))

pdf("~/Dropbox/Warbler.Molt.Migration/TDextentTree.pdf")
plotTree.wBars(TD.cm$tree,TD_ex,method="plotSimmap",colors=TD.cm$cols,scale=.05)
dev.off()

pdf("~/Dropbox/Warbler.Molt.Migration/TDextentplot.pdf")
plot(PF.cm,fsize=0)
dev.off()





Mig_ex<-Data[,6];names(Mig_ex)<-rownames(Data)
Mig_ex<-Mig_ex[complete.cases(Mig_ex)]

Mig.cm<-contMap(presTree,Mig_ex,outline=TRUE,lwd=6)
Mig.cm<-setMap(Mig.cm,colors=c("white","black"))

pdf("~/Dropbox/Warbler.Molt.Migration/MigDistTree.pdf")
plotTree.wBars(Mig.cm$tree,Mig_ex,method="plotSimmap",colors=Mig.cm$cols,scale=.05)
dev.off()

pdf("~/Dropbox/Warbler.Molt.Migration/MigDistplot.pdf")
plot(PF.cm,fsize=0)
dev.off()




Strat_ex<-Data[,14];names(Strat_ex)<-rownames(Data)
Strat_ex<-Strat_ex[complete.cases(Strat_ex)]

Strat.cm<-contMap(presTree,Strat_ex,outline=TRUE,lwd=6)
Strat.cm<-setMap(Strat.cm,colors=c("white","black"))

pdf("~/Dropbox/Warbler.Molt.Migration/StratTree.pdf")
plotTree.wBars(Strat.cm$tree,Strat_ex,method="plotSimmap",colors=Strat.cm$cols,scale=.05)
dev.off()

pdf("~/Dropbox/Warbler.Molt.migration/Stratplot.pdf")
plot(PF.cm,fsize=0)
dev.off()

pdf("~/Dropbox/Warbler.Molt.migration/Mig~PA.pdf")
ggplot(Data,aes(x=migDist,y=paextent))+geom_jitter()+theme_tufte()+theme(panel.border=element_rect(fill="transparent"))+geom_smooth(method="lm",col="black")+xlim(0,50)
dev.off()


pdf("~/Dropbox/Warbler.Molt.migration/TD~Strat.pdf")
ggplot(Data,aes(x=winter_strat,y=tdExtent))+geom_jitter()+theme_tufte()+theme(panel.border=element_rect(fill="transparent"))+geom_smooth(method="lm",col="black")+xlim(0,4)
dev.off()











###write a function or loop that makes a simmaps by column names - and do this for body region & mod unit



###make simmaps for each tract, plot together somehow; maybe in a matrix w/ tip labels only on rigthmost tree


###f'n calculates simmaps for each state, outputs posterior distibutions of gains and losses




###f'n calculates ACE for continuous characters and saves trees with plotted states



###also maybe make a plot with all the trees? 



















