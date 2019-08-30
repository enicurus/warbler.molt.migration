`#testing correlated evolution and transition rates of precense and absence of PA and SD


library(phytools)
library(geiger)

tree<-read.nexus("/Users/rterrill/Dropbox/Warbler.Molt.Migration/Parulidae_phylogeny/Trees/r8stree.tre")


Data<-read.csv("~/Dropbox/Warbler.Molt.Migration/Data_treeNames.csv",row.names=1)

matrix<-data.frame(Data$paextent,Data$migDist)

rownames(matrix)<-rownames(Data)

matrix$PA<-0
matrix$mig<-0


matrix$PA<-ifelse(matrix[,1]==0,0,1)
matrix$PA<-ifelse(is.na(matrix[,1]),0,matrix$PA)



matrix$mig<-ifelse(abs(matrix[,2])>20,1,0)
matrix$mig<-ifelse(is.na(matrix[,2]),0,matrix$mig)


PA<-matrix$PA
names(PA)<-rownames(matrix)

mig<-matrix$mig
names(mig)<-rownames(matrix)

drop<-name.check(tree,PA)$tree_not_data

tree<-drop.tip(tree,drop)



mK<-fitPagel(tree,PA,mig)
mK.x<-fitPagel(tree,PA,mig,drp.var="x")
ace<-fitPagel(tree,PA,mig,method="ace")
ace.x<-fitPagel(tree,PA,mig,method="ace",dep.var="x")
fD<-fitPagel(tree,PA,mig,method="fitDiscrete")
fD.x<-fitPagel(tree,PA,mig,method="fitDiscrete",dep.var="x")


aic<-setNames(c(mK$independent.AIC,mK$dependent.AIC,mK.x$dependent.AIC,
	ace$independent.AIC,ace$dependent.AIC,ace.x$dependent.AIC,fD$independent.AIC,
	fD$dependent.AIC,fD.x$dependent.AIC),c("mK independent","mK dependent","mK x dependent",
	"ace independent", "ace dependent", "ace x dependent", "fitDiscrete independent", "fitDiscrete dependent", 
	"fitDiscrete x dependent"))

aic.w(aic)

plot(ace.x,lwd.by.weight=TRUE)
