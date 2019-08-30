#testing correlated evolution and transition rates of precense and absence of PA and SD


library(phytools)

tree<-read.nexus("/Users/rterrill/Dropbox/Warbler.Molt.Migration/Parulidae_phylogeny/Trees/r8stree.tre")


Data<-read.csv("~/Dropbox/Warbler.Molt.Migration/Data_treeNames.csv",row.names=1)

matrix<-data.frame(Data$paextent,Data$tdExtent)

rownames(matrix)<-rownames(Data)

matrix$PA<-0
matrix$SD<-0



matrix$PA<-ifelse(is.na(matrix[,1]),NA,0)
matrix$PA<-ifelse(matrix[,1]==0,0,1)


matrix$SD<-ifelse(is.na(matrix[,2]),NA,0)
matrix$SD<-ifelse(matrix[,2]==0,0,1)

matrix<-na.omit(matrix)

PA<-matrix$PA
names(PA)<-rownames(matrix)

SD<-matrix$SD
names(SD)<-rownames(matrix)

drop<-name.check(tree,PA)$tree_not_data

tree<-drop.tip(tree,drop)



fitPagel(tree,PA,SD)

fitPagel(tree,PA,SD,method="ace")

fitPagel(tree,PA,SD,method="fitDiscrete")