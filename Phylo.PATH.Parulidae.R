#### Phylogenetic PATH analysis
#### Follows approach of Ch 8 in Modern Phylogenetic Comparative Methods and Their Application in Evolutionary Biology




require(ape)
require(caper)
require(nlme)
require(geiger)


###setting up the data Matrix###



Data<-read.csv("~/Dropbox/Warbler.Molt.Migration/Data_treeNames.csv",row.names=1)
Data$names<-rownames(Data)

warTree<-read.nexus("/Users/ryanterrill/Dropbox/Warbler.Molt.Migration/Parulidae_phylogeny/Trees/r8stree.tre")
warTree<-multi2di(warTree)
warTree<-chronos(warTree)
class(warTree)<-"phylo"


path.dat<-Data[,c(1,2,3,5,6,14,17,31)]
path.dat<-path.dat[complete.cases(path.dat),]
path.phy<-drop.tip(warTree,setdiff(warTree$tip.label,rownames(path.dat)))
path.dat$sp<-rownames(path.dat)

colnames(path.dat)<-c("pa","pf","sd","dl","mig","strat","sol","mass","sp")

com.dat<-comparative.data(path.phy,path.dat,sp,vcv=TRUE,vcv.dim=3,warn.dropped=TRUE)

#####  Directed acyclic path graphs are represented in: ~/Dropbox/Warbler.Molt.Matrix/path_hypotheses.pdf  ####

##function to get number of conditional dependencies by number of vertices (V) and edges (A)
condNum <- function(V, A) {
    (factorial(V)/(2 * factorial(V - 2))) - A
}


#function to calculate CICc of models

CICc <- function(C, q, n) {
    C + 2 * q * (n/(n - 1 - q))
}



##Notes on translating directed acyclic graphs to conditional independencies: 
#first step is to use condNum to check how many CIs there should be
#second step is to list each unique set of non-adjacent vertices
#third step is to list parent variable of EACH non-adjecent variable - if variable has no parent don't list

##Note on translating CIs into Linear models
#start with the downstream non-adjacent vertex
# on the other side of the ~, list all other vertices in the CI, 
##Coefficients vertex is the upstream non-adjacent vertex


###These p-values are the probability of the model, not the probability of rejecting the model for the null!!!!!

#Thes CI statements are based on the PATH hypotheses figure outlining the acyclic graph hypotheses

#Model 1 (9 CIs)



#1.1 (dl,sd){mig,pa}
	m1.1<-pgls(sd~mig+pa+dl,data=com.dat,lambda="ML");m1.1p<-summary(m1.1)$coefficients["dl",4]
#1.2 (dl,pf){mig}
	m1.2<-pgls(pf~mig+dl,data=com.dat,lambda="ML");m1.2p<-summary(m1.2)$coefficients["dl",4]
#1.3 (dl,strat){mig}
	m1.3<-pgls(strat~mig+dl,data=com.dat,lambda="ML");m1.3p<-summary(m1.3)$coefficients["dl",4]
#1.4 (mig,pf){\(varnothing\)}
	m1.4<-pgls(pf~mig,data=com.dat,lambda="ML");m1.4p<-summary(m1.4)$coefficients["mig",4]
#1.5 (mig,sd){pa}
	m1.5<-pgls(sd~pa+mig,data=com.dat,lambda="ML");m1.5p<-summary(m1.5)$coefficients["mig",4]
#1.6 (mig,strat){\(varnothing\)}
	m1.6<-pgls(strat~mig,data=com.dat,lambda="ML");m1.6p<-summary(m1.6)$coefficients["mig",4]
#1.7 (pf,sd){pa}
	m1.7<-pgls(sd~pa+pf,data=com.dat,lambda="ML");m1.7p<-summary(m1.7)$coefficients["pf",4]
#1.8 (pf,strat){\(varnothing\)}
	m1.8<-pgls(strat~pf,data=com.dat,lambda="ML");m1.8p<-summary(m1.8)$coefficients["pf",4]
#1.9 (pa,strat){pf,dl,mig}
	m1.9<-pgls(strat~pf+dl+mig+pa,data=com.dat,lambda="ML");m1.9p<-summary(m1.9)$coefficients["pa",4]



C1<--2*(log(m1.1p)+log(m1.2p)+log(m1.3p)+log(m1.4p)+log(m1.5p)+log(m1.6p)+log(m1.7p)+log(m1.8p)+log(m1.9p))
C1.pval<-1-pchisq(C1,2*9)
CICc1<-CICc(C1,9,49)


#Model 2 (6 CIs) 

#2.1 (dl,sd){mig,pa} = 1.1
	m2.1<-m1.1;m2.1p<-m1.1p
#2.2 (dl,strat){mig} = 1.3
	m2.2<-m1.3;m2.2p<-m1.3p
#2.3 (mig,sd){pa} = 1.5
	m2.3<-m1.5;m2.3p<-m1.5p
#2.4 (mig,strat){\(varnothing\)} = 1.6
	m2.4<-m1.6;m2.4p<-m1.6p
#2.5 (pa,strat){dl,mig}
	m2.5<-pgls(strat~dl+mig+pa,data=com.dat,lambda="ML");m2.5p<-summary(m2.5)$coefficients["pa",4]
#2.6(mig,pa){dl}
	m2.6<-pgls(pa~dl+mig,data=com.dat,lambda="ML";m2.6p<-summary(m2.6)$coefficients["mig",4]

C2<--2*(log(m2.1p)+log(m2.2p)+log(m2.3p)+log(m2.4p)+log(m2.5p)+log(m2.6p))
C2.pval<-1-pchisq(C2,2*6)
CICc2<-CICc(C2,6,49)




#Model 3 (6 CIs)

#3.1 (dl,sd){mig,pa} = 1.1
	m3.1<-m1.1;m3.1p<-m1.1p
#3.2 (pa,strat){mig}
	m3.2<-pgls(strat~mig+pa,data=com.dat,lambda="ML"),m3.2p<-summary(m3.2)$coefficients["pa",4]
#3.3 (dl,strat){mig} = 1.3
	m3.3<-m1.3;m3.3p<-m1.3p
#3.4 (mig,sd){pa} = 1.5
	m3.4<-m1.5;m3.4p<-m1.5p
#3.5 (mig,strat){\(varnothing\)} = 1.6
	m3.5<-m1.6;m3.5p<-m1.6p
#3.6 (dl,pa){mig}
	m3.6<-pgls(pa~mig+dl,data=com.dat,lambda="ML")m3.6p<-summary(m3.6)$coefficiants["dl",4]


C3<--2*(log(m3.1p)+log(m3.2p)+log(m3.3p)+log(m3.4p)+log(m3.5p)+log(m3.6p))
C3.pval<-1-pchisq(C3,2*6)
CICc3<-CICc(C3,6,49)



#Model 4 (9 CIs)

#4.1 (dl,sd){mig,pa,strat}
	m4.1<-pgls(sd~mig+pa+strat+dl,data=com.dat,lambda="ML");m4.1p<-summary(m4.1)$coefficients["dl",4]
#4.2 (pa,strat){sol,dl,mig} = 1.9
	m4.2<-m1.9;m4.2p<-m1.9p
#4.3 (dl,strat){mig} = 1.3
	m4.3<-m1.3;m4.3p<-m1.3p
#4.4 (mig,sol){dl} = 2.4
	m4.4<-m2.4;m4.4p<-m2.4p
#4.5 (mig,sd){pa,strat}
	m4.5<-pgls(sd~pa+strat+mig,data=com.dat,lambda="ML");m4.5p<-summary(m4.5)$coefficients["mig",4]
#4.6 (sol,sd){pa,strat}
	m4.6<-pgls(sd~pa+dl+strat+sol,data=com.dat,,lambda="ML");m4.6p<-summary(m4.6)$coefficients["sol",4]
#4.7 (sol,strat){/(varnothing)}
	m4.7<-pgls(strat~dl+sol,data=com.dat,lambda="ML");m4.7p<-summary(m4.7)$coefficients["sol",4]
#4.8 (strat,mig){/(varnothing)}
	m4.8<-pgls(mig~strat,data=com.dat,lambda="ML");m4.8p<-summary(m4.8)$coefficients["strat",4]

C4<--2*(log(m4.1p)+log(m4.2p)+log(m4.3p)+log(m4.4p)+log(m4.5p)+log(m4.6p)+log(m4.7p)+log(m4.8p))
C4.pval<-1-pchisq(C4,2*8)
CICc4<-CICc(C4,8,49)


#Model 5 (5 CIs)

#5.1 (dl,sd){mig,pa} = 1.1
	m5.1<-m1.1;m5.1p<-m1.1p
#5.2 (dl,strat){mig} = 1.3
	m5.2<-m1.3;m5.2p<-m1.3p
#5.3 (mig,sd){pa} = 1.5
	m5.3<-m1.5;m5.3p<-m1.5p
#5.4 (mig,strat){\(varnothing\)} = 1.6
	m5.4<-m1.6;m5.4p<-m1.6p
#5.5 (pa,strat){mig,dl}
	m5.5<-pgls(strat~sol+dl+pa,data=com.dat,lambda="ML");m5.5p<-summary(m5.5)$coefficients["pa",4]


C5<--2*(log(m5.1p)+log(m5.2p)+log(m5.3p)+log(m5.4p)+log(m5.5p))
C5.pval<-1-pchisq(C5,2*9)
CICc5<-CICc(C5,9,49)



#Model 6 (10 CIs)


#6.1 (dl,sd){mig,pa} = 1.1
	m6.1<-m1.1;m6.1p<-m1.1p
#6.2 (dl,sol){mig} = 1.2 
	m6.2<-m1.2;m6.2p<-m1.2p
#6.3 (dl,strat){mig} = 1.3
	m6.3<-m1.3;m6.3p<-m1.3p
#6.4 (mig,sol){\(varnothing\)} = 1.4
	m6.4<-m1.4;m6.4p<-m1.4p
#6.5 (mig,sd){pa} = 1.5
	m6.5<-m1.5;m6.5p<-m1.5p
#6.6 (mig,strat){\(varnothing\)} = 1.6
	m6.6<-m1.6;m6.6p<-m1.6p
#6.7 (sol,sd){pa} = 1.7
	m6.7<-m1.7;m6.7p<-m1.7p
#6.8 (sol,strat){\(varnothing\)} = 1.8
	m6.8<-m1.8;m6.8p<-m1.8p
#6.9 (pa,strat){sol,mig}
	m6.9<-pgls(strat~sol+mig+pa,data=com.dat,lambda="ML");m6.9p<-summary(m6.9)$coefficients["pa",4]
#6.10 (dl,pa){mig}
	m6.10<-pgls(pa~mig+dl,data=com.dat,lambda="ML");m6.10p<-summary(m6.10)$coefficients["dl",4]


C6<--2*(log(m6.1p)+log(m6.2p)+log(m6.3p)+log(m6.4p)+log(m6.5p)+log(m6.6p)+log(m6.7p)+log(m6.8p)+log(m6.9p)+log(m6.10p))
C6.pval<-1-pchisq(C6,2*10)
CICc6<-CICc(C6,10,49)


#Model 7 (9 CIs)



#7.1 (dl,sd){mig,pa}
	m7.1<-pgls(sd~mig+pa+dl,data=com.dat,lambda="ML");m7.1p<-summary(m7.1)$coefficients["dl",4]
#7.2 (dl,mass){mig}
	m7.2<-pgls(mass~mig+dl,data=com.dat,lambda="ML");m7.2p<-summary(m7.2)$coefficients["dl",4]
#7.3 (dl,strat){mig}
	m7.3<-pgls(strat~mig+dl,data=com.dat,lambda="ML");m7.3p<-summary(m7.3)$coefficients["dl",4]
#7.4 (mig,mass){\(varnothing\)}
	m7.4<-pgls(mass~mig,data=com.dat,lambda="ML");m7.4p<-summary(m7.4)$coefficients["mig",4]
#7.5 (mig,sd){pa}
	m7.5<-pgls(sd~pa+mig,data=com.dat,lambda="ML");m7.5p<-summary(m7.5)$coefficients["mig",4]
#7.6 (mig,strat){\(varnothing\)}
	m7.6<-pgls(strat~mig,data=com.dat,lambda="ML");m7.6p<-summary(m7.6)$coefficients["mig",4]
#7.7 (mass,sd){pa}
	m7.7<-pgls(sd~pa+mass,data=com.dat,lambda="ML");m7.7p<-summary(m7.7)$coefficients["mass",4]
#7.8 (mass,strat){\(varnothing\)}
	m7.8<-pgls(strat~mass,data=com.dat,lambda="ML");m7.8p<-summary(m7.8)$coefficients["mass",4]
#7.9 (pa,strat){mass,dl,mig}
	m7.9<-pgls(strat~mass+dl+mig+pa,data=com.dat,lambda="ML");m7.9p<-summary(m7.9)$coefficients["pa",4]



C7<--2*(log(m7.1p)+log(m7.2p)+log(m7.3p)+log(m7.4p)+log(m7.5p)+log(m7.6p)+log(m7.7p)+log(m7.8p)+log(m7.9p))
C7.pval<-1-pchisq(C7,2*9)
CICc7<-CICc(C7,9,49)



#Model 8 (6 CIs)
#8.1 (dl,sd){mig,pa} = 1.1
	m8.1<-m1.1;m8.1p<-m1.1p
#8.2 (dl,strat){mig} = 1.3
	m8.2<-m1.3;m8.2p<-m1.3p
#8.3 (mig,pa),{/(varnothing)}
	m8.3<-pgls(pa~mig,data=com.dat,lambda="ML");m8.3p<-summary(m8.3)$coefficients["mig",4]
#8.4 (strat,sd){pa,mig}
	m8.4<-pgls(sd~pa+mig+strat,data=com.dat,lambda-"ML");m8.4p<-summary(m8.4)$coefficients["strat",4]
#8.5 (strat,mig){/(varnothing)}
	m8.5<-pgls(mig~strat,data=com.dat,lambda="ML");m8.5p<-summary(m8.5)$coefficients["strat",4]
#8.6 (dl,pa){strat}
	m8.5<-pgls(pa~strat,data=com.dat,lambda="ML")m8.6p<-summary(m8.6)$coefficients["dl",4]


C8<--2*(log(m8.1p)+log(m8.2p)+log(m8.3p)+log(m8.4p)+log(m8.5p)+log(m8.6p))
C8.pval<-1-pchisq(C8,2*6)
CICc8<-CICc(C8,6,49)






#Model 9 (9 CIs)
#9.1 (dl,sd){mig,pa} = 1.1
	m9.1<-m1.1;m9.1p<-m1.1p
#9.2 (dl,strat){mig} = 1.3
	m9.2<-m1.3;m9.2p<-m1.3p
#9.3 (strat,mass){dl}
	m9.3<-pgls(mass~dl+strat,data=com.dat,lambda="ML");m9.3p<-summary(m9.3)$coefficients["strat",4]
#9.4 (strat,pa){dl}
	m9.4<-pgls(pa~dl+strat,data=com.dat,lambda="ML");m9.4p<-summary(m9.4)$coefficients["strat",4]
#9.5 (strat,sd){pa,mig}
	m9.5<-pgls(sd~pa+strat+mig,data=com.dat,lambda="ML");m9.5p<-summary(m9.5)$coefficients["strat",4]
#9.6 (mig,mass){dl,strat}
	m9.6<-pgls(mass~dl+strat+mig,data=com.dat,lambda="ML");m9.6p<-summary(m9.6)$coefficients["mig",4]
#9.7 (mig,pa){strat,dl}
	m9.7<-pgls(pa~strat+dl+mig,data=com.dat,lambda="ML");m9.7p<-summary(m9.7)$coefficients["mig",4]
#9.8 (mass,pa){dl}
	m9.9<-pgls(pa~dl+mass,data=com.dat,lambda="ML");m9.8p<-summary(m9.8)$coefficients["mass",4]
#9.9 (mass,sd){dl,pa,mig}
	m9.9<-pgls(dl~pa+sd+mass+mig,data=com.dat,lambda="ML");m9.9p<-summary(m9.9)$coefficients["sol",4]
	

C9<--2*(log(m9.1p)+log(m9.2p)+log(m9.3p)+log(m9.4p)+log(m9.5p)+log(m9.6p)+log(m9.7p)+log(m9.8p)+log(m9.9p))
C9.pval<-1-pchisq(C9,2*10)
CICc9<-CICc(C9,10,49)







#Model 10 (6 CIs)
#10.1 (dl,pa){mig} = 6.10
	m10.1<-m6.10;m10.1p<-m6.10p
#10.2 (dl,strat){mig}
	m10.2<-pgls(strat~mig+dl,data=com.dat,lambda="ML");m10.2p<-summary(m10.2)$coefficients["dl",4]
#10.3 (mig,sd){pa,dl}
	m10.3<-pgls(sd~pa+dl+mig,data=com.dat,lambda="ML");m10.3p<-summary(m10.3)$coefficients["mig",4]
#10.4 (mig,strat){(/varnothing/}
	m10.4<-pgls(strat~mig,data=com.dat,lambda="ML");m10.4p<-summary(m10.4)$coefficients["mig",4]
#10.5 (mig,pa){strat}
	m10.5<-pgls(pa~strat+mig,data=com.dat,lambda="ML");m10.5p<-summary(m10.5)$coefficients["mig",4]
#10.6 (strat,sd){pa,dl}
	m10.6<-pgls(sd~pa+strat+dl,data=com.dat,lambda="ML");m10.6p<-summary(m10.6)$coefficients["strat",4]
	


C10<--2*(log(m10.1p)+log(m10.2p)+log(m10.3p)+log(m10.4p)+log(m10.5p)+log(m10.6p))
C10.pval<-1-pchisq(C10,2*6)
CICc10<-CICc(C10,6,49)







#Model 11 (6 CIs)
#11.1 (mig,sd){pa,strat} = 4.5
	m11.1<-m4.5;m11.1p<-m4.5p
#11.2 (strat,dl){mig} = 8.6
	m11.2<-m8.6;m11.2p<-m8.6p
#11.3 (strat,pa){(/varnothing)}
	m11.3<-pgls(pa~strat,data=com.dat,lambda="ML");m11.3p<-summary(m11.5)$coefficients["strat",4]
#11.4 (strat,sd){pa} = 9.5
	m11.4<-m9.5;m11.4p<-m9.5p
#11.5 (dl,pa){mig}
	m11.5<-pgls(mig~sol+pa+dl,data=com.dat,lambda="ML");m11.5p<-summary(m11.7)$coefficients["dl",4]
#11.6 (dl,sd){mig,pa} = 1.1
	m11.6<-m1.1;m11.6p<-m1.1p

	
C11<--2*(log(m11.1p)+log(m11.2p)+log(m11.3p)+log(m11.4p)+log(m11.5p)+log(m11.6p)+log(m11.7p)+log(m11.8p)+log(m11.9p))
C11.pval<-1-pchisq(C11,2*6)
CICc11<-CICc(C11,6,49)

	
 
#Model 12 (9 CIs)
#12.1 (mig,sol){dl} = 2.4
	m12.1<-m2.4;m12.1p<-m2.4p
#12.2 (mig,pa){dl,sol}
	m12.2<-pgls(pa~dl+sol+mig,data=com.dat,lambda="ML");m12.2p<-summary(m12.2)$coefficients["mig",4]
#12.3 (mig,sd){pa} = 1.5
	m12.3<-m1.5;m12.3p<-m1.5p
#12.4 (mig,strat){\(varnothing\)} = 1.6
	m12.4<-m1.6;m12.4p<-m1.6p
#12.5 (dl,sd){mig,pa,strat} = 4.1
	m12.5<-m4.1;m12.5p<-m4.1p
#12.6 (dl,strat){mig} = 1.3
	m12.6<-m1.3;m12.6p<-m1.3p
#12.7 (sol,sd){dl,pa} = 9.10
	m12.7<-m9.10;m12.7p<-m9.10p
#12.8 (sol,strat){dl} = 2.8
	m12.8<-m2.8;m12.8p<-m2.8p
#12.9 (pa,strat){sol,dl} = 5.9
	m12.9<-m5.9;m12.9p<-m5.9p

C12<--2*(log(m12.1p)+log(m12.2p)+log(m12.3p)+log(m12.4p)+log(m12.5p)+log(m12.6p)+log(m12.7p)+log(m12.8p)+log(m12.9p))
C12.pval<-1-pchisq(C12,2*9)
CICc12<-CICc(C12,9,49)



###models 4, 8 and 10 are overweighted because few CI's - maybe change or increase CIs ####


model.all <- c("Model 1", "Model 2", "Model 3", "Model 4", "Model 5", "Model 6", "Model 7", "Model 8", "Model 9","Model 10","Model 11","Model 12")
    
C.all <- c(C1, C2, C3, C4, C5, C6, C7, C8, C9,C10,C11,C12)
C.pval.all <- c(C1.pval, C2.pval, C3.pval, C4.pval, C5.pval, C6.pval, C7.pval, C8.pval, C9.pval,C10.pval,C11.pval,C12.pval)
CICc.all <- c(CICc1, CICc2, CICc3, CICc4, CICc5, CICc6, CICc7, CICc8, CICc9, CICc10, CICc11, CICc12)
C.all <- round(C.all, digits = 3)
C.pval.all <- signif(C.pval.all, digits = 3)
CICc.all <- round(CICc.all, digits = 3)
results.dat <- data.frame(model.all, C.all, C.pval.all, CICc.all)
names(results.dat) <- c("Model", "C statistic", "p-value", "CICc")
results.dat <- results.dat[order(results.dat$CICc), ]
print(results.dat) 

write.table(results.dat,file="~/Dropbox/Warbler.Molt.Migration/PATHResults.txt",quote=FALSE)