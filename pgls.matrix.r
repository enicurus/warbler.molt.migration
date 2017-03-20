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