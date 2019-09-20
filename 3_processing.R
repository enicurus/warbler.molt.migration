library(ape)
library(caper)
library(nlme)
library(RColorBrewer)
source('~/Dropbox/Warbler.Molt.Migration/name.check.R', chdir = TRUE)

taxa.key=read.delim('~/Dropbox/Warbler.Molt.Migration/Parulidae_taxonomies.txt',stringsAsFactors=F)
phy=read.nexus('~/Dropbox/Warbler.Molt.Migration/Parulidae_phylogeny/Trees/r8stree.tre')
d = read.delim('~/Dropbox/Warbler.Molt.Migration/Parulidae_molt.txt',stringsAsFactors=F)
migdist = read.delim('~/Dropbox/Warbler.Molt.Migration/migratory_distance.txt')
solar = read.delim('~/Dropbox/Warbler.Molt.Migration/solar_daylight.txt')
migdist$Name = rownames(migdist)




#combine solar and migdist
head(solar)
head(migdist)
x = merge(migdist[,c('Shapefile','migdist1')],solar[,c('Shapefile','strategy','solar','daylength')],by='Shapefile')
rownames(x) = rownames(migdist)
migdist = x
head(migdist)


#edit the phylogney
outgroups = c('Granatellus_pelzelni','Coereba_flaveola','Icteria_virens','Microligea_palustris','Xenoligea_montana','Zeledonia_coronata','Spindalis_zena','Teretistris_fernandinae')
phy = drop.tip(phy,outgroups)

###clade ages
ages = data.frame(phy$tip.label,phy$edge.length[phy$edge[,2] %in% 1:length(phy$tip.label)])
colnames(ages) = c('AOU','age')


#check to make sure all the phy$tip.label's are in the data
name.check(phy$tip.label,d$Scientific.Name)
name.check(phy$tip.label, taxa.key$AOU)

#edit the data to only include taxa in the phylogeny
data = d[d$Scientific.Name %in% phy$tip.label,c("Scientific.Name","PA.score",'PF.score','Seasonal.score','Sex.score')]
names(data)[1] = 'AOU'
head(data)

#edit the migdist so that the species names are AOU not birdlife
migdist$Birdlife = rownames(migdist)
migdist = migdist[c('Birdlife','strategy','migdist1','solar','daylength')]
rownames(migdist) = seq(1:nrow(migdist))
name.check(migdist$Birdlife,taxa.key$Birdlife)
x = merge(taxa.key,migdist,all.y=F)
x = merge(x,ages,all.y=F)
migdist = x[order(x$AOU),c('AOU','strategy','migdist1','solar','daylength','age')]
migdist = migdist[-c(which(duplicated(migdist$AOU)),which(is.na(migdist$AOU))),]
name.check(migdist$AOU,taxa.key$AOU)
name.check(migdist$AOU,phy$tip.label)
name.check(migdist$AOU,data$AOU)



############
#Analysis
############
d = merge(data,migdist,all.y=T)
name.check(d$AOU,phy$tip.label)

#make comparative data object
comp = comparative.data(phy=phy,data=d,names.col='AOU',vcv=T)


## Add an alpha value to a colour
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
                     function(x) 
                       rgb(x[1], x[2], x[3], alpha=alpha))  
}

key = setNames(add.alpha(brewer.pal(4,'Set1'),alpha=0.75),c('m','nm','mix','par'))
bg = as.character(key[unlist(comp$data$strategy)])



fit = pgls(PA.score ~ daylength,data = comp, lambda='ML')
summary(fit)
plot(comp$data$daylength,jitter(comp$data$PA.score,.25),bg=bg,col='black',pch=21,cex=2,lwd=0.5)
#text(comp$data$daylength,comp$data$PA.score,comp$phy$tip.label)
abline(fit)
legend(13.7,.5, names(key),pch = 21,pt.cex=2,pt.=.5,pt.bg=key) # gives the legend lines the correct color and width




x = data.frame(comp$phy$tip.label,comp$data$strategy,comp$data$daylength,comp$data$PA.score)
x[order(x[,4],x[,3]),]


fit = pgls(migdist1 ~ daylength,data = comp, lambda='ML')
summary(fit)
plot(comp$data$daylength,comp$data$migdist1,bg=bg,col='black',pch=21,cex=2,lwd=0.5)
#text(comp$data$daylength,comp$data$migdist1,comp$phy$tip.label)
abline(fit)

fit = pgls(PA.score ~ solar,data = comp, lambda='ML')
summary(fit)
plot(comp$data$solar,jitter(comp$data$PA.score,.25),bg=bg,col='black',pch=21,cex=2,lwd=0.5)
abline(fit)

fit = pgls(migdist1 ~ solar,data = comp, lambda='ML')
summary(fit)
plot(comp$data$solar,comp$data$migdist1,bg=bg,col='black',pch=21,cex=2,lwd=0.5)
abline(fit)

fit = pgls(PA.score~migdist1,data = comp, lambda='ML')
summary(fit)
plot(comp$data$migdist1,jitter(comp$data$PA.score,0.25),bg=bg,col='black',pch=21,cex=2,lwd=0.5)
abline(fit)

fit = pgls(Seasonal.score~migdist1,data = comp, lambda='ML')
summary(fit)
plot(comp$data$migdist1,comp$data$Seasonal.score,bg=bg,col='black',pch=21,cex=2,lwd=0.5)
abline(fit)


fit = pgls(Seasonal.score~PA.score,data = comp, lambda='ML')
summary(fit)
plot(jitter(comp$data$PA.score,0.25),jitter(comp$data$Seasonal.score,0.25),bg=bg,col='black',pch=21,cex=2,lwd=0.5)
abline(fit)

fit = pgls(Seasonal.score~Sex.score,data = comp, lambda='ML')
summary(fit)
plot(comp$data$PA.score,comp$data$Seasonal.score,,bg=bg,col='black',pch=21,cex=2,lwd=0.5)
abline(fit)


head(comp$data)





d = merge(data,migdist,all.y=T)
name.check(d$AOU,phy$tip.label)

#make comparative data object
migs = migdist[migdist$strategy == 'm',] 

phy.migs = drop.tip(phy,setdiff(phy$tip.label,migs$AOU))



comp = comparative.data(phy=phy.migs,data=migs,names.col='AOU',vcv=T)



fit = pgls(migdist1 ~ age,data = comp, lambda='ML')
summary(fit)
plot(comp$data$age,jitter(comp$data$migdist1,.25),pch=21,cex=2,lwd=0.5)
#text(comp$data$daylength,comp$data$PA.score,comp$phy$tip.label)
abline(fit)


