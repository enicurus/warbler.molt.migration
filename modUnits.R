####Calculating parameter values for modular units

mod<-read.csv("/Users/ryanterrill/Dropbox/Warbler.Molt.Migration/mod_units.csv")

PA_tr<-read.csv("/Users/ryanterrill/Dropbox/Warbler.Molt.Migration/Molt\ extent/Parulidae_PA.csv",row.names=1)
PF_tr<-read.csv("/Users/ryanterrill/Dropbox/Warbler.Molt.Migration/Molt\ extent/Parulidae_PF.csv",row.names=1)
SD_tr<-read.csv("~/Dropbox/Warbler.Molt.Migration/Parulidae_sexual_dimorphism.csv",row.names=1)
TD_tr<-read.csv("~/Dropbox/Warbler.Molt.Migration/Parulidae_temporal_dimorphism.csv",row.names=1)

colnames(PA_tr)<-sub("PA_","",names(PA_tr))
colnames(PF_tr)<-sub("PF_","",names(PF_tr))
colnames(SD_tr)<-sub("SD_","",names(SD_tr))
colnames(TD_tr)<-sub("TD_","",names(TD_tr))


PA1<-as.character(mod[mod$var=="PA"&mod$mod_unit==1,]$tract)
PA2<-as.character(mod[mod$var=="PA"&mod$mod_unit==2,]$tract)
PA3<-as.character(mod[mod$var=="PA"&mod$mod_unit==3,]$tract)
PA4<-as.character(mod[mod$var=="PA"&mod$mod_unit==4,]$tract)

PF1<-as.character(mod[mod$var=="PF"&mod$mod_unit==1,]$tract)
PF2<-as.character(mod[mod$var=="PF"&mod$mod_unit==2,]$tract)
PF3<-as.character(mod[mod$var=="PF"&mod$mod_unit==3,]$tract)
PF4<-as.character(mod[mod$var=="PF"&mod$mod_unit==4,]$tract)
PF5<-as.character(mod[mod$var=="PF"&mod$mod_unit==5,]$tract)
PF6<-as.character(mod[mod$var=="PF"&mod$mod_unit==6,]$tract)

SD1<-as.character(mod[mod$var=="SD"&mod$mod_unit==1,]$tract)
SD2<-as.character(mod[mod$var=="SD"&mod$mod_unit==2,]$tract)
SD3<-as.character(mod[mod$var=="SD"&mod$mod_unit==3,]$tract)
SD4<-as.character(mod[mod$var=="SD"&mod$mod_unit==4,]$tract)

TD1<-as.character(mod[mod$var=="TD"&mod$mod_unit==1,]$tract)
TD2<-as.character(mod[mod$var=="TD"&mod$mod_unit==2,]$tract)
TD3<-as.character(mod[mod$var=="TD"&mod$mod_unit==3,]$tract)
TD4<-as.character(mod[mod$var=="TD"&mod$mod_unit==4,]$tract)





PA1_val<-rowMeans(subset(PA_tr[complete.cases(PA_tr),],select=PA1))
PA2_val<-rowMeans(subset(PA_tr[complete.cases(PA_tr),],select=PA2))
PA3_val<-rowMeans(subset(PA_tr[complete.cases(PA_tr),],select=PA3))
PA4_val<-rowMeans(subset(PA_tr[complete.cases(PA_tr),],select=PA4))

PA_mods<-data.frame(PA1_val,PA2_val,PA3_val,PA4_val)
names(PA_mods)<-c("mod_1","mod_2","mod_3","mod_4")

PF1_val<-rowMeans(subset(PF_tr[complete.cases(PF_tr),],select=PF1))
PF2_val<-rowMeans(subset(PF_tr[complete.cases(PF_tr),],select=PF2))
PF3_val<-rowMeans(subset(PF_tr[complete.cases(PF_tr),],select=PF3))
PF4_val<-rowMeans(subset(PF_tr[complete.cases(PF_tr),],select=PF4))
PF5_val<-rowMeans(subset(PF_tr[complete.cases(PF_tr),],select=PF5))
PF6_val<-rowMeans(subset(PF_tr[complete.cases(PF_tr),],select=PF6))

PF_mods<-data.frame(PF1_val,PF2_val,PF3_val,PF4_val,PF5_val,PF6_val)
names(PF_mods)<-c("mod_1","mod_2","mod_3","mod_4""mod_5""mod_6")

SD1_val<-rowMeans(subset(SD_tr[complete.cases(SD_tr),],select=SD1))
SD2_val<-rowMeans(subset(SD_tr[complete.cases(SD_tr),],select=SD2))
SD3_val<-rowMeans(subset(SD_tr[complete.cases(SD_tr),],select=SD3))
SD4_val<-rowMeans(subset(SD_tr[complete.cases(SD_tr),],select=SD4))

SD_mods<-data.frame(SD1_val,SD2_val,SD3_val,SD4_val)
names(SD_mods)<-c("mod_1","mod_2","mod_3","mod_4")

TD1_val<-rowMeans(subset(TD_tr[complete.cases(TD_tr),],select=TD1))
TD2_val<-rowMeans(subset(TD_tr[complete.cases(TD_tr),],select=TD2))
TD3_val<-rowMeans(subset(TD_tr[complete.cases(TD_tr),],select=TD3))
TD4_val<-rowMeans(subset(TD_tr[complete.cases(TD_tr),],select=TD4))

TD_mods<-data.frame(TD1_val,TD2_val,TD3_val,TD4_val)
names(TD_mods)<-c("mod_1","mod_2","mod_3","mod_4")







