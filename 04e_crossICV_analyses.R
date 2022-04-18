library(tidyverse)
library(mgcv)
library(readxl)
library(writexl)
library(patchwork)
library(lme4)



adj_df<-read_xlsx('/Users/mariajalbrzikowski/Dropbox (BCH)/22q_dup_del_sMRI/data/20210628_aftercombat_data.xlsx')
totalICV<-read.csv('/Users/mariajalbrzikowski/Dropbox (BCH)/22q_dup_del_sMRI/data/total_cranial_volumes_cross_sectional_v2.csv')
adj_df$Scan_ID<-as.factor(adj_df$Scan_ID)
totalICV$Scan_ID<-as.factor(totalICV$session)

adj_df$Mean_thickness.combat<-adj_df$lh_MeanThickness_thickness.combat+adj_df$rh_MeanThickness_thickness.combat
adj_df$TotalSA_area.combat<-adj_df$lh_TotalSA_area.combat+adj_df$rh_TotalSA_area.combat

comb<-merge(adj_df,totalICV,by="Scan_ID",all.x=T)
comb$Group<-as.factor(comb$SUBJECT_IDENTITY)
comb$PTID<-as.factor(comb$PTID)
comb$CONVERTEDVISITNUM<-as.factor(comb$CONVERTEDVISITNUM)


#CNP controls were cross-sectional only so take their data from the longitudinal sheet
CNP<-grep("sub_",comb$Scan_ID)
comb$volume[CNP]<-comb$EstimatedTotalIntraCranialVol[CNP]

comb<-comb[order(comb$PTID,comb$TESTDATE...9),]

try<-comb[comb$PTID %in% comb$PTID[duplicated(comb$PTID)],]


#try<-comb
#make visit variable a factor
try2<-try %>% group_by(PTID) %>% mutate(visit=rank(CONVERTEDVISITNUM))

try2$visit<-as.factor(try2$visit)
try2<-try2[try2$visit<6,]
try2$visit<-as.factor(try2$visit)
#Volume is the estimated ICV variable that Chris & Leila gave me from the cross-sectional output
mod<-gam(try2$volume~s(Age,by=Group,k=4)+Group+
           visit+SEX+SCANNER+s(PTID,bs="re"),selection=TRUE,method="REML",data=try2)
summary(mod)
sum_mod<-as.data.frame(sum_mod)

#get partial residuals for plotting
mod<-gam(try2$volume~s(Age,by=Group,k=4)+Group+
           SEX+SCANNER+s(PTID,bs="re"),selection=TRUE,method="REML",data=try2)
try2$vol_resid<-mod$residuals

try2$visit<-as.numeric(try2$visit)

  try2$visit=as.numeric(try2$visit)
  try2$Group<-factor(try2$Group,levels=c("PATIENT-DEL","CONTROL","PATIENT-DUP"))
ggplot(try2,aes(x=visit,y=vol_resid)) + geom_point(aes(group=as.factor(Group),color=as.factor(Group)),alpha=0.2)+
  geom_line(data = try2,aes_string(x = "visit", y = "vol_resid",group="PTID",color = "Group"),size = 1,alpha=0.2)+
  geom_smooth(aes(color=as.factor(Group)),method="loess",span=2)+theme_bw()+ylab("Partial residuals of Estimated ICV-cross sectional")+
  scale_fill_viridis(discrete=TRUE)+scale_color_viridis(discrete = TRUE)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("Visit")
  
  
ggplot(try2,aes(x=visit,y=volume)) + geom_point(aes(group=as.factor(Group),color=as.factor(Group)),alpha=0.2)+
  geom_line(data = try2,aes_string(x = "visit", y = "volume",group="PTID",color = "Group"),size = 1,alpha=0.2)+
  geom_smooth(aes(color=as.factor(Group)),method="loess",span=2)+theme_bw()+ylab("Estimated ICV-cross sectional")+
  scale_fill_viridis(discrete=TRUE)+scale_color_viridis(discrete = TRUE)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("Visit")

+theme(legend.position=c(.9,.75))

  
ggplot(try2, aes(x=visit, y=vol_resid)) +
  geom_point() + geom_smooth(method="loess") 

  theme_bw() 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("Visit") + 
  ylab("Est ICV")
  
  

  geom_point(data = try2,aes_string(x = "visit", y = "vol_resid",color = "Group"),size = 2)+

  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  scale_fill_viridis(discrete=TRUE)+scale_color_viridis(discrete = TRUE)

ggplot(comb,aes(x=visit,y=vol_resid )) +  geom_point(aes(group=Group,color=Group))+
  geom_line(aes(group=paste(comb$Group,comb$PTID),color=adj_df$Group,fill=adj_df$Group))




ggplot(data = comb, aes_string(x = "CONVERTEDVISITNUM",y = "volume") )+
  geom_line(data = comb,aes_string(x = "CONVERTEDVISITNUM", y = "volume",group="PTID",color = "Group"),size = 1)+
geom_point(data = comb,aes_string(x = "CONVERTEDVISITNUM", y = "volume",color = "Group"),size = 1)+
  geom_smooth(data = comb,aes_string(x = "CONVERTEDVISITNUM", y = "volume",color = "Group"))
  
  


  
  
 