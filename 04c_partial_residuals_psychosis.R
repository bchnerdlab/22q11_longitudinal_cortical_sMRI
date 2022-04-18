##make sure you've run the following:
#library(remotes)
#remotes::install_version("Rttf2pt1", version = "1.3.8")

library(ggplot2)
library(patchwork)
library(viridis)
library(tidyverse)
library(extrafont)
library(gratia)
library(readxl)
library(mgcv)

adj_df<-read_xlsx("/Users/jalbrzikowskime/Box/22q_dup_del_sMRI/data/20210628_aftercombat_data.xlsx")
adj_df$PTID<-factor(adj_df$PTID)

adj_df$Group <- factor(adj_df$SUBJECT_IDENTITY)
adj_df$PTID<-factor(adj_df$PTID)
adj_df$Mean_thickness.combat<-adj_df$lh_MeanThickness_thickness.combat+adj_df$rh_MeanThickness_thickness.combat
adj_df$TotalSA_area.combat<-adj_df$lh_TotalSA_area.combat+adj_df$rh_TotalSA_area.combat

#22q-del psychosis specific
del22q<-adj_df[adj_df$SUBJECT_IDENTITY=="PATIENT-DEL",]#| adj_df$SUBJECT_IDENTITY=="CONTROL",]
#according to SIPSPSY.code.definition, 3 means "not assessed"
#removed those with 3 and 22qdel
rem<-which(del22q$SIPSPSY.code==3 & del22q$SUBJECT_IDENTITY=="PATIENT-DEL")
del22q<-del22q[-rem,]
#prodromal/hr=1, psychotic =2
#create psychosis spectrum
del22q$psychosis_spect<-gsub(2,1,del22q$SIPSPSY.code)
del22q$psychosis_spect<-gsub(1,"psychosis_spect",del22q$psychosis_spect)
del22q$psychosis_spect<-gsub(0,"nopsychosis",del22q$psychosis_spect)
#del22q$psychosis_spect<-as.factor(del22q$psychosis_spect)
cont<-which(del22q$SUBJECT_IDENTITY=="CONTROL")
del22q$psychosis_spect[cont]<-"control"
#check consistency of psychosis spectrum
table(del22q$psychosis_spect,del22q$PTID)
#psychosis_spect1= keep it the original way that the data is coded
del22q$psychosis_spect1<-del22q$psychosis_spect
#psychosis_spect2= recode those who are psychosis spectrum at one point as psychosis spectrum at ALL time points
del22q$psychosis_spect2<-del22q$psychosis_spect


#these are individauls who are sometimes psychosis spectrum and sometimes not psychosis spectrum
#we are recoding them as psychosis_spect (2) for all visits
recode<-c("q_0001","q_0005","q_0016","q_0021","q_0130","q_0234",
          "q_0238","q_0260","q_0263")
rec<-which(del22q$PTID=="q_0001")
del22q$psychosis_spect2[rec]<-"psychosis_spect"
rec<-which(del22q$PTID=="q_0005")
del22q$psychosis_spect2[rec]<-"psychosis_spect"
rec<-which(del22q$PTID=="q_0016")
del22q$psychosis_spect2[rec]<-"psychosis_spect"
rec<-which(del22q$PTID=="q_0021")
del22q$psychosis_spect2[rec]<-"psychosis_spect"
rec<-which(del22q$PTID=="q_0130")
del22q$psychosis_spect2[rec]<-"psychosis_spect"
rec<-which(del22q$PTID=="q_0234")
del22q$psychosis_spect2[rec]<-"psychosis_spect"
rec<-which(del22q$PTID=="q_0238")
del22q$psychosis_spect2[rec]<-"psychosis_spect"
rec<-which(del22q$PTID=="q_0260")
del22q$psychosis_spect2[rec]<-"psychosis_spect"
rec<-which(del22q$PTID=="q_0263")
del22q$psychosis_spect2[rec]<-"psychosis_spect"

del22q$psychosis_spect2<-as.factor(del22q$psychosis_spect2)

#del22q<-del22q[del22q$Age<26,]

brain_regions<-c("Mean_thickness.combat","lh_frontal_thickness.combat","rh_frontal_thickness.combat","lh_temporal_thickness.combat",
                 "rh_temporal_thickness.combat","lh_parietal_thickness.combat","rh_parietal_thickness.combat",
                 "lh_occipital_thickness.combat","rh_occipital_thickness.combat","TotalSA_area.combat","lh_frontal_area.combat","rh_frontal_area.combat",
                 "lh_temporal_area.combat","rh_temporal_area.combat","lh_parietal_area.combat",
                 "rh_parietal_area.combat","lh_occipital_area.combat","rh_occipital_area.combat")
plotnames<-c("Mean CT", "LH Frontal CT","RH Frontal CT","LH Temporal CT","RH Temporal CT","LH Parietal CT","RH Parietal CT","LH Occipital CT","RH Occipital CT",
             "Total SA","LH Frontal SA","RH Frontal SA","LH Temporal SA","RH Temporal SA","LH Parietal SA","RH Parietal SA","LH Occipital SA","RH Occipital SA")

pdf('/Users/jalbrzikowskime/Dropbox (BCH)/22q_dup_del_sMRI/results/Psychosis_partial_residuals.pdf',width=7,height=4)
for (i in c(1:length(brain_regions))){
  
  var=brain_regions[i]
  mod<-gam(del22q[[var]]~s(Age,by=psychosis_spect2,k=4)+psychosis_spect2+SEX+s(PTID,bs="re",k=4)+EstimatedTotalIntraCranialVol.combat,selection=TRUE,method="REML",data=del22q)
  p<-draw(mod,residuals = T,scale="fixed",select=c(1,2))
  p2<-p+plot_annotation(title=plotnames[[i]])
  print(p2)
}
dev.off()
