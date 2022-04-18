knitr::opts_chunk$set(echo = TRUE)
library(tidyverse)
library(mgcv)
library(gratia)
library(readxl)
library(writexl)
library(patchwork)
library(cowplot)
library(extrafont)

font_import(paths = "/Users/jalbrzikowskime/Dropbox (BCH)/22q_dup_del_sMRI/code/Arialfont",prompt=T)

adj_df<-read_xlsx('/Users/jalbrzikowskime/Dropbox (BCH)/22q_dup_del_sMRI/data/20210628_aftercombat_data.xlsx')
adj_df$Mean_thickness.combat<-adj_df$lh_MeanThickness_thickness.combat+adj_df$rh_MeanThickness_thickness.combat
adj_df$TotalSA_area.combat<-adj_df$lh_TotalSA_area.combat+adj_df$rh_TotalSA_area.combat


#22q-del ASD specific
del22q<-adj_df[adj_df$SUBJECT_IDENTITY=="PATIENT-DEL",]

#ASD was only formally assessed at baseline , so usually we just carried that diagnosis forward if they had it at baseline but sometimes it was sort of ‘iffy’ (like they met  it on ADOS but not ADI but were also really quiet, or something, then at followup they seemed very not ASD ). So if they got ASD at the majority of visits I would give that a ‘yes’ , and that just leaves q_0200 as a ? so for that one I’d say if they got it at baseline they should get the ASD dx. does that sound reasonable? we couldn’t really be consistent if we wanted to bc we assess ASD spectrum at each visit, and not ASD , so – there you go. hope that helps!
#q_0260: 1x = 0, 3x =1 1x=3 (*x refers to numbers of visits they received that dx)
#q_0001: 1x =0; 2x =1
#q_0200: 1x=0, 1x=1
#q_200 v1= 0, v2=3, v3=1
q0001_recode<-which(del22q$PTID=="q_0001")
del22q$ASD.DIAGNOS[q0001_recode]<-1

q0200_recode<-which(del22q$PTID=="q_0200")
del22q$ASD.DIAGNOS[q0200_recode]<-0

q0260_recode<-which(del22q$PTID=="q_0260")
del22q$ASD.DIAGNOS[q0260_recode]<-1

#check for any labeled 3
check<-table(del22q$PTID,del22q$ASD.DIAGNOS)


#columns 0 1 3
# q_0114 0 1 2
q0114_recode<-which(del22q$PTID=="q_0114")
del22q$ASD.DIAGNOS[q0114_recode]<-1

# q_0124 0 1 2
q0124_recode<-which(del22q$PTID=="q_0124")
del22q$ASD.DIAGNOS[q0124_recode]<-0

#q_0190 1 0 2
q0190_recode<-which(del22q$PTID=="q_0190")
del22q$ASD.DIAGNOS[q0190_recode]<-0


#q_0217 4 0 1
q_0217_recode<-which(del22q$PTID=="q_0217")
del22q$ASD.DIAGNOS[q_0217_recode]<-0


# q_0244 0 1 1
q_0244_recode<-which(del22q$PTID=="q_0244")
del22q$ASD.DIAGNOS[q_0244_recode]<-1

# q_0246	2	0	2
q_0246_recode<-which(del22q$PTID=="q_0246")
del22q$ASD.DIAGNOS[q_0246_recode]<-0

#q_0355 0 1 1
q_0355_recode<-which(del22q$PTID=="q_0355")
del22q$ASD.DIAGNOS[q_0355_recode]<-1



check2<-table(del22q$PTID,del22q$ASD.DIAGNOS)
del22q<-del22q[del22q$ASD.DIAGNOS!=3,]

del22q$PTID<-as.factor(del22q$PTID)
del22q$ASD.DIAGNOS<-as.factor(del22q$ASD.DIAGNOS)

un<-del22q[!duplicated(del22q$PTID),]


del22q<-del22q[del22q$Age<26,]

brain_regions<-c("Mean_thickness.combat","lh_frontal_thickness.combat","rh_frontal_thickness.combat","lh_temporal_thickness.combat",
                 "rh_temporal_thickness.combat","lh_parietal_thickness.combat","rh_parietal_thickness.combat",
                 "lh_occipital_thickness.combat","rh_occipital_thickness.combat","TotalSA_area.combat","lh_frontal_area.combat","rh_frontal_area.combat",
                 "lh_temporal_area.combat","rh_temporal_area.combat","lh_parietal_area.combat",
                 "rh_parietal_area.combat","lh_occipital_area.combat","rh_occipital_area.combat")
plotnames<-c("Mean CT", "LH Frontal CT","RH Frontal CT","LH Temporal CT","RH Temporal CT","LH Parietal CT","RH Parietal CT","LH Occipital CT","RH Occipital CT",
             "Total SA","LH Frontal SA","RH Frontal SA","LH Temporal SA","RH Temporal SA","LH Parietal SA","RH Parietal SA","LH Occipital SA","RH Occipital SA")

pdf('/Users/jalbrzikowskime/Dropbox (BCH)/22q_dup_del_sMRI/22qdelASD_partial_residuals.pdf',width=7,height=4)
for (i in c(1:length(brain_regions))){
  
  var=brain_regions[i]
  mod<-gam(del22q[[var]]~s(Age,by=ASD.DIAGNOS,k=4)+
             ASD.DIAGNOS+SEX+SCANNER+EstimatedTotalIntraCranialVol.combat+s(PTID,bs="re"),selection=TRUE,method="REML",data=del22q)
  p<-draw(mod,residuals = T,scale="fixed",select=c(1,2))
  p2<-p+plot_annotation(title=plotnames[[i]])
  print(p2)
}
dev.off()

dup22q<-adj_df[adj_df$SUBJECT_IDENTITY=="PATIENT-DUP",]
table(dup22q$PTID,dup22q$ASD.DIAGNOS)
#q_0402 had a 3 for not assessed for ASD, need to remove
rem=which(dup22q$PTID=="q_0402")
dup22q=dup22q[-rem,]
dup22q<-dup22q[dup22q$Age<26,]





dup22q$PTID<-as.factor(dup22q$PTID)
dup22q$ASD.DIAGNOS<-as.factor(dup22q$ASD.DIAGNOS)
pdf('/Users/jalbrzikowskime/Dropbox (BCH)/22q_dup_del_sMRI/22qdupASD_partial_residuals.pdf',width=7,height=4)
for (i in c(1:length(brain_regions))){
  
  var=brain_regions[i]
  mod<-gam(dup22q[[var]]~s(Age,by=ASD.DIAGNOS,k=4)+
             ASD.DIAGNOS+SEX+SCANNER+EstimatedTotalIntraCranialVol.combat+s(PTID,bs="re"),selection=TRUE,method="REML",data=dup22q)
  p<-draw(mod,residuals = T,scale="fixed",select=c(1,2))
  p2<-p+plot_annotation(title=plotnames[[i]])
  print(p2)
}

dev.off()
