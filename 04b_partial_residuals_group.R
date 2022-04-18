
library(mgcv)
library(gratia)
library(ggplot2)
library(patchwork)
adj_df<-read_xlsx("/Users/jalbrzikowskime/Dropbox (BCH)/22q_dup_del_sMRI/data/20210628_aftercombat_data.xlsx")
adj_df$PTID<-factor(adj_df$PTID)


adj_df$Group <- factor(adj_df$SUBJECT_IDENTITY)
adj_df$PTID<-factor(adj_df$PTID)
adj_df$Mean_thickness.combat<-adj_df$lh_MeanThickness_thickness.combat+adj_df$rh_MeanThickness_thickness.combat
adj_df$TotalSA_area.combat<-adj_df$lh_TotalSA_area.combat+adj_df$rh_TotalSA_area.combat


brain_regions<-c("Mean_thickness.combat","lh_frontal_thickness.combat","rh_frontal_thickness.combat","lh_temporal_thickness.combat",
                 "rh_temporal_thickness.combat","lh_parietal_thickness.combat","rh_parietal_thickness.combat",
                 "lh_occipital_thickness.combat","rh_occipital_thickness.combat","TotalSA_area.combat","lh_frontal_area.combat","rh_frontal_area.combat",
                 "lh_temporal_area.combat","rh_temporal_area.combat","lh_parietal_area.combat",
                 "rh_parietal_area.combat","lh_occipital_area.combat","rh_occipital_area.combat")
plotnames<-c("Mean CT", "LH Frontal CT","RH Frontal CT","LH Temporal CT","RH Temporal CT","LH Parietal CT","RH Parietal CT","LH Occipital CT","RH Occipital CT",
             "Total SA","LH Frontal SA","RH Frontal SA","LH Temporal SA","RH Temporal SA","LH Parietal SA","RH Parietal SA","LH Occipital SA","RH Occipital SA")

pdf('/Users/jalbrzikowskime/Dropbox (BCH)/22q_dup_del_sMRI/results/Group_partial_residuals.pdf',width=7,height=8)
for (i in c(1:length(brain_regions))){
  
  var=brain_regions[i]
  mod<-gam(adj_df[[var]]~s(Age,by=Group,k=4)+Group+SEX+s(PTID,bs="re",k=4)+EstimatedTotalIntraCranialVol.combat,selection=TRUE,method="REML",data=adj_df)
  p<-draw(mod,residuals = T,scale="fixed",select=c(1,2,3))
  p2<-p+plot_annotation(title=plotnames[[i]])
  print(p2)
}
dev.off()