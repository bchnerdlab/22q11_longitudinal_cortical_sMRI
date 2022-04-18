library(readxl)
library(writexl)
library(tidyverse)
old_df<-read_excel('/Users/mariaj/Box/22q_dup_del_sMRI/data/REVISED_new_merged_data_WithSIPSratings_03022021.xlsx')
old_df2<-old_df%>%select("Scan_ID","PTID","SUBJECT_IDENTITY","SCANNER","AGEMONTH","SEX")

lh_area<-read_excel('/Users/mariaj/Box/22q_dup_del_sMRI/data/ROIs_06222021_mj.xlsx',sheet="lh_area")
rh_area<-read_excel('/Users/mariaj/Box/22q_dup_del_sMRI/data/ROIs_06222021_mj.xlsx',sheet="rh_area")
lh_thick<-read_excel('/Users/mariaj/Box/22q_dup_del_sMRI/data/ROIs_06222021_mj.xlsx',sheet="lh_thick")
rh_thick<-read_excel('/Users/mariaj/Box/22q_dup_del_sMRI/data/ROIs_06222021_mj.xlsx',sheet="rh_thick")

ind_roi<- lh_area %>% merge(rh_area,by="Scan_ID") %>% merge(lh_thick,by="Scan_ID") %>% merge(rh_thick,by="Scan_ID") 
ind_roi<- ind_roi %>% select(-PTID.x,-PTID.y)
#to merge with old df, need to remove anything after15 characters for scan ID
ind_roi$Scan_ID<-gsub("Q_0051_020302017","Q_0051_02032017",ind_roi$Scan_ID)
ind_roi$Scan_ID<-substr(ind_roi$Scan_ID,1,15)


ind_roi<-merge(ind_roi,old_df2,by="Scan_ID")
#remove excluded subjects
exclude<-which(ind_roi$`EXCLUDE?`=="Y")
ind_roi<-ind_roi[-exclude,]
write.table(ind_roi,file='/Users/mariaj/Box/22q_dup_del_sMRI/data/ind_roi_data_for_combat.txt',row.names = F,quote=F)

combat <- 
  read.table('/Users/mariaj/Box/22q_dup_del_sMRI/data/ind_roi_data_for_combat.txt',header=T) %>%  
  longCombat::longCombat(data=.,idvar='PTID',timevar='AGEMONTH',batchvar='SCANNER',features=c("lh_bankssts_area","lh_caudalanteriorcingulate_area","lh_caudalmiddlefrontal_area","lh_cuneus_area","lh_entorhinal_area", "lh_fusiform_area","lh_inferiorparietal_area","lh_inferiortemporal_area","lh_isthmuscingulate_area", "lh_lateraloccipital_area","lh_lateralorbitofrontal_area","lh_lingual_area","lh_medialorbitofrontal_area","lh_middletemporal_area","lh_parahippocampal_area", "lh_paracentral_area","lh_parsopercularis_area","lh_parsorbitalis_area","lh_parstriangularis_area","lh_pericalcarine_area","lh_postcentral_area","lh_posteriorcingulate_area",   "lh_precentral_area","lh_precuneus_area","lh_rostralanteriorcingulate_area","lh_rostralmiddlefrontal_area","lh_superiorfrontal_area","lh_superiorparietal_area", "lh_superiortemporal_area", "lh_supramarginal_area", "lh_frontalpole_area","lh_temporalpole_area","lh_transversetemporal_area","lh_insula_area","lh_WhiteSurfArea_area", "rh_bankssts_area","rh_caudalanteriorcingulate_area","rh_caudalmiddlefrontal_area","rh_cuneus_area","rh_entorhinal_area","rh_fusiform_area",   
                                                                                              "rh_inferiorparietal_area","rh_inferiortemporal_area","rh_isthmuscingulate_area", "rh_lateraloccipital_area","rh_lateralorbitofrontal_area", "rh_lingual_area",
                                                                                              "rh_medialorbitofrontal_area","rh_middletemporal_area","rh_parahippocampal_area","rh_paracentral_area","rh_parsopercularis_area","rh_parsorbitalis_area","rh_parstriangularis_area","rh_pericalcarine_area", "rh_postcentral_area","rh_posteriorcingulate_area","rh_precentral_area", "rh_precuneus_area","rh_rostralanteriorcingulate_area","rh_rostralmiddlefrontal_area", 
                                                                                              "rh_superiorfrontal_area","rh_superiorparietal_area","rh_superiortemporal_area","rh_supramarginal_area","rh_frontalpole_area","rh_temporalpole_area","rh_transversetemporal_area","rh_insula_area","rh_WhiteSurfArea_area","lh_bankssts_thickness", "lh_caudalanteriorcingulate_thickness","lh_caudalmiddlefrontal_thickness","lh_cuneus_thickness","lh_entorhinal_thickness", "lh_fusiform_thickness","lh_inferiorparietal_thickness","lh_inferiortemporal_thickness","lh_isthmuscingulate_thickness","lh_lateraloccipital_thickness","lh_lateralorbitofrontal_thickness",
                                                                                              "lh_lingual_thickness","lh_medialorbitofrontal_thickness","lh_middletemporal_thickness",  "lh_parahippocampal_thickness", "lh_paracentral_thickness","lh_parsopercularis_thickness", "lh_parsorbitalis_thickness","lh_parstriangularis_thickness","lh_pericalcarine_thickness","lh_postcentral_thickness","lh_posteriorcingulate_thickness","lh_precentral_thickness","lh_precuneus_thickness","lh_rostralanteriorcingulate_thickness","lh_rostralmiddlefrontal_thickness","lh_superiorfrontal_thickness","lh_superiorparietal_thickness","lh_superiortemporal_thickness","lh_supramarginal_thickness", "lh_frontalpole_thickness","lh_temporalpole_thickness", "lh_transversetemporal_thickness","lh_insula_thickness","lh_MeanThickness_thickness","rh_bankssts_thickness", "rh_caudalanteriorcingulate_thickness","rh_caudalmiddlefrontal_thickness","rh_cuneus_thickness","rh_entorhinal_thickness","rh_fusiform_thickness","rh_inferiorparietal_thickness","rh_inferiortemporal_thickness","rh_isthmuscingulate_thickness","rh_lateraloccipital_thickness","rh_lateralorbitofrontal_thickness","rh_lingual_thickness","rh_medialorbitofrontal_thickness","rh_middletemporal_thickness","rh_parahippocampal_thickness","rh_paracentral_thickness","rh_parsopercularis_thickness", "rh_parsorbitalis_thickness", "rh_parstriangularis_thickness","rh_pericalcarine_thickness","rh_postcentral_thickness","rh_posteriorcingulate_thickness","rh_precentral_thickness","rh_precuneus_thickness","rh_rostralanteriorcingulate_thickness", "rh_rostralmiddlefrontal_thickness","rh_superiorfrontal_thickness","rh_superiorparietal_thickness","rh_superiortemporal_thickness","rh_supramarginal_thickness", "rh_frontalpole_thickness","rh_temporalpole_thickness","rh_transversetemporal_thickness","rh_insula_thickness","rh_MeanThickness_thickness"), formula='SUBJECT_IDENTITY+AGEMONTH+SEX',ranef='(1|PTID)',niter=100)
                                    
                                     

adj_df<-combat$data_combat
adj_df<-cbind(ind_roi$Scan_ID,adj_df)
colnames(adj_df)[1]<-"Scan_ID"
write_xlsx(adj_df,'/Users/mariaj/Desktop/22q_dup_del_sMRI/data/20210628_aftercombat_ind_roi.xlsx')
