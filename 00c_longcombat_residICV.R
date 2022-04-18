
library(magrittr)
library(readxl)
library(longCombat)
library(writexl)



df<-read_excel('/Users/jalbrzikowskime/Dropbox (BCH)/22q_dup_del_sMRI/data/REVISED_new_merged_data_WithSIPSratings_03022021.xlsx')


df$SUBJECT_IDENTITY<-as.factor(df$SUBJECT_IDENTITY)
df$SCANNER<-as.factor(df$SCANNER)
df$PTID<-as.factor(df$PTID)
#remove excluded subjects
exclude<-which(df$`EXCLUDE?`=="Y")
df<-df[-exclude,]

#for some unknown reason, you need to select only the variables you are going to use in longcombat, save it, then read it back in, this is the only way I got it to work
keep<-c("PTID","SUBJECT_IDENTITY","SCANNER","AGEMONTH","SEX","lh_frontal_thickness","lh_occipital_thickness",
        "lh_parietal_thickness","lh_temporal_thickness","rh_frontal_thickness","rh_occipital_thickness",
        "rh_parietal_thickness","rh_temporal_thickness","lh_frontal_area" ,"lh_occipital_area","lh_parietal_area",
        "lh_temporal_area","rh_frontal_area","rh_occipital_area","rh_parietal_area","rh_temporal_area","lh_TotalSA_area",
        "rh_TotalSA_area","lh_MeanThickness_thickness","rh_MeanThickness_thickness","EstimatedTotalIntraCranialVol")
keepcol<-match(keep,colnames(df))
df<-df[,keepcol]
write.table(df,file='/Users/jalbrzikowskime/Dropbox (BCH)/22q_dup_del_sMRI/data/data_for_combat_ICVresid.txt',row.names = F,quote=F)

combat <- 
  read.table('/Users/jalbrzikowskime/Dropbox (BCH)/22q_dup_del_sMRI/data/data_for_combat_ICVresid.txt',header=T) %>%  
  longCombat::longCombat(data=.,idvar='PTID',timevar='AGEMONTH',batchvar='SCANNER',features=c("lh_frontal_thickness","lh_occipital_thickness",
                                                                                              "lh_parietal_thickness","lh_temporal_thickness","rh_frontal_thickness","rh_occipital_thickness",
                                                                                              "rh_parietal_thickness","rh_temporal_thickness","lh_frontal_area" ,"lh_occipital_area","lh_parietal_area",
                                                                                              "lh_temporal_area","rh_frontal_area","rh_occipital_area","rh_parietal_area","rh_temporal_area","lh_TotalSA_area",
                                                                                              "rh_TotalSA_area","lh_MeanThickness_thickness","rh_MeanThickness_thickness","EstimatedTotalIntraCranialVol"),
                         formula='SUBJECT_IDENTITY+AGEMONTH+SEX',ranef='(1|PTID)',niter=100)

adj_df<-combat$data_combat
adj_df<-adj_df[,4:24]

adj_df<-cbind(df,adj_df)
adj_df$Age<-adj_df$AGEMONTH/12
write_xlsx(adj_df,'/Users/jalbrzikowskime/Dropbox (BCH)/22q_dup_del_sMRI/data/20220413_aftercombat_data.xlsx')

