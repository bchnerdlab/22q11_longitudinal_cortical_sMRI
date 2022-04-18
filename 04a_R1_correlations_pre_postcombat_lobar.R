

adj_df<-read_xlsx('/Users/mariajalbrzikowski/Dropbox (BCH)/22q_dup_del_sMRI/data/20210628_aftercombat_data.xlsx')

brain_reg<-c("lh_frontal_thickness","lh_occipital_thickness",
"lh_parietal_thickness","lh_temporal_thickness","rh_frontal_thickness","rh_occipital_thickness",
"rh_parietal_thickness","rh_temporal_thickness","lh_frontal_area" ,"lh_occipital_area","lh_parietal_area",
"lh_temporal_area","rh_frontal_area","rh_occipital_area","rh_parietal_area","rh_temporal_area","lh_TotalSA_area",
"rh_TotalSA_area","lh_MeanThickness_thickness","rh_MeanThickness_thickness","EstimatedTotalIntraCranialVol")

for (i in c(1:length(brain_reg)))
{
  var1=brain_reg[i]
  var2<-paste(brain_reg[i],"combat",sep=".")
corres<-cor.test(adj_df[[var1]],adj_df[[var2]])#cor 0.9707882 
print(corres$estimate)
}


