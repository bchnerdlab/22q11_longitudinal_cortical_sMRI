library(patchwork)
library(viridis)
library(tidyverse)
library(gratia)


adj_df<-read_xlsx("/Users/jalbrzikowskime/Dropbox (BCH)/22q_dup_del_sMRI/data/20210628_aftercombat_data.xlsx")
adj_df$PTID<-factor(adj_df$PTID)

adj_df$Group <- factor(adj_df$SUBJECT_IDENTITY)
adj_df$PTID<-factor(adj_df$PTID)
adj_df$Mean_thickness.combat<-adj_df$lh_MeanThickness_thickness.combat+adj_df$rh_MeanThickness_thickness.combat
adj_df$TotalSA_area.combat<-adj_df$lh_TotalSA_area.combat+adj_df$rh_TotalSA_area.combat



mod<-gam(adj_df$Mean_thickness.combat~s(Age,by=Group,k=4)+Group+SEX+s(PTID,bs="re",k=4)+EstimatedTotalIntraCranialVol.combat,selection=TRUE,method="REML",data=adj_df)
lmod<-lm(adj_df$Mean_thickness.combat~Group+SEX+EstimatedTotalIntraCranialVol.combat,data=adj_df)
adj_df$Mean_thickness.combat.res<-lmod$residuals
try<-evaluate_smooth(mod,"s(Age)",n=44)

try$selo <- try$est - try$se
try$sehi <- try$est + try$se

line_size=2
try$Group<-factor(try$Group,levels=c("PATIENT-DEL","CONTROL","PATIENT-DUP"))
p1<-ggplot(data = try, aes_string(x = "Age",y = "est") )+
  geom_ribbon(data = try,aes_string(x = "Age", ymin = "selo",ymax = "sehi", fill = "Group"),alpha = .18, linetype = 0)+
  geom_line(data = try,aes_string(x = "Age", y = "est",color = "Group"),size = 1)+theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  geom_point(data = adj_df,aes_string(x = "Age", y = "Mean_thickness.combat.res",color="Group"),size = 1,alpha=0.4)+
  geom_line(data = adj_df,aes_string(x = "Age", y = "Mean_thickness.combat.res",group="PTID",color = "Group"),size = 0.2,alpha=0.4)+
  theme(legend.title = element_blank())+ylab("Mean Cortical Thickness")+scale_fill_viridis(discrete=TRUE)+scale_color_viridis(discrete = TRUE)+ggtitle("A")

derv<-derivatives(mod,term="s(Age)",n=44,partial_match = T)
derv$t<-derv$derivative/derv$se
cont=derv[derv$smooth=="s(Age):GroupCONTROL",]
cont$t<-cont$derivative/cont$se
del=derv[derv$smooth=="s(Age):GroupPATIENT-DEL",]
del$t<-del$derivative/del$se

dup=derv[derv$smooth=="s(Age):GroupPATIENT-DUP",]
dup$t<-dup$derivative/dup$se

tile_con <- ggplot(cont) + aes_string(x = "data", y = 1,   
                                      fill = "derivative") + geom_raster(interpolate = FALSE)+scale_fill_gradient2(high="#dd4124",low="#00496f", mid="white", limits = c(min(derv$derivative),max(derv$derivative))) + xlab("")+ ylab("Control")+theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +theme(legend.position = "none")+ theme(axis.ticks.y = element_blank(),
                                                                                                                        axis.text.y = element_blank())

tile_del <- ggplot(del) + aes_string(x = "data", y = 1,   
                                     fill = "derivative") + geom_raster(interpolate = FALSE)+scale_fill_gradient2(high="#dd4124",low="#00496f", mid="white", limits = c(min(derv$derivative),max(derv$derivative))) + xlab("")+ ylab("22qDel")+theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme(axis.ticks.y = element_blank(),
                                                                                       axis.text.y = element_blank())

tile_dup <- ggplot(dup) + aes_string(x = "data", y = 1,   
                                     fill = "derivative") + geom_raster(interpolate = FALSE)+scale_fill_gradient2(high="#dd4124",low="#00496f", mid="white", limits = c(min(derv$derivative),max(derv$derivative))) + xlab("")+ ylab("22qDup")+theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +theme(legend.position = "none")+ theme(axis.ticks.y = element_blank(),
                                                                                                                        axis.text.y = element_blank())
pdf('/Users/jalbrzikowskime/Dropbox (BCH)/22q_dup_del_sMRI/results/Figure02_age_CT_wpoints.pdf',width=7,height=8)
print((p1|tile_con|tile_del|tile_dup)+plot_layout(nrow=4,heights = c(5,1,1,1)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
dev.off()




mod<-gam(adj_df$TotalSA_area.combat~s(Age,by=Group,k=4)+Group+SEX+s(PTID,bs="re",k=4)+EstimatedTotalIntraCranialVol.combat,selection=TRUE,method="REML",data=adj_df)
lmod<-lm(adj_df$TotalSA_area.combat~Group+SEX+EstimatedTotalIntraCranialVol.combat,data=adj_df)
adj_df$TotalSA_area.combat.res<-lmod$residuals

try<-evaluate_smooth(mod,"s(Age)",n=44)

try$selo <- try$est - try$se
try$sehi <- try$est + try$se



line_size=2
try$Group<-factor(try$Group,levels=c("PATIENT-DEL","CONTROL","PATIENT-DUP"))
p2<-ggplot(data = try, aes_string(x = "Age",y = "est") )+
  geom_ribbon(data = try,aes_string(x = "Age", ymin = "selo",ymax = "sehi", fill = "Group"),alpha = .18, linetype = 0)+
  geom_line(data = try,aes_string(x = "Age", y = "est",color = "Group"),size = 1)+theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  geom_point(data = adj_df,aes_string(x = "Age", y = "TotalSA_area.combat.res",color="Group"),size = 1,alpha=0.4)+
  geom_line(data = adj_df,aes_string(x = "Age", y = "TotalSA_area.combat.res",group="PTID",color = "Group"),size = 0.2,alpha=0.4)+
  theme(legend.title = element_blank())+ylab("Total Surface Area")+scale_fill_viridis(discrete=TRUE)+scale_color_viridis(discrete = TRUE)+ggtitle("B")

derv<-derivatives(mod,term="Age",n=44,partial_match = T)
derv$t<-derv$derivative/derv$se
cont=derv[derv$smooth=="s(Age):GroupCONTROL",]
cont$t<-cont$derivative/cont$se
del=derv[derv$smooth=="s(Age):GroupPATIENT-DEL",]
del$t<-del$derivative/del$se

dup=derv[derv$smooth=="s(Age):GroupPATIENT-DUP",]
dup$t<-dup$derivative/dup$se

tile_con2 <- ggplot(cont) + aes_string(x = "data", y = 1,   
                                      fill = "derivative") + geom_raster(interpolate = FALSE)+scale_fill_gradient2(high="#dd4124",low="#00496f", mid="white", limits = c(min(derv$derivative),max(derv$derivative))) + xlab("")+ ylab("Control")+theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +theme(legend.position = "none")+ theme(axis.ticks.y = element_blank(),
                                                                                                                        axis.text.y = element_blank())

tile_del2 <- ggplot(del) + aes_string(x = "data", y = 1,   
                                     fill = "derivative") + geom_raster(interpolate = FALSE)+scale_fill_gradient2(high="#dd4124",low="#00496f", mid="white", limits = c(min(derv$derivative),max(derv$derivative))) + xlab("")+ ylab("22qDel")+theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme(axis.ticks.y = element_blank(),
                                                                                       axis.text.y = element_blank())

tile_dup2 <- ggplot(dup) + aes_string(x = "data", y = 1,   
                                     fill = "derivative") + geom_raster(interpolate = FALSE)+scale_fill_gradient2(high="#dd4124",low="#00496f", mid="white", limits = c(min(derv$derivative),max(derv$derivative))) + xlab("")+ ylab("22qDup")+theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +theme(legend.position = "none")+ theme(axis.ticks.y = element_blank(),
                                                                                                                        axis.text.y = element_blank())


pdf('/Users/jalbrzikowskime/Dropbox (BCH)/22q_dup_del_sMRI/results/Figure02_age_SA_wpoints.pdf',width=7,height=8)
print((p2|tile_con2|tile_del2|tile_dup2)+plot_layout(nrow=4,heights = c(5,1,1,1)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
dev.off()