library(magrittr)
library(readxl)
library(writexl)
library(ggplot2)
library(mgcv)
library(gratia)
library(patchwork)
library(viridis)
library(tidyverse)
library(extrafont)
library(cowplot)
library(egg)


adj_df<-read_xlsx("/Users/rahay/Box/22q_dup_del_sMRI/data/20210628_aftercombat_data.xlsx")
adj_df$PTID<-factor(adj_df$PTID)

brain_regions<-c("lh_frontal_thickness.combat","rh_frontal_thickness.combat","lh_temporal_thickness.combat","rh_temporal_thickness.combat","lh_parietal_thickness.combat","rh_parietal_thickness.combat","lh_occipital_thickness.combat","rh_occipital_thickness.combat","lh_frontal_area.combat","rh_frontal_area.combat",
                 "lh_temporal_area.combat","rh_temporal_area.combat","lh_parietal_area.combat","rh_parietal_area.combat","lh_occipital_area.combat","rh_occipital_area.combat")
plotnames<-c("LH Frontal","RH Frontal","LH Temporal","RH Temporal","LH Parietal","RH Parietal","LH Occipital","RH Occipital","LH Frontal","RH Frontal","LH Temporal","RH Temporal","LH Parietal","RH Parietal","LH Occipital","RH Occipital")
axislab<-c("mm","mm","mm","mm","mm","mm","mm","mm",expression(mm^2),expression(mm^2),expression(mm^2),expression(mm^2),expression(mm^2),expression(mm^2),expression(mm^2),expression(mm^2))


adj_df$Group <- factor(adj_df$SUBJECT_IDENTITY)
adj_df$PTID<-factor(adj_df$PTID)
adj_df$Mean_thickness.combat<-adj_df$lh_MeanThickness_thickness.combat+adj_df$rh_MeanThickness_thickness.combat
adj_df$TotalSA_area.combat<-adj_df$lh_TotalSA_area.combat+adj_df$rh_TotalSA_area.combat

plots<-list()

for (i in c(1:length(brain_regions))){
  
  var=brain_regions[i]
mod<-gam(adj_df[[var]]~s(Age,by=Group,k=4)+Group+SEX+s(PTID,bs="re",k=4)+EstimatedTotalIntraCranialVol.combat,selection=TRUE,method="REML",data=adj_df)


try<-evaluate_smooth(mod,"s(Age)",n=44)

try$selo <- try$est - try$se
try$sehi <- try$est + try$se

line_size=2
try$Group<-factor(try$Group,levels=c("PATIENT-DEL","CONTROL","PATIENT-DUP"))
p2<-ggplot(data = try, aes_string(x = "Age",y = "est"))+
  geom_ribbon(data = try,aes_string(x = "Age", ymin = "selo",ymax = "sehi", fill = "Group"),alpha = .18, linetype = 0)+
  geom_line(data = try,aes_string(x = "Age", y = "est",color = "Group"),size = 1)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  theme(legend.title = element_blank())+ylab(axislab[[i]])+xlab("")+
  scale_fill_viridis(discrete=TRUE)+scale_color_viridis(discrete = TRUE)+
  theme(legend.position = "None") + theme(plot.margin = margin(t = 0,  r = 0, b = 0, l = 0)) +
  ggtitle(plotnames[[i]]) +theme(text = element_text(family = "Arial"))+coord_cartesian(xlim = c(5, 50), expand = T)
derv<-derivatives(mod,term="Age",n=44,partial_match = TRUE)
#derv1<-derivatives(mod,term=c("s(Age):GroupCONTROL","s(Age):GroupPATIENT-DEL","s(Age):GroupPATIENT-DUP"),n=44)


derv$t<-derv$derivative/derv$se
cont=derv[derv$smooth=="s(Age):GroupCONTROL",]
cont$t<-cont$derivative/cont$se
del=derv[derv$smooth=="s(Age):GroupPATIENT-DEL",]
del$t<-del$derivative/del$se

dup=derv[derv$smooth=="s(Age):GroupPATIENT-DUP",]
dup$t<-dup$derivative/dup$se

tile_con2 <- ggplot(cont) + aes_string(x = "data", y = 1,   
                                      fill = "derivative") + geom_raster(interpolate = FALSE)+scale_fill_gradient2(high="#dd4124",low="#00496f", mid="white", limits = c(min(derv$derivative),max(derv$derivative)),guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) + xlab("")+ ylab("")+theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_y_continuous(breaks = seq(1), labels = "Control") +
  theme(legend.position = "none")+ theme(axis.ticks.y = element_blank(),axis.text.x = element_blank()) + 
  theme(plot.margin = margin(t = 0,  # Top margin
                             r = -.1,  # Right margin
                             b = 0,  # Bottom margin
                             l = -.1))+theme(text = element_text(family = "Arial"))+labs(fill = "Change\nper year")+coord_cartesian(xlim = c(5, 50), expand = T)

tile_del2 <- ggplot(del) + aes_string(x = "data", y = 1,   
                                     fill = "derivative") + geom_raster(interpolate = FALSE)+scale_fill_gradient2(high="#dd4124",low="#00496f", mid="white", limits = c(min(derv$derivative),max(derv$derivative)),guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) + xlab("")+ ylab("")+theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_y_continuous(breaks = seq(1), labels = "22qDel") +
  theme(axis.ticks.y = element_blank(),axis.text.x = element_blank()) +
  theme(plot.margin = margin(t = 0,  # Top margin
                             r = 0,  # Right margin
                             b = 0,  # Bottom margin
                             l = 0))+theme(text = element_text(family = "Arial"))+labs(fill = "Change\nper year")+coord_cartesian(xlim = c(5, 50), expand = T)

tile_dup2 <- ggplot(dup) + aes_string(x = "data", y = 1,   
                                     fill = "derivative") + geom_raster(interpolate = FALSE)+scale_fill_gradient2(high="#dd4124",low="#00496f", mid="white", limits = c(min(derv$derivative),max(derv$derivative)),guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) + xlab("")+ ylab("")+theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +theme(legend.position = "none") + scale_y_continuous(breaks = seq(1), labels = "22qDup") +
  theme(axis.ticks.y = element_blank()) + 
  theme(plot.margin = margin(t = 0,  # Top margin
                             r = 0,  # Right margin
                             b = 0,  # Bottom margin
                             l = 0))+theme(text = element_text(family = "Arial"))+labs(fill = "Change\nper year")+coord_cartesian(xlim = c(5, 50), expand = T)


filename=sprintf('/Users/rahay/Box/22q_dup_del_sMRI/results/SFig_%s_FONT.pdf',var)



plots[[i]] <-  (p2/ (tile_con2/tile_del2/tile_dup2+plot_layout(guides='collect')) + plot_layout(nrow=2,heights=c(6,3))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())) +theme(text = element_text(family = "Arial"))




#ggsave(filename,tiled,width=5,height=4)
#Sys.setenv(R_GSCMD = "/usr/local/bin/gs")
#Sys.setenv(R_GSCMD = "C:\\Program Files\\gs\\gs9.54.0\\bin\\gswin64.exe")
#embed_fonts(filename,format="pdfwrite",options="-sFONTPATH=/Users/rahay/Box/22q_dup_del_sMRI/code/arialfont/arial.ttf")
        
#dev.off()
}

tiledplots<- (plots[[1]] | plots[[2]] | plots[[3]] | plots[[4]])/(plots[[5]]|plots[[6]]|plots[[7]]|plots[[8]]) + plot_annotation(title="Cortical Thickness")
ggsave('/Users/rahay/Box/22q_dup_del_sMRI/results/tiledplots_TEST.pdf',tiledplots,width=20,height=10.5)

tiledplots2<-(plots[[9]] | plots[[10]] | plots[[11]] | plots[[12]])/(plots[[13]]|plots[[14]]|plots[[15]]|plots[[16]]) + plot_annotation(title="Surface Area")
ggsave('/Users/rahay/Box/22q_dup_del_sMRI/results/tiledplots_SA_TEST.pdf',tiledplots2,width=20,height=10.5)
