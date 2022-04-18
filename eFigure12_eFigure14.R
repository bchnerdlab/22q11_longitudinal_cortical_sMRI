##make sure you've run the following:
#library(remotes)
#remotes::install_version("Rttf2pt1", version = "1.3.8")


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
library(extrafont)
library(equatiomatic)



font_import(paths = "/Users/jalbrzikowskime/Dropbox (BCH)/22q_dup_del_sMRI/code/arialfont",prompt=T)

adj_df<-read_xlsx('/Users/jalbrzikowskime/Dropbox (BCH)/22q_dup_del_sMRI/data/20210628_aftercombat_data.xlsx')
adj_df$Mean_thickness.combat<-adj_df$lh_MeanThickness_thickness.combat+adj_df$rh_MeanThickness_thickness.combat
adj_df$TotalSA_area.combat<-adj_df$lh_TotalSA_area.combat+adj_df$rh_TotalSA_area.combat
#22q-del ASD specific
del22q<-adj_df[adj_df$SUBJECT_IDENTITY=="PATIENT-DEL"| adj_df$SUBJECT_IDENTITY=="CONTROL",]

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

cont<-which(del22q$SUBJECT_IDENTITY=="CONTROL")
del22q$ASD.DIAGNOS[cont]<-"control"


check2<-table(del22q$PTID,del22q$ASD.DIAGNOS)
del22q<-del22q[del22q$ASD.DIAGNOS!=3,]

del22q$PTID<-as.factor(del22q$PTID)
del22q$ASD.DIAGNOS<-as.factor(del22q$ASD.DIAGNOS)

un<-del22q[!duplicated(del22q$PTID),]


del22q<-del22q[del22q$Age<26,]


brain_regions<-c("lh_frontal_thickness.combat","rh_frontal_thickness.combat","lh_temporal_thickness.combat","rh_temporal_thickness.combat","lh_parietal_thickness.combat","rh_parietal_thickness.combat","lh_occipital_thickness.combat","rh_occipital_thickness.combat","lh_frontal_area.combat","rh_frontal_area.combat",
                 "lh_temporal_area.combat","rh_temporal_area.combat","lh_parietal_area.combat","rh_parietal_area.combat","lh_occipital_area.combat","rh_occipital_area.combat")
plotnames<-c("LH Frontal","RH Frontal","LH Temporal","RH Temporal","LH Parietal","RH Parietal","LH Occipital","RH Occipital","LH Frontal","RH Frontal","LH Temporal","RH Temporal","LH Parietal","RH Parietal","LH Occipital","RH Occipital")
axislab<-c("mm","mm","mm","mm","mm","mm","mm","mm",expression(mm^2),expression(mm^2),expression(mm^2),expression(mm^2),expression(mm^2),expression(mm^2),expression(mm^2),expression(mm^2))



adj_df$Mean_thickness.combat<-adj_df$lh_MeanThickness_thickness.combat+adj_df$rh_MeanThickness_thickness.combat
adj_df$TotalSA_area.combat<-adj_df$lh_TotalSA_area.combat+adj_df$rh_TotalSA_area.combat

plots<-list()

for (i in c(1:length(brain_regions))){
  
  var=brain_regions[i]
  mod<-gam(del22q[[var]]~s(Age,by=ASD.DIAGNOS,k=4)+
             ASD.DIAGNOS+SEX+SCANNER+EstimatedTotalIntraCranialVol.combat+s(PTID,bs="re"),selection=TRUE,method="REML",data=del22q)
  lmod<-lm(del22q[[var]]~ASD.DIAGNOS+SEX+EstimatedTotalIntraCranialVol.combat,data=del22q)
  del22q$yvar.res<-lmod$residuals

  try<-evaluate_smooth(mod,"s(Age)",n=20)
  
  try$selo <- try$est - try$se
  try$sehi <- try$est + try$se

line_size=2
try$ASD.DIAGNOS<-factor(try$ASD.DIAGNOS,levels=c("0","1","control"))
p2<-ggplot(data = try, aes_string(x = "Age",y = "est"))+
  geom_ribbon(data = try,aes_string(x = "Age", ymin = "selo",ymax = "sehi", fill = "ASD.DIAGNOS"),alpha = .18, linetype = 0)+
  geom_line(data = try,aes_string(x = "Age", y = "est",color = "ASD.DIAGNOS"),size = 1)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  theme(legend.title = element_blank())+ylab(axislab[[i]])+xlab("Age (years)")+
  geom_point(data = del22q,aes_string(x = "Age", y = "yvar.res",color="ASD.DIAGNOS"),size = 1,alpha=0.4)+
  geom_line(data = del22q,aes_string(x = "Age", y = "yvar.res",group="PTID",color = "ASD.DIAGNOS"),size = 0.2,alpha=0.4)+
  scale_fill_manual(values=c("#334B4A","#AF4974","#21908C"),labels = c("22qDel-ASD-", "22qDel-ASD+","Control"))+
  scale_color_manual(values=c("#334B4A","#AF4974","#21908C"),labels = c("22qDel-ASD-", "22qDel-ASD+","Control"))+
  theme(plot.margin = margin(t = 0,  r = 0, b = 0, l = 0)) +ggtitle(plotnames[[i]]) +theme(text = element_text(family = "Arial"))+coord_cartesian(xlim = c(6, 25), expand = T)
derv<-derivatives(mod,term=c("s(Age):ASD.DIAGNOS0","s(Age):ASD.DIAGNOS1","s(Age):ASD.DIAGNOScontrol"),n=20)

del22q_noASD=derv[derv$smooth=="s(Age):ASD.DIAGNOS0",]
del22q_ASD=derv[derv$smooth=="s(Age):ASD.DIAGNOS1",]
cont=derv[derv$smooth=="s(Age):ASD.DIAGNOScontrol",]

tile_con2 <- ggplot(cont) + aes_string(x = "data", y = 1,   
                                      fill = "derivative") + geom_raster(interpolate = FALSE)+scale_fill_gradient2(high="#dd4124",low="#00496f", mid="white", limits = c(min(derv$derivative),max(derv$derivative)),guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) + xlab("")+ ylab("")+theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_y_continuous(breaks = seq(1), labels = "Control") +
  theme(legend.position = "none")+ theme(axis.ticks.y = element_blank(),axis.text.x = element_blank()) + 
  theme(plot.margin = margin(t = 0,  # Top margin
                             r = -.1,  # Right margin
                             b = 0,  # Bottom margin
                             l = -.1))+theme(text = element_text(family = "Arial"))+labs(fill = "Chang per\nyear (mm)")+coord_cartesian(xlim = c(6, 25), expand = T)

tile_del22q_noASD <- ggplot(del22q_noASD) + aes_string(x = "data", y = 1,   
                                     fill = "derivative") + geom_raster(interpolate = FALSE)+scale_fill_gradient2(high="#dd4124",low="#00496f", mid="white", limits = c(min(derv$derivative),max(derv$derivative)),guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) + xlab("")+ ylab("")+theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_y_continuous(breaks = seq(1), labels = "22qDel-no ASD") +
  theme(axis.ticks.y = element_blank(),axis.text.x = element_blank()) +
  theme(plot.margin = margin(t = 0,  # Top margin
                             r = 0,  # Right margin
                             b = 0,  # Bottom margin
                             l = 0))+theme(text = element_text(family = "Arial"))+labs(fill = "Chang per\nyear (mm)")+coord_cartesian(xlim = c(6, 25), expand = T)

tile_del22q_ASD <- ggplot(del22q_ASD) + aes_string(x = "data", y = 1,   
                                     fill = "derivative") + geom_raster(interpolate = FALSE)+scale_fill_gradient2(high="#dd4124",low="#00496f", mid="white", limits = c(min(derv$derivative),max(derv$derivative)),guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) + xlab("")+ ylab("")+theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +theme(legend.position = "none") + scale_y_continuous(breaks = seq(1), labels = "22qDel-ASD") +
  theme(axis.ticks.y = element_blank()) + 
  theme(plot.margin = margin(t = 0,  # Top margin
                             r = 0,  # Right margin
                             b = 0,  # Bottom margin
                             l = 0))+theme(text = element_text(family = "Arial"))+labs(fill = "Chang per\nyear (mm)")+coord_cartesian(xlim = c(6, 25), expand = T)


filename=sprintf('/Users/jalbrzikowskime/Dropbox (BCH)/22q_dup_del_sMRI/results/SFig_%s_FONT.pdf',var)



plots[[i]] <-  (p2/ (tile_con2/tile_del22q_noASD/tile_del22q_ASD+plot_layout(guides='collect')) + plot_layout(nrow=2,heights=c(6,3))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())) +theme(text = element_text(family = "Arial"))




#ggsave(filename,tiled,width=5,height=4)
#Sys.setenv(R_GSCMD = "/usr/local/bin/gs")
#Sys.setenv(R_GSCMD = "C:\\Program Files\\gs\\gs9.54.0\\bin\\gswin64.exe")
#embed_fonts(filename,format="pdfwrite",options="-sFONTPATH=/Users/jalbrzikowskime/Dropbox (BCH)/22q_dup_del_sMRI/code/arialfont/arial.ttf")
        
#dev.off()



}

tiledplots<- (plots[[1]] | plots[[2]] | plots[[3]] | plots[[4]])/(plots[[5]]|plots[[6]]|plots[[7]]|plots[[8]]) + plot_annotation(title="Cortical Thickness")
ggsave('/Users/jalbrzikowskime/Dropbox (BCH)/22q_dup_del_sMRI/results/tiledplots_CT_ASD_del.pdf',tiledplots,width=20,height=10.5,device=cairo_pdf)
Sys.setenv(R_GSCMD = "/usr/local/bin/gs")
embed_fonts("/Users/jalbrzikowskime/Dropbox (BCH)/22q_dup_del_sMRI/results/tiledplots_CT_ASD_del.pdf")
dev.off()

tiledplots2<-(plots[[9]] | plots[[10]] | plots[[11]] | plots[[12]])/(plots[[13]]|plots[[14]]|plots[[15]]|plots[[16]]) + plot_annotation(title="Surface Area")
ggsave('/Users/jalbrzikowskime/Dropbox (BCH)/22q_dup_del_sMRI/results/tiledplots_SA_ASD_del_wpoints.pdf',tiledplots2,width=20,height=10.5,device=cairo_pdf)
Sys.setenv(R_GSCMD = "/usr/local/bin/gs")
embed_fonts("/Users/jalbrzikowskime/Dropbox (BCH)/22q_dup_del_sMRI/results/tiledplots_SA_ASD_wpoints.pdf")
dev.off()
