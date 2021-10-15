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



font_import(paths = "/Users/jalbrzikowskime/Box/22q_dup_del_sMRI/code/arialfont",prompt=T)


adj_df<-read_xlsx("/Users/jalbrzikowskime/Box/22q_dup_del_sMRI/data/20210628_aftercombat_data.xlsx")
adj_df$PTID<-factor(adj_df$PTID)

#22q-del psychosis specific
del22q<-adj_df[adj_df$SUBJECT_IDENTITY=="PATIENT-DEL"| adj_df$SUBJECT_IDENTITY=="CONTROL",]
#according to SIPSPSY.code.definition, 3 means "not assessed"
#removed those with 3 and 22qdel
rem<-which(del22q$SIPSPSY.code==3 & del22q$SUBJECT_IDENTITY=="PATIENT-DEL")
del22q<-del22q[-rem,]
#prodromal/hr=1, psychotic =2
#create psychosis spectrum
del22q$psychosis_spect<-gsub(2,1,del22q$SIPSPSY.code)
del22q$psychosis_spect<-gsub(1,"psychosis_spect",del22q$psychosis_spect)
del22q$psychosis_spect<-gsub(0,"nopsychosis",del22q$psychosis_spect)
cont<-which(del22q$SUBJECT_IDENTITY=="CONTROL")
del22q$psychosis_spect[cont]<-"control"
del22q$psychosis_spect<-as.factor(del22q$psychosis_spect)




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








brain_regions<-c("lh_frontal_thickness.combat","rh_frontal_thickness.combat","lh_temporal_thickness.combat","rh_temporal_thickness.combat","lh_parietal_thickness.combat","rh_parietal_thickness.combat","lh_occipital_thickness.combat","rh_occipital_thickness.combat","lh_frontal_area.combat","rh_frontal_area.combat",
                 "lh_temporal_area.combat","rh_temporal_area.combat","lh_parietal_area.combat","rh_parietal_area.combat","lh_occipital_area.combat","rh_occipital_area.combat")
plotnames<-c("LH Frontal","RH Frontal","LH Temporal","RH Temporal","LH Parietal","RH Parietal","LH Occipital","RH Occipital","LH Frontal","RH Frontal","LH Temporal","RH Temporal","LH Parietal","RH Parietal","LH Occipital","RH Occipital")
axislab<-c("mm","mm","mm","mm","mm","mm","mm","mm",expression(mm^2),expression(mm^2),expression(mm^2),expression(mm^2),expression(mm^2),expression(mm^2),expression(mm^2),expression(mm^2))


adj_df$psychosis_spect2 <- factor(adj_df$SUBJECT_IDENTITY)
adj_df$PTID<-factor(adj_df$PTID)
adj_df$Mean_thickness.combat<-adj_df$lh_MeanThickness_thickness.combat+adj_df$rh_MeanThickness_thickness.combat
adj_df$TotalSA_area.combat<-adj_df$lh_TotalSA_area.combat+adj_df$rh_TotalSA_area.combat

plots<-list()

for (i in c(1:length(brain_regions))){
  
  var=brain_regions[i]
  mod<-gam(del22q[[var]]~s(Age,by=psychosis_spect2,k=4)+psychosis_spect2+SEX+SCANNER+EstimatedTotalIntraCranialVol.combat+s(PTID,bs="re"),selection=TRUE,method="REML",data=del22q)


try<-evaluate_smooth(mod,"s(Age)",n=44)

try$selo <- try$est - try$se
try$sehi <- try$est + try$se

line_size=2
try$psychosis_spect2<-factor(try$psychosis_spect2,levels=c("control","nopsychosis","psychosis_spect"))
p2<-ggplot(data = try, aes_string(x = "Age",y = "est"))+
  geom_ribbon(data = try,aes_string(x = "Age", ymin = "selo",ymax = "sehi", fill = "psychosis_spect2"),alpha = .18, linetype = 0)+
  geom_line(data = try,aes_string(x = "Age", y = "est",color = "psychosis_spect2"),size = 1)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  theme(legend.title = element_blank())+ylab(axislab[[i]])+xlab("")+
  scale_fill_manual(values=c("#1B9E77","#7570B3","#E7298A"))+scale_color_manual(values=c("#1B9E77","#7570B3","#E7298A"))+
  theme(legend.position = "None") +theme(plot.margin = margin(t = 0,  r = 0, b = 0, l = 0)) +ggtitle(plotnames[[i]]) +theme(text = element_text(family = "Arial"))+coord_cartesian(xlim = c(5, 50), expand = T)
derv<-derivatives(mod,term="Age",n=44,partial_match = TRUE)
#derv1<-derivatives(mod,term=c("s(Age):psychosis_spect2control","s(Age):psychosis_spect2nopsychosis","s(Age):psychosis_spect2psychosis_spect"),n=44)


derv$t<-derv$derivative/derv$se
cont=derv[derv$smooth=="s(Age):psychosis_spect2control",]
cont$t<-cont$derivative/cont$se
del=derv[derv$smooth=="s(Age):psychosis_spect2nopsychosis",]
del$t<-del$derivative/del$se

dup=derv[derv$smooth=="s(Age):psychosis_spect2psychosis_spect",]
dup$t<-dup$derivative/dup$se

tile_con2 <- ggplot(cont) + aes_string(x = "data", y = 1,   
                                      fill = "derivative") + geom_raster(interpolate = FALSE)+scale_fill_gradient2(high="#dd4124",low="#00496f", mid="white", limits = c(min(derv$derivative),max(derv$derivative)),guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) + xlab("")+ ylab("")+theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_y_continuous(breaks = seq(1), labels = "Control") +
  theme(legend.position = "none")+ theme(axis.ticks.y = element_blank(),axis.text.x = element_blank()) + 
  theme(plot.margin = margin(t = 0,  # Top margin
                             r = -.1,  # Right margin
                             b = 0,  # Bottom margin
                             l = -.1))+theme(text = element_text(family = "Arial"))+labs(fill = "Change per\n year (mm)")+coord_cartesian(xlim = c(5, 50), expand = T)

tile_del2 <- ggplot(del) + aes_string(x = "data", y = 1,   
                                     fill = "derivative") + geom_raster(interpolate = FALSE)+scale_fill_gradient2(high="#dd4124",low="#00496f", mid="white", limits = c(min(derv$derivative),max(derv$derivative)),guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) + xlab("")+ ylab("")+theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_y_continuous(breaks = seq(1), labels = "22qDel-PS-") +
  theme(axis.ticks.y = element_blank(),axis.text.x = element_blank()) +
  theme(plot.margin = margin(t = 0,  # Top margin
                             r = 0,  # Right margin
                             b = 0,  # Bottom margin
                             l = 0))+theme(text = element_text(family = "Arial"))+labs(fill = "Change per\n year (mm)")+coord_cartesian(xlim = c(5, 50), expand = T)

tile_dup2 <- ggplot(dup) + aes_string(x = "data", y = 1,   
                                     fill = "derivative") + geom_raster(interpolate = FALSE)+scale_fill_gradient2(high="#dd4124",low="#00496f", mid="white", limits = c(min(derv$derivative),max(derv$derivative)),guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) + xlab("")+ ylab("")+theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +theme(legend.position = "none") + scale_y_continuous(breaks = seq(1), labels = "22qDel-PS+") +
  theme(axis.ticks.y = element_blank()) + 
  theme(plot.margin = margin(t = 0,  # Top margin
                             r = 0,  # Right margin
                             b = 0,  # Bottom margin
                             l = 0))+theme(text = element_text(family = "Arial"))+labs(fill = "Change per\n year (mm)")+coord_cartesian(xlim = c(5, 50), expand = T)


filename=sprintf('/Users/jalbrzikowskime/Box/22q_dup_del_sMRI/results/SFig_%s_FONT.pdf',var)



plots[[i]] <-  (p2/ (tile_con2/tile_del2/tile_dup2+plot_layout(guides='collect')) + plot_layout(nrow=2,heights=c(6,3))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())) +theme(text = element_text(family = "Arial"))




#ggsave(filename,tiled,width=5,height=4)
#Sys.setenv(R_GSCMD = "/usr/local/bin/gs")
#Sys.setenv(R_GSCMD = "C:\\Program Files\\gs\\gs9.54.0\\bin\\gswin64.exe")
#embed_fonts(filename,format="pdfwrite",options="-sFONTPATH=/Users/jalbrzikowskime/Box/22q_dup_del_sMRI/code/arialfont/arial.ttf")
        
#dev.off()



}

tiledplots<- (plots[[1]] | plots[[2]] | plots[[3]] | plots[[4]])/(plots[[5]]|plots[[6]]|plots[[7]]|plots[[8]]) + plot_annotation(title="Cortical Thickness")
ggsave('/Users/jalbrzikowskime/Box/22q_dup_del_sMRI/results/tiledplots_CT_psychosis.pdf',tiledplots,width=20,height=10.5)
Sys.setenv(R_GSCMD = "/usr/local/bin/gs")
embed_fonts("/Users/jalbrzikowskime/Box/22q_dup_del_sMRI/results/tiledplots_CT_psychosis.pdf")
dev.off()

tiledplots2<-(plots[[9]] | plots[[10]] | plots[[11]] | plots[[12]])/(plots[[13]]|plots[[14]]|plots[[15]]|plots[[16]]) + plot_annotation(title="Surface Area")
ggsave('/Users/jalbrzikowskime/Box/22q_dup_del_sMRI/results/tiledplots_SA_psychosis.pdf',tiledplots2,width=20,height=10.5)
Sys.setenv(R_GSCMD = "/usr/local/bin/gs")
embed_fonts("/Users/jalbrzikowskime/Box/22q_dup_del_sMRI/results/tiledplots_SA_psychosis.pdf")
dev.off()