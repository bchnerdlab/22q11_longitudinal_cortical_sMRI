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






mod<-gam(del22q$Mean_thickness.combat~s(Age,by=ASD.DIAGNOS,k=4)+
           ASD.DIAGNOS+SEX+SCANNER+EstimatedTotalIntraCranialVol.combat+s(PTID,bs="re"),selection=TRUE,method="REML",data=del22q)
lmod<-lm(del22q$Mean_thickness.combat~ASD.DIAGNOS+SEX+EstimatedTotalIntraCranialVol.combat,data=del22q)
draw(mod,residuals = T,scale="fixed")
del22q$Mean_thickness.combat.res<-lmod$residuals
derv<-derivatives(mod,term=c("s(Age):ASD.DIAGNOS0","s(Age):ASD.DIAGNOS1","s(Age):ASD.DIAGNOScontrol"),n=20)
noASD=derv[derv$smooth=="s(Age):ASD.DIAGNOS0",]


del=derv[derv$smooth=="s(Age):ASD.DIAGNOS1",]

try<-evaluate_smooth(mod,"s(Age)")

try$selo <- try$est - try$se
try$sehi <- try$est + try$se

line_size=2

derv$t<-derv$derivative/derv$se

# old legend position is theme(legend.position = c(0.65, 0.8) or 0.35,0.25
p1<-ggplot(data = try, aes_string(x = "Age",y = "est") )+
  geom_ribbon(data = try,aes_string(x = "Age", ymin = "selo",ymax = "sehi", fill = "ASD.DIAGNOS"),alpha = .18, linetype = 0)+
  scale_fill_manual(values=c("#334B4A","#AF4974","#21908C"),labels = c("22qdel-no ASD (N=35)", "22qdel-ASD (N=40)","Controls (N=114)"))+
  geom_line(data = try,aes_string(x = "Age", y = "est",color = "ASD.DIAGNOS"),size = 2)+
  scale_color_manual(values=c("#334B4A","#AF4974","#21908C"),labels = c("22qdel-no ASD (N=35)", "22qdel-ASD (N=40)","Controls (N=114)"))+theme_bw()+
  geom_point(data = del22q,aes_string(x = "Age", y = "Mean_thickness.combat.res",color="ASD.DIAGNOS"),size = 1,alpha=0.4)+
  geom_line(data = del22q,aes_string(x = "Age", y = "Mean_thickness.combat.res",group="PTID",color = "ASD.DIAGNOS"),size = 0.2,alpha=0.4)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  theme(legend.title = element_blank(),axis.text.x = element_blank())+ylab("mm")+
  theme(text = element_text(family = "ArialMT")) +theme(legend.position = c(0.275, 0.175),legend.background = element_rect(fill='transparent') )+
  ggtitle("Mean Cortical Thickness")+xlab("")+coord_cartesian(xlim = c(5.5, 25.25), expand = T)
del22q_noASD=derv[derv$smooth=="s(Age):ASD.DIAGNOS0",]
#del22q_noASD$t<-del22q_noASD$derivative/del22q_noASD$se
del22q_ASD=derv[derv$smooth=="s(Age):ASD.DIAGNOS1",]
cont=derv[derv$smooth=="s(Age):ASD.DIAGNOScontrol",]
tile_cont <- ggplot(cont) + aes_string(x = "data", y = 1,   
                                       fill = "derivative") + geom_raster(interpolate = FALSE)+scale_fill_gradient2(high="#dd4124",low="#00496f", mid="white", limits = c(min(derv$derivative),max(derv$derivative)),guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) + xlab("")+ ylab("")+theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_y_continuous(breaks = seq(1), labels = "Control") +
  theme(legend.position = "none")+ theme(axis.ticks.y = element_blank(),axis.text.x = element_blank()) + 
  theme(plot.margin = margin(t = 0,  # Top margin
                             r = -.1,  # Right margin
                             b = 0,  # Bottom margin
                             l = -.1))+theme(text = element_text(family = "Arial"))+labs(fill = "Change per\nyear (mm)")+coord_cartesian(xlim = c(5.5, 25.25), expand = T)

tile_del22q_noASD <- ggplot(del22q_noASD) + aes_string(x = "data", y = 1,   
                                                       fill = "derivative") + geom_raster(interpolate = FALSE)+scale_fill_gradient2(high="#dd4124",low="#00496f", mid="white", limits = c(min(derv$derivative),max(derv$derivative)),guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) + xlab("")+ ylab("")+theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_y_continuous(breaks = seq(1), labels = "22qdel-\nno ASD") +
  theme(axis.ticks.y = element_blank(),axis.text.x = element_blank()) + 
  theme(plot.margin = margin(t = 0,  # Top margin
                             r = 0,  # Right margin
                             b = 0,  # Bottom margin
                             l = 0))+theme(text = element_text(family = "Arial"))+labs(fill = "Change per\nyear (mm)")+coord_cartesian(xlim = c(5.5, 25.25), expand = T)

tile_del22q_ASD <- ggplot(del22q_ASD) + aes_string(x = "data", y = 1,   
                                                   fill = "derivative") + geom_raster(interpolate = FALSE)+scale_fill_gradient2(high="#dd4124",low="#00496f", mid="white", limits = c(min(derv$derivative),max(derv$derivative)),guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) + xlab("Age (years)")+ ylab("")+theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_y_continuous(breaks = seq(1), labels = "22qdel-\nASD") +
  theme(legend.position = "none")+ theme(axis.ticks.y = element_blank()) + 
  theme(plot.margin = margin(t = 0,  # Top margin
                             r = -5,  # Right margin
                             b = 0,  # Bottom margin
                             l = -5))+theme(text = element_text(family = "Arial"))+labs(fill = "Change per\nyear (mm)")+coord_cartesian(xlim = c(5.5, 25.25), expand = T)


fig4a_left<-(p1/ (tile_cont/tile_del22q_noASD/tile_del22q_ASD+plot_layout(guides='collect')) + plot_layout(nrow=2,heights=c(5,3))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())) +plot_annotation(title="A.")+theme(text = element_text(family = "Arial"))
#ggsave("/Users/jalbrzikowskime/Dropbox (BCH)/22q_dup_del_sMRI/results/Fig4A_left.pdf",fig4b_left,width = 7, height = 8)
#Sys.setenv(R_GSCMD = "/usr/local/bin/gs")
#embed_fonts("/Users/jalbrzikowskime/Dropbox (BCH)/22q_dup_del_sMRI/results/Fig4A_left.pdf")




mod<-gam(del22q$TotalSA_area.combat~s(Age,by=ASD.DIAGNOS,k=4)+
           ASD.DIAGNOS+SEX+SCANNER+EstimatedTotalIntraCranialVol.combat+s(PTID,bs="re"),selection=TRUE,method="REML",data=del22q)
lmod<-lm(del22q$TotalSA_area.combat~ASD.DIAGNOS+SEX+EstimatedTotalIntraCranialVol.combat,data=del22q)
del22q$TotalSA_area.combat.res<-lmod$residuals
derv<-derivatives(mod,term=c("s(Age):ASD.DIAGNOS0","s(Age):ASD.DIAGNOS1","s(Age):ASD.DIAGNOScontrol"),n=20)
noASD=derv[derv$smooth=="s(Age):ASD.DIAGNOS0",]


del=derv[derv$smooth=="s(Age):ASD.DIAGNOS1",]

try<-evaluate_smooth(mod,"s(Age)")

try$selo <- try$est - try$se
try$sehi <- try$est + try$se

line_size=2


p1<-ggplot(data = try, aes_string(x = "Age",y = "est") )+
  geom_ribbon(data = try,aes_string(x = "Age", ymin = "selo",ymax = "sehi", fill = "ASD.DIAGNOS"),alpha = .18, linetype = 0)+
  scale_fill_manual(values=c("#334B4A","#AF4974","#21908C"),labels = c("22qdel-no ASD (N=35)", "22qdel-ASD (N=40)","Controls (N=114)"))+
  geom_line(data = try,aes_string(x = "Age", y = "est",color = "ASD.DIAGNOS"),size = 2)+
  scale_color_manual(values=c("#334B4A","#AF4974","#21908C"),labels = c("22qdel-no ASD (N=35)", "22qdel-ASD (N=40)","Controls (N=114)"))+theme_bw()+
  geom_point(data = del22q,aes_string(x = "Age", y = "TotalSA_area.combat.res",color="ASD.DIAGNOS"),size = 1,alpha=0.4)+
  geom_line(data = del22q,aes_string(x = "Age", y = "TotalSA_area.combat.res",group="PTID",color = "ASD.DIAGNOS"),size = 0.2,alpha=0.4)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  theme(legend.title = element_blank(),axis.text.x = element_blank())+ylab(expression(mm^2))+
  theme(text = element_text(family = "Arial")) +theme(legend.position = c( 0.35, 0.175),legend.background = element_rect(fill='transparent'))+
  ggtitle("Total Surface Area")+xlab("")+coord_cartesian(xlim = c(5.5, 25.25), expand = T)
del22q_noASD=derv[derv$smooth=="s(Age):ASD.DIAGNOS0",]
del22q_noASD$t<-del22q_noASD$derivative/del22q_noASD$se
del22q_ASD=derv[derv$smooth=="s(Age):ASD.DIAGNOS1",]
cont=derv[derv$smooth=="s(Age):ASD.DIAGNOScontrol",]
mm2<-bquote(atop(Change~per~phantom(),year~(mm^2)))
tile_cont <- ggplot(cont) + aes_string(x = "data", y = 1,   
                                       fill = "derivative") + geom_raster(interpolate = FALSE)+scale_fill_gradient2(high="#dd4124",low="#00496f", mid="white", limits = c(min(derv$derivative),max(derv$derivative)),guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) + xlab("")+ ylab("")+theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_y_continuous(breaks = seq(1), labels = "Control") +
  theme(legend.position = "none")+ theme(axis.ticks.y = element_blank(),axis.text.x = element_blank()) + 
  theme(plot.margin = margin(t = 0,  # Top margin
                             r = -.1,  # Right margin
                             b = 0,  # Bottom margin
                             l = -.1))+theme(text = element_text(family = "Arial"))+labs(fill = mm2)+coord_cartesian(xlim = c(5.5, 25.25), expand = T)

tile_del22q_noASD <- ggplot(del22q_noASD) + aes_string(x = "data", y = 1,   
                                                       fill = "derivative") + geom_raster(interpolate = FALSE)+scale_fill_gradient2(high="#dd4124",low="#00496f", mid="white", limits = c(min(derv$derivative),max(derv$derivative)),guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) + xlab("")+ ylab("")+theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_y_continuous(breaks = seq(1), labels = "22qdel-\nno ASD") +
  theme(axis.ticks.y = element_blank(),axis.text.x = element_blank()) + 
  theme(plot.margin = margin(t = 0,  # Top margin
                             r = 0,  # Right margin
                             b = 0,  # Bottom margin
                             l = 0))+theme(text = element_text(family = "Arial"))+labs(fill = mm2)+coord_cartesian(xlim = c(5.5, 25.25), expand = T)

tile_del22q_ASD <- ggplot(del22q_ASD) + aes_string(x = "data", y = 1,   
                                                   fill = "derivative") + geom_raster(interpolate = FALSE)+scale_fill_gradient2(high="#dd4124",low="#00496f", mid="white", limits = c(min(derv$derivative),max(derv$derivative)),guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) + xlab("Age (years)")+ ylab("")+theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_y_continuous(breaks = seq(1), labels = "22qdel-\nASD") +
  theme(legend.position = "none")+ theme(axis.ticks.y = element_blank()) + 
  theme(plot.margin = margin(t = 0,  # Top margin
                             r = -5,  # Right margin
                             b = 0,  # Bottom margin
                             l = -5))+theme(text = element_text(family = "Arial"))+labs(fill = mm2)+coord_cartesian(xlim = c(5.5, 25.25), expand = T)


fig4a_right<-(p1/ (tile_cont/tile_del22q_noASD/tile_del22q_ASD+plot_layout(guides='collect')) + plot_layout(nrow=2,heights=c(5,3))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())) +plot_annotation(title="B.")+theme(text = element_text(family = "Arial"))






#ggsave("/Users/jalbrzikowskime/Dropbox (BCH)/22q_dup_del_sMRI/results/Fig4A_right.pdf",fig4b_left,width = 7, height = 8)
#Sys.setenv(R_GSCMD = "/usr/local/bin/gs")
#embed_fonts("/Users/jalbrzikowskime/Dropbox (BCH)/22q_dup_del_sMRI/results/Fig4A_right.pdf")#22q-dup ASD specific



dup22q<-adj_df[adj_df$SUBJECT_IDENTITY=="PATIENT-DUP"| adj_df$SUBJECT_IDENTITY=="CONTROL",]
table(dup22q$PTID,dup22q$ASD.DIAGNOS)
#q_0402 had a 3 for not assessed for ASD, need to remove
rem=which(dup22q$PTID=="q_0402")
dup22q=dup22q[-rem,]
dup22q<-dup22q[dup22q$Age<26,]

cont<-which(dup22q$SUBJECT_IDENTITY=="CONTROL")
dup22q$ASD.DIAGNOS[cont]<-"control"



dup22q$PTID<-as.factor(dup22q$PTID)
dup22q$ASD.DIAGNOS<-as.factor(dup22q$ASD.DIAGNOS)

un<-dup22q[!duplicated(dup22q$PTID),]

mod<-gam(dup22q$Mean_thickness.combat~s(Age,by=ASD.DIAGNOS,k=4)+
ASD.DIAGNOS+SEX+SCANNER+EstimatedTotalIntraCranialVol.combat+s(PTID,bs="re"),selection=TRUE,method="REML",data=dup22q)
summary(mod,freq=T)
lmod<-lm(dup22q$Mean_thickness.combat~ASD.DIAGNOS+SEX+EstimatedTotalIntraCranialVol.combat,data=dup22q)
dup22q$Mean_thickness.combat.res<-lmod$residuals
 
derv<-derivatives(mod,term=c("s(Age):ASD.DIAGNOS0","s(Age):ASD.DIAGNOS1","s(Age):ASD.DIAGNOScontrol"),n=20)
    noASD=derv[derv$smooth=="s(Age):ASD.DIAGNOS0",]
    
    dup=derv[derv$smooth=="s(Age):ASD.DIAGNOS1",]
    try<-evaluate_smooth(mod,"s(Age)")
    try$selo <- try$est - try$se
    try$sehi <- try$est + try$se
    
    line_size=2
 
    p1<-ggplot(data = try, aes_string(x = "Age",y = "est") )+
      geom_ribbon(data = try,aes_string(x = "Age", ymin = "selo",ymax = "sehi", fill = "ASD.DIAGNOS"),alpha = .18, linetype = 0)+
      scale_fill_manual(values=c("#6C80B6","#FF89ED","#21908C"),labels = c("22qdup-no ASD (N=12)", "22qdup-ASD (N=14)","Controls (N=114)"))+
      geom_line(data = try,aes_string(x = "Age", y = "est",color = "ASD.DIAGNOS"),size = 2)+
      geom_point(data = dup22q,aes_string(x = "Age", y = "Mean_thickness.combat.res",color="ASD.DIAGNOS"),size = 1,alpha=0.4)+
      geom_line(data = dup22q,aes_string(x = "Age", y = "Mean_thickness.combat.res",group="PTID",color = "ASD.DIAGNOS"),size = 0.2,alpha=0.4)+
      scale_color_manual(values=c("#6C80B6","#FF89ED","#21908C"),labels = c("22qdup-no ASD (N=12)", "22qdup-ASD (N=14)","Controls (N=114)"))+
      theme_bw()+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
      theme(legend.title = element_blank(),axis.text.x = element_blank())+ylab("mm")+
      theme(text = element_text(family = "Arial")) +theme(legend.position = c(0.275, 0.175),legend.background = element_rect(fill='transparent') )+
      ggtitle("Mean Cortical Thickness")+xlab("")+coord_cartesian(xlim = c(5.5, 25.25), expand = T)
    dup22q_noASD=derv[derv$smooth=="s(Age):ASD.DIAGNOS0",]
    dup22q_noASD$t<-dup22q_noASD$derivative/dup22q_noASD$se
   dup22q_ASD=derv[derv$smooth=="s(Age):ASD.DIAGNOS1",]
    cont=derv[derv$smooth=="s(Age):ASD.DIAGNOScontrol",]
    tile_cont <- ggplot(cont) + aes_string(x = "data", y = 1,   
                                           fill = "derivative") + geom_raster(interpolate = FALSE)+scale_fill_gradient2(high="#dd4124",low="#00496f", mid="white", limits = c(min(derv$derivative),max(derv$derivative)),guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) + xlab("")+ ylab("")+theme_bw()+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_y_continuous(breaks = seq(1), labels = "Control") +
      theme(legend.position = "none")+ theme(axis.ticks.y = element_blank(),axis.text.x = element_blank()) + 
      theme(plot.margin = margin(t = 0,  # Top margin
                                 r = -.1,  # Right margin
                                 b = 0,  # Bottom margin
                                 l = -.1))+theme(text = element_text(family = "Arial"))+labs(fill = "Change per\nyear(mm)")+coord_cartesian(xlim = c(5.5, 25.25), expand = T)
    
    
    tile_dup22q_noASD <- ggplot(dup22q_noASD) + aes_string(x = "data", y = 1,   
                                                           fill = "derivative") + geom_raster(interpolate = FALSE)+scale_fill_gradient2(high="#dd4124",low="#00496f", mid="white", limits = c(min(derv$derivative),max(derv$derivative)),guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) + xlab("")+ ylab("")+theme_bw()+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_y_continuous(breaks = seq(1), labels = "22qdup-\nno ASD") +
     theme(axis.ticks.y = element_blank(),axis.text.x = element_blank()) + 
      theme(plot.margin = margin(t = 0,  # Top margin
                                 r = 0,  # Right margin
                                 b = 0,  # Bottom margin
                                 l = 0))+theme(text = element_text(family = "Arial"))+labs(fill = "Change per\nyear (mm)")+coord_cartesian(xlim = c(5.5, 25.25), expand = T)
    
    
    tile_dup22q_ASD <- ggplot(dup22q_ASD) + aes_string(x = "data", y = 1,   
                                                       fill = "derivative") + geom_raster(interpolate = FALSE)+scale_fill_gradient2(high="#dd4124",low="#00496f", mid="white", limits = c(min(derv$derivative),max(derv$derivative)),guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) + xlab("Age (years)")+ ylab("")+theme_bw()+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_y_continuous(breaks = seq(1), labels = "22qdup-\nASD") +
      theme(legend.position = "none")+ theme(axis.ticks.y = element_blank()) + 
      theme(plot.margin = margin(t = 0,  # Top margin
                                 r = -5,  # Right margin
                                 b = 0,  # Bottom margin
                                 l = -5))+theme(text = element_text(family = "Arial"))+labs(fill = "Change per\nyear (mm)")+coord_cartesian(xlim = c(5.5, 25.25), expand = T)
        
   fig4b_left<-(p1/(tile_cont/tile_dup22q_noASD/tile_dup22q_ASD+plot_layout(guides='collect'))+plot_layout(nrow=2,heights = c(5,3))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))+plot_annotation(title="C.") +theme(text = element_text(family = "Arial"))
   
#  ggsave("/Users/jalbrzikowskime/Dropbox (BCH)/22q_dup_del_sMRI/results/Fig4B_left.pdf",fig4b_left,width = 7, height = 8)
#  Sys.setenv(R_GSCMD = "/usr/local/bin/gs")
#  embed_fonts("/Users/jalbrzikowskime/Dropbox (BCH)/22q_dup_del_sMRI/results/Fig4B_left.pdf")

#dev.off()


mod<-gam(dup22q$TotalSA_area.combat~s(Age,by=ASD.DIAGNOS,k=4)+
           ASD.DIAGNOS+SEX+SCANNER+EstimatedTotalIntraCranialVol.combat+s(PTID,bs="re"),selection=TRUE,method="REML",data=dup22q)
summary(mod,freq=T)
lmod<-lm(dup22q$TotalSA_area.combat~ASD.DIAGNOS+SEX+EstimatedTotalIntraCranialVol.combat,data=dup22q)
dup22q$TotalSA_area.combat.res<-lmod$residuals

derv<-derivatives(mod,term=c("s(Age):ASD.DIAGNOS0","s(Age):ASD.DIAGNOS1","s(Age):ASD.DIAGNOScontrol"),n=20)
noASD=derv[derv$smooth=="s(Age):ASD.DIAGNOS0",]

dup=derv[derv$smooth=="s(Age):ASD.DIAGNOS1",]
try<-evaluate_smooth(mod,"s(Age)")
try$selo <- try$est - try$se
try$sehi <- try$est + try$se

line_size=2

p1<-ggplot(data = try, aes_string(x = "Age",y = "est") )+
  geom_ribbon(data = try,aes_string(x = "Age", ymin = "selo",ymax = "sehi", fill = "ASD.DIAGNOS"),alpha = .18, linetype = 0)+
  scale_fill_manual(values=c("#6C80B6","#FF89ED","#21908C"),labels = c("22qdup-no ASD (N=12)", "22qdup-ASD (N=14)","Controls (N=114)"))+
  geom_line(data = try,aes_string(x = "Age", y = "est",color = "ASD.DIAGNOS"),size = 2)+
  geom_point(data = dup22q,aes_string(x = "Age", y = "TotalSA_area.combat.res",color="ASD.DIAGNOS"),size = 1,alpha=0.4)+
  geom_line(data = dup22q,aes_string(x = "Age", y = "TotalSA_area.combat.res",group="PTID",color = "ASD.DIAGNOS"),size = 0.2,alpha=0.4)+
  scale_color_manual(values=c("#6C80B6","#FF89ED","#21908C"),labels = c("22qdup-no ASD (N=12)", "22qdup-ASD (N=14)","Controls (N=114)"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  theme(legend.title = element_blank(),axis.text.x = element_blank())+ylab(expression(mm^2))+
  theme(text = element_text(family = "Arial")) +theme(legend.position = c(0.35, 0.175),legend.background = element_rect(fill='transparent'))+
  ggtitle("Total Surface Area") + xlab("")+coord_cartesian(xlim = c(5.5, 25.25),ylim=c(-40000,20000), expand = T)

dup22q_noASD=derv[derv$smooth=="s(Age):ASD.DIAGNOS0",]
dup22q_noASD$t<-dup22q_noASD$derivative/dup22q_noASD$se
dup22q_ASD=derv[derv$smooth=="s(Age):ASD.DIAGNOS1",]
cont=derv[derv$smooth=="s(Age):ASD.DIAGNOScontrol",]
cont$t<-cont$derivative/cont$se
tile_cont <- ggplot(cont) + aes_string(x = "data", y = 1,   
                                       fill = "derivative") + geom_raster(interpolate = FALSE)+scale_fill_gradient2(high="#dd4124",low="#00496f", mid="white", limits = c(min(derv$derivative),max(derv$derivative)),guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) + xlab("")+ ylab("")+theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_y_continuous(breaks = seq(1), labels = "Control") +
  theme(legend.position = "none")+ theme(axis.ticks.y = element_blank(),axis.text.x = element_blank()) + 
  theme(plot.margin = margin(t = 0,  # Top margin
                             r = -.1,  # Right margin
                             b = 0,  # Bottom margin
                             l = -.1))+theme(text = element_text(family = "Arial"))+labs(fill = mm2)+coord_cartesian(xlim = c(5.5, 25.25), expand = T)
   

tile_dup22q_noASD <- ggplot(dup22q_noASD) + aes_string(x = "data", y = 1,   
                                                                      fill = "derivative") + geom_raster(interpolate = FALSE)+scale_fill_gradient2(high="#dd4124",low="#00496f", mid="white", limits = c(min(derv$derivative),max(derv$derivative)),guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) + xlab("")+ ylab("")+theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_y_continuous(breaks = seq(1), labels = "22qdup-\nno ASD") +
  theme(axis.ticks.y = element_blank(),axis.text.x = element_blank()) + 
  theme(plot.margin = margin(t = 0,  # Top margin
                             r = 0,  # Right margin
                             b = 0,  # Bottom margin
                             l = 0))+theme(text = element_text(family = "Arial"))+labs(fill = mm2)+coord_cartesian(xlim = c(5.5, 25.25), expand = T)


tile_dup22q_ASD <- ggplot(dup22q_ASD) + aes_string(x = "data", y = 1,   
                                                   fill = "derivative") + geom_raster(interpolate = FALSE)+scale_fill_gradient2(high="#dd4124",low="#00496f", mid="white", limits = c(min(derv$derivative),max(derv$derivative)),guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) + xlab("Age (years)")+ ylab("")+theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_y_continuous(breaks = seq(1), labels = "22qdup-\nASD") +
  theme(legend.position = "none")+ theme(axis.ticks.y = element_blank()) + 
  theme(plot.margin = margin(t = 0,  # Top margin
                             r = -5,  # Right margin
                             b = 0,  # Bottom margin
                             l = -5))+theme(text = element_text(family = "Arial"))+labs(fill = mm2)+coord_cartesian(xlim = c(5.5, 25.25), expand = T)



fig4b_right<-(p1/(tile_cont/tile_dup22q_noASD/tile_dup22q_ASD+plot_layout(guides = 'collect'))+plot_layout(nrow=2,heights = c(5,3))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))+plot_annotation(title="D.") +theme(text = element_text(family = "Arial"))
#ggsave("/Users/jalbrzikowskime/Dropbox (BCH)/22q_dup_del_sMRI/results/Fig4B_right.pdf",fig4b_left,width = 7, height = 8)
#Sys.setenv(R_GSCMD = "/usr/local/bin/gs")
#embed_fonts("/Users/jalbrzikowskime/Dropbox (BCH)/22q_dup_del_sMRI/results/Fig4B_right.pdf")

#dev.off()

tiled<-(wrap_elements(fig4a_left)|wrap_elements(fig4a_right))/(wrap_elements(fig4b_left)|wrap_elements(fig4b_right))
tiled

ggsave('/Users/jalbrzikowskime/Dropbox (BCH)/22q_dup_del_sMRI/results/fig4_wpoints.pdf',tiled,width=11,height=11,device = cairo_pdf)
        