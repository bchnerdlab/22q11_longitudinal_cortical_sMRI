##make sure you've run the following:
#library(remotes)
#remotes::install_version("Rttf2pt1", version = "1.3.8")

library(ggplot2)
library(patchwork)
library(viridis)
library(tidyverse)
library(extrafont)
library(gratia)
library(readxl)
library(mgcv)

adj_df<-read_xlsx("/Users/jalbrzikowskime/Box/22q_dup_del_sMRI/data/20210628_aftercombat_data.xlsx")
adj_df$PTID<-factor(adj_df$PTID)

adj_df$Group <- factor(adj_df$SUBJECT_IDENTITY)
adj_df$PTID<-factor(adj_df$PTID)
adj_df$Mean_thickness.combat<-adj_df$lh_MeanThickness_thickness.combat+adj_df$rh_MeanThickness_thickness.combat
adj_df$TotalSA_area.combat<-adj_df$lh_TotalSA_area.combat+adj_df$rh_TotalSA_area.combat

#22q-del psychosis specific
del22q<-adj_df[adj_df$SUBJECT_IDENTITY=="PATIENT-DEL",]#| adj_df$SUBJECT_IDENTITY=="CONTROL",]
#according to SIPSPSY.code.definition, 3 means "not assessed"
#removed those with 3 and 22qdel
rem<-which(del22q$SIPSPSY.code==3 & del22q$SUBJECT_IDENTITY=="PATIENT-DEL")
del22q<-del22q[-rem,]
#prodromal/hr=1, psychotic =2
#create psychosis spectrum
del22q$psychosis_spect<-gsub(2,1,del22q$SIPSPSY.code)
del22q$psychosis_spect<-gsub(1,"psychosis_spect",del22q$psychosis_spect)
del22q$psychosis_spect<-gsub(0,"nopsychosis",del22q$psychosis_spect)
#del22q$psychosis_spect<-as.factor(del22q$psychosis_spect)
cont<-which(del22q$SUBJECT_IDENTITY=="CONTROL")
del22q$psychosis_spect[cont]<-"control"
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

del22q$psychosis_spect2<-as.factor(del22q$psychosis_spect2)

#del22q<-del22q[del22q$Age<26,]


mod<-gam(del22q$Mean_thickness.combat~s(Age,by=psychosis_spect2,k=4)+psychosis_spect2+SEX+SCANNER+
           EstimatedTotalIntraCranialVol.combat+s(PTID,bs="re"),selection=TRUE,method="REML",data=del22q)
try<-evaluate_smooth(mod,"s(Age)",n=44)

try$selo <- try$est - try$se
try$sehi <- try$est + try$se

line_size=2

test<-difference_smooths(model=mod,"s(Age)")

##This is the font import command that makes the font change work.  Hopefully.
#font_import(paths = "/Users/rahay/Box/22q_dup_del_sMRI/code/arialfont",prompt=T)

p1<-ggplot(data = try, aes_string(x = "Age",y = "est") )+
  geom_ribbon(data = try,aes_string(x = "Age", ymin = "selo",ymax = "sehi", fill = "psychosis_spect2"),alpha = .18, linetype = 0)+scale_fill_manual(values=c("#7570B3","#E7298A"))+geom_line(data = try,aes_string(x = "Age", y = "est",color = "psychosis_spect2"),size = 1)+
  scale_color_manual(values=c("#7570B3","#E7298A"))+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme(legend.title = element_blank())+ylab("mm")+xlab("")+theme(text = element_text(family = "Arial"))
#ggsave("p1.pdf")
#Sys.setenv(R_GSCMD = "/usr/local/bin/gs")
#embed_fonts("p1.pdf")

derv<-derivatives(mod,term="Age",n=44)
derv$t<-derv$derivative/derv$se
del22q_nopsychosis=derv[derv$smooth=="s(Age):psychosis_spect2nopsychosis",]
del22q_yessychosis=derv[derv$smooth=="s(Age):psychosis_spect2psychosis_spect",]
del22q_yessychosis$t<-del22q_yessychosis$derivative/del22q_yessychosis$se
#cont=derv[derv$smooth=="s(Age):psychosis_spect2control",]

tile_del22q_nopsychosis <- ggplot(del22q_nopsychosis) + aes_string(x = "data", y = 1,   
                                                                   fill = "t") + geom_raster(interpolate = FALSE)+scale_fill_gradient2(high="#dd4124",low="#00496f", mid="white", limits = c(min(derv$t),max(derv$t))) + xlab("")+ ylab("22qdel-nopsychosis")+theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +theme(legend.position = "none")+ theme(axis.ticks.y = element_blank(),
                                                                                                                        axis.text.y = element_blank())

tile_del22q_psychosis <- ggplot(del22q_yessychosis) + aes_string(x = "data", y = 1,   
                                                                 fill = "t") + geom_raster(interpolate = FALSE)+scale_fill_gradient2(high="#dd4124",low="#00496f", mid="white",limits = c(min(derv$t),max(derv$t))) + xlab("")+ ylab("22qdel-psychosis")+theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme(axis.ticks.y = element_blank(),axis.text.y = element_blank())+labs(fill="Change \nper year")

#tile_cont <- ggplot(cont) + aes_string(x = "data", y = 1,   
                              #         fill = "derivative") + geom_raster(interpolate = FALSE)+scale_fill_gradient2(high="#dd4124",low="#00496f", mid="white", limits = c(min(derv$derivative),max(derv$derivative))) + xlab("")+ ylab("Control")+theme_bw()+
#  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme(axis.ticks.y = element_blank(),axis.text.y = element_blank())+labs(fill="Change \nper year")



pdf('/Users/jalbrzikowskime/Box/22q_dup_del_sMRI/results/Figure03_age_CT_psychosistest.pdf',width=7,height=8)
print((p1|tile_del22q_nopsychosis|tile_del22q_psychosis)+plot_layout(nrow=4,heights = c(5,1,1,1)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
dev.off()



mod<-gam(del22q$TotalSA_area.combat~s(Age,by=psychosis_spect2,k=4)+psychosis_spect2+SEX+SCANNER+
           EstimatedTotalIntraCranialVol.combat+s(PTID,bs="re"),selection=TRUE,method="REML",data=del22q)

try<-evaluate_smooth(mod,"s(Age)",n=44)



try<-evaluate_smooth(mod,"s(Age)",n=44)

try$selo <- try$est - try$se
try$sehi <- try$est + try$se

line_size=2

p2<-ggplot(data = try, aes_string(x = "Age",y = "est") )+
  geom_ribbon(data = try,aes_string(x = "Age", ymin = "selo",ymax = "sehi", fill = "psychosis_spect2"),alpha = .18, linetype = 0)+scale_fill_manual(values=c("#7570B3","#E7298A"))+geom_line(data = try,aes_string(x = "Age", y = "est",color = "psychosis_spect2"),size = 1)+
  scale_color_manual(values=c("#7570B3","#E7298A"))+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme(legend.title = element_blank())+ylab("mm")+xlab("")
derv<-derivatives(mod,term="Age",n=44)
derv$t<-derv$derivative/derv$se
del22q_nopsychosis=derv[derv$smooth=="s(Age):psychosis_spect2nopsychosis",]
del22q_yessychosis=derv[derv$smooth=="s(Age):psychosis_spect2psychosis_spect",]
del22q_yessychosis$t<-del22q_yessychosis$derivative/del22q_yessychosis$se
#cont=derv[derv$smooth=="s(Age):psychosis_spect2control",]

tile_del22q_nopsychosis2 <- ggplot(del22q_nopsychosis) + aes_string(x = "data", y = 1,   
                                                                   fill = "t") + geom_raster(interpolate = FALSE)+scale_fill_gradient2(high="#dd4124",low="#00496f", mid="white", limits = c(min(derv$t),max(derv$t))) + xlab("")+ ylab("22qdel-nopsychosis")+theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +theme(legend.position = "none")+ theme(axis.ticks.y = element_blank(),
                                                                                                                        axis.text.y = element_blank())

tile_del22q_psychosis2 <- ggplot(del22q_yessychosis) + aes_string(x = "data", y = 1,   
                                                                 fill = "t") + geom_raster(interpolate = FALSE)+scale_fill_gradient2(high="#dd4124",low="#00496f", mid="white",limits = c(min(derv$t),max(derv$t))) + xlab("")+ ylab("22qdel-psychosis")+theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme(axis.ticks.y = element_blank(),axis.text.y = element_blank())+labs(fill="Change \nper year")

#tile_cont2 <- ggplot(cont) + aes_string(x = "data", y = 1,   
                                       fill = "derivative") + geom_raster(interpolate = FALSE)+scale_fill_gradient2(high="#dd4124",low="#00496f", mid="white", limits = c(min(derv$derivative),max(derv$derivative))) + xlab("")+ ylab("Control")+theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme(axis.ticks.y = element_blank(),axis.text.y = element_blank())+labs(fill="Change \nper year")

pdf('/Users/jalbrzikowskime/Box/22q_dup_del_sMRI/results/Figure02_age_SA_psychosistest.pdf',width=7,height=8)
print((p2|tile_del22q_nopsychosis2 |tile_del22q_psychosis2)+plot_layout(nrow=4,heights = c(5,1,1,1)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
dev.off()


un<-del22q[!duplicated(del22q$PTID),]

