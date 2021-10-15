
library(readxl)
library(plyr)
library(ggplot2)
library(PNWColors)
library(viridis)
pnw_palette(name="Starfish",n=4,type="discrete")
adj_df<-read_xlsx('/Users/jalbrzikowskime/Box/22q_dup_del_sMRI/data/20210628_aftercombat_data.xlsx')

summary(cars)
#plot visists and & ages
#true figure is in rmd
ages<-adj_df
minAgeLookup <- ddply(ages, .(PTID), function(x){min(x$Age)} )
names(minAgeLookup)[2] <- 'startage'
ages <- merge(minAgeLookup, ages)
#ages2<-ages[order(ages$startage),]
#ages2$count<-
ages$SUBJECT_IDENTITY<-factor(ages$SUBJECT_IDENTITY,levels=c("PATIENT-DEL","CONTROL","PATIENT-DUP"))
p <- ggplot(ages,aes(x=Age,y=as.factor(startage),group=PTID)) +
  geom_line(aes(color=SUBJECT_IDENTITY)) + geom_point(aes(fill=SUBJECT_IDENTITY),shape=21,size=1.6) + theme_bw(base_size=12) +
  labs(x='Age (years) at scan', y='')+ #,title='Longitudinal Recordings') + 
  theme(axis.text.y=element_blank(),  panel.grid.minor.y=element_blank(),panel.grid.major.y=element_blank(), axis.ticks.y=element_blank(),panel.grid = element_blank())+
  scale_x_continuous(breaks=seq(0,51,by=4)) + coord_cartesian(ylim=c(-2,170))+
  theme(panel.border = element_blank(),axis.line = element_line(color="black"))+theme(legend.position = "top")+
  scale_fill_viridis(discrete=TRUE)+scale_color_viridis(discrete = TRUE)
p  
  
