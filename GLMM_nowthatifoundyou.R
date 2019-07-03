library("bbmle") 
library(glmmTMB)
setwd("~/Dropbox/Madison/Biotron")
newspores<-read.csv("NEWSPORES_FINAL.csv", head=T, stringsAsFactors = T)

newspores$RH<-as.factor(newspores$RH)
newspores$temp<-as.factor(newspores$temp)
newspores$lambda<-as.factor(newspores$lambda)
newspores$experiment<-as.factor(newspores$experiment)
newspores$cart<-as.factor(newspores$cart)
newspores$base<-as.factor(newspores$base)
newspores$flux<-as.numeric(newspores$flux)
newspores$flux.cat<-as.character(newspores$flux)
newspores$flux.cat<-as.factor(newspores$flux.cat)

ATLyo<-newspores[newspores$species == "ALT",]
SOLyo<-newspores[newspores$species == "SOL",]

ALTlim<-ATLyo[ATLyo$TotalSpores >=100,]
SOLlim<-SOLyo[SOLyo$TotalSpores >=40,]

View(head(ALTlim[ALTlim$hour == 0,], n=100))

hist(ALTlim[ALTlim$hour != 0,]$TotalSpores)
hist(SOLlim[SOLlim$hour != 0,]$TotalSpores)

hist(ALTlim[ALTlim$experiment %in% c("Exp3") & ALTlim$hour.rescale %in% c(2) & ALTlim$base %in% c("b4"),]$TotalSpores)


View(as.data.frame(list.files("~/Dropbox/Madison/Biotron/Genoa")))

##########################################################################################################
# sploring
library(gridExtra)
library(ggplot2)
grid.arrange(
ggplot(ALTlim[ALTlim$lambda %in% c("UVA") & ALTlim$hour.rescale %in% c(1,2,3,4),], aes(x=RH, y=AliveSpores/TotalSpores, group=RH))+geom_boxplot(width=.8, outlier.alpha = 0)+theme_bw()+labs(x="RH", y="Prop. Germinated", title="ALT:UVA"),
ggplot(SOLlim[SOLlim$lambda %in% c("UVA") & SOLlim$hour.rescale %in% c(1,2,3,4),], aes(x=RH, y=AliveSpores/TotalSpores, group=RH))+geom_boxplot(width=.8, outlier.alpha = 0)+theme_bw()+labs(x="RH", y="Prop. Germinated", title="SOL:UVA"),
ggplot(ALTlim[ALTlim$lambda %in% c("UVB") & ALTlim$hour.rescale %in% c(1,2,3,4),], aes(x=RH, y=AliveSpores/TotalSpores, group=RH))+geom_boxplot(width=.8, outlier.alpha = 0)+theme_bw()+labs(x="RH", y="Prop. Germinated", title="ALT:UVB"),
ggplot(SOLlim[SOLlim$lambda %in% c("UVB") & SOLlim$hour.rescale %in% c(1,2,3,4),], aes(x=RH, y=AliveSpores/TotalSpores, group=RH))+geom_boxplot(width=.8, outlier.alpha = 0)+theme_bw()+labs(x="RH", y="Prop. Germinated", title="SOL:UVB"),
ggplot(ALTlim[ALTlim$lambda %in% c("b21") & ALTlim$hour.rescale %in% c(1,2,3,4),], aes(x=RH, y=AliveSpores/TotalSpores, group=RH))+geom_boxplot(width=.8, outlier.alpha = 0)+theme_bw()+labs(x="RH", y="Prop. Germinated", title="ALT:dark"),
ggplot(SOLlim[SOLlim$lambda %in% c("b21") & SOLlim$hour.rescale %in% c(1,2,3,4),], aes(x=RH, y=AliveSpores/TotalSpores, group=RH))+geom_boxplot(width=.8, outlier.alpha = 0)+theme_bw()+labs(x="RH", y="Prop. Germinated", title="SOL:dark"),
nrow=3, ncol=2)

grid.arrange(
ggplot(ALTlim[ALTlim$lambda %in% c("UVA") & ALTlim$hour.rescale %in% c(1,2,3,4),], aes(x=temp, y=AliveSpores/TotalSpores, group=temp))+geom_boxplot(width=.8, outlier.alpha = 0)+theme_bw()+labs(x="T", y="Prop. Germinated", title="ALT:UVA"),
ggplot(SOLlim[SOLlim$lambda %in% c("UVA") & SOLlim$hour.rescale %in% c(1,2,3,4),], aes(x=temp, y=AliveSpores/TotalSpores, group=temp))+geom_boxplot(width=.8, outlier.alpha = 0)+theme_bw()+labs(x="T", y="Prop. Germinated", title="SOL:UVA"),
ggplot(ALTlim[ALTlim$lambda %in% c("UVB") & ALTlim$hour.rescale %in% c(1,2,3,4),], aes(x=temp, y=AliveSpores/TotalSpores, group=temp))+geom_boxplot(width=.8, outlier.alpha = 0)+theme_bw()+labs(x="T", y="Prop. Germinated", title="ALT:UVB"),
ggplot(SOLlim[SOLlim$lambda %in% c("UVB") & SOLlim$hour.rescale %in% c(1,2,3,4),], aes(x=temp, y=AliveSpores/TotalSpores, group=temp))+geom_boxplot(width=.8, outlier.alpha = 0)+theme_bw()+labs(x="T", y="Prop. Germinated", title="SOL:UVB"),
ggplot(ALTlim[ALTlim$lambda %in% c("b21") & ALTlim$hour.rescale %in% c(1,2,3,4),], aes(x=temp, y=AliveSpores/TotalSpores, group=temp))+geom_boxplot(width=.8, outlier.alpha = 0)+theme_bw()+labs(x="T", y="Prop. Germinated", title="ALT:dark"),
ggplot(SOLlim[SOLlim$lambda %in% c("b21") & SOLlim$hour.rescale %in% c(1,2,3,4),], aes(x=temp, y=AliveSpores/TotalSpores, group=temp))+geom_boxplot(width=.8, outlier.alpha = 0)+theme_bw()+labs(x="T", y="Prop. Germinated", title="SOL:dark"),
nrow=3, ncol=2)

ggplot(SOLlim[SOLlim$lambda %in% c("UVA") & SOLlim$hour.rescale %in% c(1,2,3,4),], aes(x=temp, y=AliveSpores/TotalSpores, fill=experiment))+
  geom_boxplot(width=.8, outlier.alpha = 0, position = position_dodge2(preserve = "single"))+theme_bw()+labs(x="T", y="Prop. Germinated", title="SOL:UVA")+
  scale_fill_discrete(limits=c("Exp1.1", "Exp2.1", "Exp3", "Exp4", "Exp5", "Exp6", "Exp7", 'Exp8', "Exp9", "Exp10"))+scale_colour_discrete(palette="Paired")


library(magrittr)
library(dplyr)
# summarize by BLOCK
ALTlimRH <- ALTlim %>%
  group_by(hour.rescale, lambda, BLOCK, RH) %>%
  summarise( 
    n=n(),
    mean.RH=mean(AliveSpores/TotalSpores),
    sd.RH=sd(AliveSpores/TotalSpores)
  ) %>%
  mutate (se.RH=sd.RH/sqrt(n)) 

SOLlimRH <- SOLlim %>%
  group_by(hour.rescale, lambda, BLOCK, RH) %>%
  summarise( 
    n=n(),
    mean.RH=mean(AliveSpores/TotalSpores),
    sd.RH=sd(AliveSpores/TotalSpores)
  ) %>%
  mutate (se.RH=sd.RH/sqrt(n)) 


ggplot(ALTlimRH[ALTlimRH$lambda %in% c("UVA") & ALTlimRH$hour.rescale %in% c(1,2,3,4),], aes(x=RH, y=mean.RH, group=RH))+geom_boxplot(width=.8, outlier.alpha = 0)+theme_bw()+labs(x="RH", y="Prop. Germinated", title="ALT:UVA")+ylim(0,1)
ggplot(ALTlimRH[ALTlimRH$lambda %in% c("UVB") & ALTlimRH$hour.rescale %in% c(1,2,3,4),], aes(x=RH, y=mean.RH, group=RH))+geom_boxplot(width=.8, outlier.alpha = 0)+theme_bw()+labs(x="RH", y="Prop. Germinated", title="ALT:UVB")+ylim(0,1)
ggplot(ALTlimRH[ALTlimRH$lambda %in% c("b21") & ALTlimRH$hour.rescale %in% c(1,2,3,4),], aes(x=RH, y=mean.RH, group=RH))+geom_boxplot(width=.8, outlier.alpha = 0)+theme_bw()+labs(x="RH", y="Prop. Germinated", title="ALT:dark")+ylim(0,1)

ggplot(SOLlimRH[SOLlimRH$lambda %in% c("UVA") & SOLlimRH$hour.rescale %in% c(1,2,3,4),], aes(x=RH, y=mean.RH, group=RH))+geom_boxplot(width=.8, outlier.alpha = 0)+theme_bw()+labs(x="RH", y="Prop. Germinated", title="SOL:UVA")+ylim(0,1)
ggplot(SOLlimRH[SOLlimRH$lambda %in% c("UVB") & SOLlimRH$hour.rescale %in% c(1,2,3,4),], aes(x=RH, y=mean.RH, group=RH))+geom_boxplot(width=.8, outlier.alpha = 0)+theme_bw()+labs(x="RH", y="Prop. Germinated", title="SOL:UVB")+ylim(0,1)
ggplot(SOLlimRH[SOLlimRH$lambda %in% c("b21") & SOLlimRH$hour.rescale %in% c(1,2,3,4),], aes(x=RH, y=mean.RH, group=RH))+geom_boxplot(width=.8, outlier.alpha = 0)+theme_bw()+labs(x="RH", y="Prop. Germinated", title="SOL:dark")+ylim(0,1)

grid.arrange(
ggplot(ALTlimRH[ALTlimRH$lambda %in% c("UVA") & ALTlimRH$hour.rescale %in% c(1,2,3,4),], aes(x=RH, y=mean.RH, group=RH))+geom_boxplot(width=.8, outlier.alpha = 0)+theme_bw()+labs(x="RH", y="Prop. Germinated", title="ALT:UVA")+ylim(0,1),
ggplot(SOLlimRH[SOLlimRH$lambda %in% c("UVA") & SOLlimRH$hour.rescale %in% c(1,2,3,4),], aes(x=RH, y=mean.RH, group=RH))+geom_boxplot(width=.8, outlier.alpha = 0)+theme_bw()+labs(x="RH", y="Prop. Germinated", title="SOL:UVA")+ylim(0,1),
ggplot(ALTlimRH[ALTlimRH$lambda %in% c("UVB") & ALTlimRH$hour.rescale %in% c(1,2,3,4),], aes(x=RH, y=mean.RH, group=RH))+geom_boxplot(width=.8, outlier.alpha = 0)+theme_bw()+labs(x="RH", y="Prop. Germinated", title="ALT:UVB")+ylim(0,1),
ggplot(SOLlimRH[SOLlimRH$lambda %in% c("UVB") & SOLlimRH$hour.rescale %in% c(1,2,3,4),], aes(x=RH, y=mean.RH, group=RH))+geom_boxplot(width=.8, outlier.alpha = 0)+theme_bw()+labs(x="RH", y="Prop. Germinated", title="SOL:UVB")+ylim(0,1),
ggplot(ALTlimRH[ALTlimRH$lambda %in% c("b21") & ALTlimRH$hour.rescale %in% c(1,2,3,4),], aes(x=RH, y=mean.RH, group=RH))+geom_boxplot(width=.8, outlier.alpha = 0)+theme_bw()+labs(x="RH", y="Prop. Germinated", title="ALT:dark")+ylim(0,1),
ggplot(SOLlimRH[SOLlimRH$lambda %in% c("b21") & SOLlimRH$hour.rescale %in% c(1,2,3,4),], aes(x=RH, y=mean.RH, group=RH))+geom_boxplot(width=.8, outlier.alpha = 0)+theme_bw()+labs(x="RH", y="Prop. Germinated", title="SOL:dark")+ylim(0,1),
nrow=3, ncol=2)



ALTlimT <- ALTlim %>%
  group_by(hour.rescale, lambda, BLOCK, temp) %>%
  summarise( 
    n=n(),
    mean.T=mean(AliveSpores/TotalSpores),
    sd.T=sd(AliveSpores/TotalSpores)
  ) %>%
  mutate (se.T=sd.T/sqrt(n)) 

SOLlimT <- SOLlim %>%
  group_by(hour.rescale, lambda, BLOCK, temp) %>%
  summarise( 
    n=n(),
    mean.T=mean(AliveSpores/TotalSpores),
    sd.T=sd(AliveSpores/TotalSpores)
  ) %>%
  mutate (se.T=sd.T/sqrt(n)) 

ggplot(ALTlimT[ALTlimT$lambda %in% c("UVA") & ALTlimT$hour.rescale %in% c(1,2,3,4),], aes(x=temp, y=mean.T, group=temp))+geom_boxplot(width=.8, outlier.alpha = 0)+theme_bw()+labs(x="T", y="Prop. Germinated", title="ALT:UVA")+ylim(0,1)
ggplot(ALTlimT[ALTlimT$lambda %in% c("UVB") & ALTlimT$hour.rescale %in% c(1,2,3,4),], aes(x=temp, y=mean.T, group=temp))+geom_boxplot(width=.8, outlier.alpha = 0)+theme_bw()+labs(x="T", y="Prop. Germinated", title="ALT:UVB")+ylim(0,1)
ggplot(ALTlimT[ALTlimT$lambda %in% c("b21") & ALTlimT$hour.rescale %in% c(1,2,3,4),], aes(x=temp, y=mean.T, group=temp))+geom_boxplot(width=.8, outlier.alpha = 0)+theme_bw()+labs(x="T", y="Prop. Germinated", title="ALT:dark")+ylim(0,1)

ggplot(SOLlimT[SOLlimT$lambda %in% c("UVA") & SOLlimT$hour.rescale %in% c(1,2,3,4),], aes(x=temp, y=mean.T, group=temp))+geom_boxplot(width=.8, outlier.alpha = 0)+theme_bw()+labs(x="T", y="Prop. Germinated", title="SOL:UVA")+ylim(0,1)
ggplot(SOLlimT[SOLlimT$lambda %in% c("UVB") & SOLlimT$hour.rescale %in% c(1,2,3,4),], aes(x=temp, y=mean.T, group=temp))+geom_boxplot(width=.8, outlier.alpha = 0)+theme_bw()+labs(x="T", y="Prop. Germinated", title="SOL:UVB")+ylim(0,1)
ggplot(SOLlimT[SOLlimT$lambda %in% c("b21") & SOLlimT$hour.rescale %in% c(1,2,3,4),], aes(x=temp, y=mean.T, group=temp))+geom_boxplot(width=.8, outlier.alpha = 0)+theme_bw()+labs(x="T", y="Prop. Germinated", title="SOL:dark")+ylim(0,1)

grid.arrange(
ggplot(ALTlimT[ALTlimT$lambda %in% c("UVA") & ALTlimT$hour.rescale %in% c(1,2,3,4),], aes(x=temp, y=mean.T, group=temp))+geom_boxplot(width=.8, outlier.alpha = 0)+theme_bw()+labs(x="T", y="Prop. Germinated", title="ALT:UVA")+ylim(0,1),
ggplot(SOLlimT[SOLlimT$lambda %in% c("UVA") & SOLlimT$hour.rescale %in% c(1,2,3,4),], aes(x=temp, y=mean.T, group=temp))+geom_boxplot(width=.8, outlier.alpha = 0)+theme_bw()+labs(x="T", y="Prop. Germinated", title="SOL:UVA")+ylim(0,1),
ggplot(ALTlimT[ALTlimT$lambda %in% c("UVB") & ALTlimT$hour.rescale %in% c(1,2,3,4),], aes(x=temp, y=mean.T, group=temp))+geom_boxplot(width=.8, outlier.alpha = 0)+theme_bw()+labs(x="T", y="Prop. Germinated", title="ALT:UVB")+ylim(0,1),
ggplot(SOLlimT[SOLlimT$lambda %in% c("UVB") & SOLlimT$hour.rescale %in% c(1,2,3,4),], aes(x=temp, y=mean.T, group=temp))+geom_boxplot(width=.8, outlier.alpha = 0)+theme_bw()+labs(x="T", y="Prop. Germinated", title="SOL:UVB")+ylim(0,1),
ggplot(ALTlimT[ALTlimT$lambda %in% c("b21") & ALTlimT$hour.rescale %in% c(1,2,3,4),], aes(x=temp, y=mean.T, group=temp))+geom_boxplot(width=.8, outlier.alpha = 0)+theme_bw()+labs(x="T", y="Prop. Germinated", title="ALT:dark")+ylim(0,1),
ggplot(SOLlimT[SOLlimT$lambda %in% c("b21") & SOLlimT$hour.rescale %in% c(1,2,3,4),], aes(x=temp, y=mean.T, group=temp))+geom_boxplot(width=.8, outlier.alpha = 0)+theme_bw()+labs(x="T", y="Prop. Germinated", title="SOL:dark")+ylim(0,1),
nrow=3, ncol=2)

# Boxplot by experiment BLOCK mean AND colored by T

ALTlimTex <- ALTlim %>%
  group_by(hour.rescale, lambda, BLOCK, temp, experiment, RH) %>%
  summarise( 
    n=n(),
    mean.T=mean(AliveSpores/TotalSpores),
    sd.T=sd(AliveSpores/TotalSpores)
  ) %>%
  mutate (se.T=sd.T/sqrt(n)) 

SOLlimTex <- SOLlim %>%
  group_by(hour.rescale, lambda, BLOCK, temp, experiment, RH) %>%
  summarise( 
    n=n(),
    mean.T=mean(AliveSpores/TotalSpores),
    sd.T=sd(AliveSpores/TotalSpores)
  ) %>%
  mutate (se.T=sd.T/sqrt(n)) 

#colorful box plots, not t=0
grid.arrange(
ggplot(ALTlimTex[ALTlimTex$lambda %in% c("UVA") & ALTlimTex$hour.rescale %in% c(1,2,3,4),], aes(x=temp, y=mean.T, fill=as.factor(RH)))+
  geom_boxplot(width=.8, outlier.alpha = 0, position = position_dodge2(preserve = "single"))+theme_bw()+labs(x="T", y="Prop. Germinated", title="ALT:UVA")+
  ylim(0,1)+theme(legend.position = "none")+scale_fill_manual(values=c("#68287e","#6666ff", "#a95fa6", "#555a9a")),
ggplot(SOLlimTex[SOLlimTex$lambda %in% c("UVA") & SOLlimTex$hour.rescale %in% c(1,2,3,4),], aes(x=temp, y=mean.T, fill=as.factor(RH)))+
  geom_boxplot(width=.8, outlier.alpha = 0, position = position_dodge2(preserve = "single"))+theme_bw()+labs(x="T", y="Prop. Germinated", title="SOL:UVA")+
  ylim(0,1)+theme(legend.position = "none")+scale_fill_manual(values=c("#dc4048","#f6821f", "#feb913", "#ed9fa3")),
ggplot(ALTlimTex[ALTlimTex$lambda %in% c("UVB") & ALTlimTex$hour.rescale %in% c(1,2,3,4),], aes(x=temp, y=mean.T, fill=as.factor(RH)))+
  geom_boxplot(width=.8, outlier.alpha = 0, position = position_dodge2(preserve = "single"))+theme_bw()+labs(x="T", y="Prop. Germinated", title="ALT:UVB")+
  ylim(0,1)+theme(legend.position = "none")+scale_fill_manual(values=c("#68287e","#6666ff", "#a95fa6", "#555a9a")),
ggplot(SOLlimTex[SOLlimTex$lambda %in% c("UVB") & SOLlimTex$hour.rescale %in% c(1,2,3,4),], aes(x=temp, y=mean.T, fill=as.factor(RH)))+
  geom_boxplot(width=.8, outlier.alpha = 0, position = position_dodge2(preserve = "single"))+theme_bw()+labs(x="T", y="Prop. Germinated", title="SOL:UVB")+
  ylim(0,1)+theme(legend.position = "none")+scale_fill_manual(values=c("#dc4048","#f6821f", "#feb913", "#ed9fa3")),
ggplot(ALTlimTex[ALTlimTex$lambda %in% c("b21") & ALTlimTex$hour.rescale %in% c(1,2,3,4),], aes(x=temp, y=mean.T, fill=as.factor(RH)))+
  geom_boxplot(width=.8, outlier.alpha = 0, position = position_dodge2(preserve = "single"))+theme_bw()+labs(x="T", y="Prop. Germinated", title="ALT:dark", fill="RH")+
  ylim(0,1)+theme(legend.position = "none")+scale_fill_manual(values=c("#68287e","#6666ff", "#a95fa6", "#555a9a")),
ggplot(SOLlimTex[SOLlimTex$lambda %in% c("b21") & SOLlimTex$hour.rescale %in% c(1,2,3,4),], aes(x=temp, y=mean.T, fill=as.factor(RH)))+
  geom_boxplot(width=.8, outlier.alpha = 0, position = position_dodge2(preserve = "single"))+theme_bw()+labs(x="T", y="Prop. Germinated", title="SOL:dark", fill="RH")+
  ylim(0,1)+theme(legend.position = "none")+scale_fill_manual(values=c("#dc4048","#f6821f", "#feb913", "#ed9fa3")),
nrow=3, ncol=2)

#look at just t=72
grid.arrange(
  ggplot(ALTlimTex[ALTlimTex$lambda %in% c("UVA") & ALTlimTex$hour.rescale %in% c(3),], aes(x=temp, y=mean.T, fill=as.factor(RH)))+
    geom_boxplot(width=.8, outlier.alpha = 0, position = position_dodge2(preserve = "single"))+theme_bw()+labs(x="T", y="Prop. Germinated", title="ALT:UVA")+
    ylim(0,1)+theme(legend.position = "none")+scale_fill_manual(values=c("#68287e","#6666ff", "#a95fa6", "#555a9a")),
  ggplot(SOLlimTex[SOLlimTex$lambda %in% c("UVA") & SOLlimTex$hour.rescale %in% c(3),], aes(x=temp, y=mean.T, fill=as.factor(RH)))+
    geom_boxplot(width=.8, outlier.alpha = 0, position = position_dodge2(preserve = "single"))+theme_bw()+labs(x="T", y="Prop. Germinated", title="SOL:UVA")+
    ylim(0,1)+theme(legend.position = "none")+scale_fill_manual(values=c("#dc4048","#f6821f", "#feb913", "#ed9fa3")),
  ggplot(ALTlimTex[ALTlimTex$lambda %in% c("UVB") & ALTlimTex$hour.rescale %in% c(3),], aes(x=temp, y=mean.T, fill=as.factor(RH)))+
    geom_boxplot(width=.8, outlier.alpha = 0, position = position_dodge2(preserve = "single"))+theme_bw()+labs(x="T", y="Prop. Germinated", title="ALT:UVB")+
    ylim(0,1)+theme(legend.position = "none")+scale_fill_manual(values=c("#68287e","#6666ff", "#a95fa6", "#555a9a")),
  ggplot(SOLlimTex[SOLlimTex$lambda %in% c("UVB") & SOLlimTex$hour.rescale %in% c(3),], aes(x=temp, y=mean.T, fill=as.factor(RH)))+
    geom_boxplot(width=.8, outlier.alpha = 0, position = position_dodge2(preserve = "single"))+theme_bw()+labs(x="T", y="Prop. Germinated", title="SOL:UVB")+
    ylim(0,1)+theme(legend.position = "none")+scale_fill_manual(values=c("#dc4048","#f6821f", "#feb913", "#ed9fa3")),
  ggplot(ALTlimTex[ALTlimTex$lambda %in% c("b21") & ALTlimTex$hour.rescale %in% c(3),], aes(x=temp, y=mean.T, fill=as.factor(RH)))+
    geom_boxplot(width=.8, outlier.alpha = 0, position = position_dodge2(preserve = "single"))+theme_bw()+labs(x="T", y="Prop. Germinated", title="ALT:dark", fill="RH")+
    ylim(0,1)+theme(legend.position = "none")+scale_fill_manual(values=c("#68287e","#6666ff", "#a95fa6", "#555a9a")),
  ggplot(SOLlimTex[SOLlimTex$lambda %in% c("b21") & SOLlimTex$hour.rescale %in% c(3),], aes(x=temp, y=mean.T, fill=as.factor(RH)))+
    geom_boxplot(width=.8, outlier.alpha = 0, position = position_dodge2(preserve = "single"))+theme_bw()+labs(x="T", y="Prop. Germinated", title="SOL:dark", fill="RH")+
    ylim(0,1)+theme(legend.position = "none")+scale_fill_manual(values=c("#dc4048","#f6821f", "#feb913", "#ed9fa3")),
  nrow=3, ncol=2)

#just look at 90% at UVA
grid.arrange(
  ggplot(ALTlimTex[ALTlimTex$lambda %in% c("UVA") & ALTlimTex$hour.rescale %in% c(1,2,3,4) & ALTlimTex$RH %in% c(90),], aes(x=temp, y=mean.T, fill=as.factor(RH)))+
    geom_boxplot(width=.8, outlier.alpha = 0, position = position_dodge2(preserve = "single"))+theme_bw()+labs(x="T", y="Prop. Germinated", title="ALT:UVA")+
    ylim(0,1)+theme(legend.position = "none")+scale_fill_manual(values=c("#68287e")),
  ggplot(SOLlimTex[SOLlimTex$lambda %in% c("UVA") & SOLlimTex$hour.rescale %in% c(1,2,3,4) & SOLlimTex$RH %in% c(90),], aes(x=temp, y=mean.T, fill=as.factor(RH)))+
    geom_boxplot(width=.8, outlier.alpha = 0, position = position_dodge2(preserve = "single"))+theme_bw()+labs(x="T", y="Prop. Germinated", title="SOL:UVA")+
    ylim(0,1)+theme(legend.position = "none")+scale_fill_manual(values=c("#dc4048")),
  nrow=1, ncol=2)

#just look at 90% at UVA at t=96
grid.arrange(
  ggplot(ALTlimTex[ALTlimTex$lambda %in% c("UVA") & ALTlimTex$hour.rescale %in% c(4) & ALTlimTex$RH %in% c(90),], aes(x=temp, y=mean.T, fill=as.factor(RH)))+
    geom_boxplot(width=.8, outlier.alpha = 0, position = position_dodge2(preserve = "single"))+theme_bw()+labs(x="T", y="Prop. Germinated", title="ALT:UVA - 90%RH - Day4")+
    ylim(0,1)+theme(legend.position = "none")+scale_fill_manual(values=c("#6666ff")),
  ggplot(SOLlimTex[SOLlimTex$lambda %in% c("UVA") & SOLlimTex$hour.rescale %in% c(4) & SOLlimTex$RH %in% c(90),], aes(x=temp, y=mean.T, fill=as.factor(RH)))+
    geom_boxplot(width=.8, outlier.alpha = 0, position = position_dodge2(preserve = "single"))+theme_bw()+labs(x="T", y="Prop. Germinated", title="SOL:UVA - 90%RH - Day4")+
    ylim(0,1)+theme(legend.position = "none")+scale_fill_manual(values=c("#ff4c4c")),
  nrow=1, ncol=2)

mean(ALTlimTex[ALTlimTex$lambda %in% c("UVA") & ALTlimTex$hour.rescale %in% c(3) & ALTlimTex$RH %in% c(90) & ALTlimTex$temp %in% c(15),]$mean.T)

#just look at 90% at UVA and show times
grid.arrange(
ggplot(ALTlimTex[ALTlimTex$lambda %in% c("UVA") & ALTlimTex$hour.rescale %in% c(1,2,3,4) & ALTlimTex$RH %in% c(90) ,], aes(x=temp, y=mean.T, fill=as.factor(hour.rescale)))+
  geom_boxplot(width=.8, outlier.alpha = 0, position = position_dodge2(preserve = "single"))+theme_bw()+labs(x="T", y="Prop. Germinated", title="ALT:UVA")+
  ylim(0,1)+theme(legend.position = "none"),
ggplot(SOLlimTex[SOLlimTex$lambda %in% c("UVA") & SOLlimTex$hour.rescale %in% c(1,2,3,4) & SOLlimTex$RH %in% c(90) ,], aes(x=temp, y=mean.T, fill=as.factor(hour.rescale)))+
  geom_boxplot(width=.8, outlier.alpha = 0, position = position_dodge2(preserve = "single"))+theme_bw()+labs(x="T", y="Prop. Germinated", title="ALT:UVA")+
  ylim(0,1)+theme(legend.position = "none"),
nrow=1, ncol=2)

############################################ CLEAN BOXPLOTS by TIME

ALTlimTexbase <- ALTlim %>%
  group_by(hour.rescale, lambda, BLOCK, temp, experiment, RH, base) %>%
  summarise( 
    n=n(),
    mean.T=mean(AliveSpores/TotalSpores),
    sd.T=sd(AliveSpores/TotalSpores)
  ) %>%
  mutate (se.T=sd.T/sqrt(n)) 

SOLlimTexbase <- SOLlim %>%
  group_by(hour.rescale, lambda, BLOCK, temp, experiment, RH, base) %>%
  summarise( 
    n=n(),
    mean.T=mean(AliveSpores/TotalSpores),
    sd.T=sd(AliveSpores/TotalSpores)
  ) %>%
  mutate (se.T=sd.T/sqrt(n)) 


grid.arrange(
  ggplot(ALTlimTexbase[ALTlimTexbase$lambda %in% c("UVA") & ALTlimTexbase$hour.rescale %in% c(3) & ALTlimTexbase$RH %in% c(90),], aes(x=temp, y=mean.T, fill=as.factor(RH)))+
    geom_boxplot(width=.8, outlier.alpha = 0, position = position_dodge2(preserve = "single"))+theme_bw()+labs(x="T", y="Prop. Germinated", title="ALT:UVA - 90%RH - Day3")+
    ylim(0,1)+theme(legend.position = "none")+scale_fill_manual(values=c("#6666ff")),
  ggplot(SOLlimTexbase[SOLlimTexbase$lambda %in% c("UVA") & SOLlimTexbase$hour.rescale %in% c(3) & SOLlimTexbase$RH %in% c(90),], aes(x=temp, y=mean.T, fill=as.factor(RH)))+
    geom_boxplot(width=.8, outlier.alpha = 0, position = position_dodge2(preserve = "single"))+theme_bw()+labs(x="T", y="Prop. Germinated", title="SOL:UVA - 90%RH - Day3")+
    ylim(0,1)+theme(legend.position = "none")+scale_fill_manual(values=c("#ff4c4c")),
  ggplot(ALTlimTexbase[ALTlimTexbase$lambda %in% c("UVA") & ALTlimTexbase$hour.rescale %in% c(4) & ALTlimTexbase$RH %in% c(90),], aes(x=temp, y=mean.T, fill=as.factor(RH)))+
    geom_boxplot(width=.8, outlier.alpha = 0, position = position_dodge2(preserve = "single"))+theme_bw()+labs(x="T", y="Prop. Germinated", title="ALT:UVA - 90%RH - Day4")+
    ylim(0,1)+theme(legend.position = "none")+scale_fill_manual(values=c("#6666ff")),
  ggplot(SOLlimTexbase[SOLlimTexbase$lambda %in% c("UVA") & SOLlimTexbase$hour.rescale %in% c(4) & SOLlimTexbase$RH %in% c(90) & SOLlimTexbase$base %in% c("b1", "b2", "b3", "b4", "b5", "b6"),], aes(x=temp, y=mean.T, fill=as.factor(RH)))+
    geom_boxplot(width=.8, outlier.alpha = 0, position = position_dodge2(preserve = "single"))+theme_bw()+labs(x="T", y="Prop. Germinated", title="SOL:UVA - 90%RH - Day4")+
    ylim(0,1)+theme(legend.position = "none")+scale_fill_manual(values=c("#ff4c4c")),
  nrow=2, ncol=2)

grid.arrange(
ggplot(ALTlimTexbase[ALTlimTexbase$lambda %in% c("UVA") & ALTlimTexbase$hour.rescale %in% c(3, 4) & ALTlimTexbase$RH %in% c(90),], aes(x=temp, y=mean.T, fill=as.factor(hour.rescale)))+
  geom_boxplot(width=.8, outlier.alpha = 0, position = position_dodge2(preserve = "single"))+theme_bw()+labs(x="T", y="Prop. Germinated", title="ALT:UVA - 90%RH - Day3-4")+
  ylim(0,1)+theme(legend.position = "none")+scale_fill_manual(values=c("#6666ff", "#9999ff")),
ggplot(SOLlimTexbase[SOLlimTexbase$lambda %in% c("UVA") & SOLlimTexbase$hour.rescale %in% c(3, 4) & SOLlimTexbase$RH %in% c(90) & SOLlimTexbase$base %in% c("b1", "b2", "b3", "b4", "b5", "b6"),], aes(x=temp, y=mean.T, fill=as.factor(hour.rescale)))+
  geom_boxplot(width=.8, outlier.alpha = 0, position = position_dodge2(preserve = "single"))+theme_bw()+labs(x="T", y="Prop. Germinated", title="ALT:UVA - 90%RH - Day3-4")+
  ylim(0,1)+theme(legend.position = "none")+scale_fill_manual(values=c("#ff4c4c", "#ff9999")),
nrow=1, ncol=2)


# library("plot3D")
# x<-as.double(as.character(ALTlimTex[ALTlimTex$lambda %in% c("UVA") & ALTlimTex$hour.rescale %in% c(1,2,3,4),]$RH))
# y<-as.double(as.character(ALTlimTex[ALTlimTex$lambda %in% c("UVA") & ALTlimTex$hour.rescale %in% c(1,2,3,4),]$temp))
# z<-ALTlimTex[ALTlimTex$lambda %in% c("UVA") & ALTlimTex$hour.rescale %in% c(1,2,3,4),]$mean.T
#   
#   
# scatter3D(x,y,z,
#           phi = 0, bty ="g",
#           xlab = "RH",
#           ylab ="T", 
#           zlab = "Prop. Alive",
#           ticktype = "detailed")
# 
# ALTUVA <- ALTlimTex[ALTlimTex$lambda %in% c("UVA") & ALTlimTex$hour.rescale %in% c(1,2,3,4),] %>%
#   group_by(RH, temp) %>%
#   summarise( 
#     n=n(),
#     mean.pr=mean(mean.T),
#     sd.pr=sd(mean.T)
#   ) %>%
#   mutate (se.pr=sd.pr/sqrt(n)) 
# 
# #HOW TO MAKE A LIST INTO A SYM MATRIX FINALLLLLLY
# raw<-as.data.frame(cbind(cha$comparisons, cha$P.adjusted))
# library(stringr)
# foo <- data.frame(do.call('rbind', strsplit(as.character(raw[,1]),'-',fixed=TRUE)))
# foo$padj<-raw[,2]
# 
# foo<-as.data.frame(ALTUVA[,c(1,2,4)])
# colnames(foo)<-c("RH", "temp", "prop")
# 
# fuck <- reshape(foo, direction="wide", idvar="temp", timevar="RH")
# rownames(fuck)<-fuck[,1]
# fuck<-fuck[,-1]
# fuckmat<-as.matrix(fuck)
# 
# hist3D (x = 1:4, y = 1:4, z = fuckmat,
#         bty = "g", phi = 20,  theta = -60,
#         xlab = "", ylab = "", zlab = "", main = "VADeaths",
#         col = "#0072B2", border = "black", shade = 0.8,
#         ticktype = "detailed", space = 0.15, d = 2, cex.axis = 1e-9)
# text3D(x = 1:4, y = rep(0.5, 4), z = rep(3, 4),
#        labels = rownames(fuckmat),
#        add = TRUE, adj = 0)


# REPRODUCE VIOLIN LINE PLOT EXTRAVANGANZA ******
ALTlimmury <- ALTlim %>%
  group_by(hour.rescale, lambda) %>%
  summarise( 
    n=n(),
    mean.pr=mean(AliveSpores/TotalSpores),
    sd.pr=sd(AliveSpores/TotalSpores)
  ) %>%
  mutate (se.pr=sd.pr/sqrt(n)) 

ALTlimmury.norm <- ALTlimmury %>%
  group_by(hour.rescale) %>%
  mutate(norm.mean.pr=mean.pr/mean(ALTlimmury$mean.pr[ALTlimmury$hour.rescale==0])) %>%
  mutate(norm.sd.pr=sd.pr/mean(ALTlimmury$mean.pr[ALTlimmury$hour.rescale==0])) %>%
  mutate(norm.se.pr=se.pr/mean(ALTlimmury$mean.pr[ALTlimmury$hour.rescale==0]))

ALTlimad<-ALTlim
ALTlimad$Perc<-ALTlimad$AliveSpores/ALTlimad$TotalSpores
ggplot()+
  geom_violin(data=ALTlimad, aes(x=hour.rescale, y=Perc, group=hour.rescale, fill=as.factor(hour.rescale)), alpha=.2, adjust=2.5, trim=T, bw=.01)+
  geom_errorbar(data=ALTlimmury.norm, aes(x=hour.rescale, ymin=mean.pr-(se.pr), ymax=mean.pr+(se.pr)), width=.2)+
  geom_line(data=ALTlimmury.norm, aes(x=hour.rescale, y=mean.pr, linetype=lambda),cex=1)+
  scale_fill_manual(values=c("red", "orange", "yellow", "darkgreen", "blue"))+theme_bw()+theme(legend.position = "none")

#violin by base mean
ALTlimbase <- ALTlim %>%
  group_by(hour.rescale, lambda, base) %>%
  summarise( 
    n=n(),
    mean.pr=mean(AliveSpores/TotalSpores),
    sd.pr=sd(AliveSpores/TotalSpores)
  ) %>%
  mutate (se.pr=sd.pr/sqrt(n)) 

ggplot()+
  geom_violin(data=ALTlimbase, aes(x=hour.rescale, y=mean.pr, group=hour.rescale, fill=as.factor(hour.rescale)), alpha=.2, adjust=2, trim=F, bw=.01)+
  geom_errorbar(data=ALTlimmury.norm, aes(x=hour.rescale, ymin=mean.pr-(se.pr), ymax=mean.pr+(se.pr)), width=.2)+
  geom_line(data=ALTlimmury.norm, aes(x=hour.rescale, y=mean.pr, linetype=lambda),cex=1)+
  scale_fill_manual(values=c("red", "orange", "yellow", "darkgreen", "blue"))+theme_bw()+theme(legend.position = "none")+
  scale_linetype_manual(values=c("solid", "dotdash", "dotted"))+labs(x="Day", y="Prop. Germinated", title="ALT")


#BAUM didnt like density plot, try box plot

ggplot(ALTlimbase, aes(x=factor(hour.rescale), y=mean.pr, fill=lambda))+geom_boxplot()

ALTlimbase.sanb21<-ALTlimbase[ALTlimbase$lambda != "b21",]
ALTlimmury.normsanb21<-ALTlimmury.norm[ALTlimmury.norm$lambda != "b21",]

# ggplot()+geom_boxplot(data=ALTlimbase.sanb21, aes(x=factor(hour.rescale), y=mean.pr, fill=lambda))+
#   geom_errorbar(data=ALTlimmury.normsanb21, aes(x=hour.rescale, ymin=mean.pr-(se.pr), ymax=mean.pr+(se.pr)), width=.2)+
#   geom_line(data=ALTlimmury.normsanb21, aes(x=hour.rescale, y=mean.pr, linetype=lambda),cex=1)+
#   theme_bw()+theme(legend.position = "none")+
#   scale_linetype_manual(values=c("solid", "dotdash", "dotted"))+labs(x="Day", y="Prop. Germinated", title="ALT")


ggplot(ALTlimbase.sanb21, aes(x=factor(hour.rescale), y=mean.pr, fill=lambda)) + 
  geom_boxplot() + 
  stat_summary(fun.y=mean, geom="line", aes(group=lambda))  + 
  stat_summary(fun.y=mean, geom="point")+
  theme_bw()+theme(legend.position = "none")+
  labs(x="Day", y="Prop. Germinated", title="ALT")





# ALTlimmury.exp <- ALTlim %>%
#   group_by(hour.rescale, lambda, experiment) %>%
#   summarise( 
#     n=n(),
#     mean.pr=mean(AliveSpores/TotalSpores),
#     sd.pr=sd(AliveSpores/TotalSpores)
#   ) %>%
#   mutate (se.pr=sd.pr/sqrt(n)) 
# 
# ALTlimmury.norm.exp <- ALTlimmury.exp %>%
#   group_by(hour.rescale) %>%
#   mutate(norm.mean.pr=mean.pr/mean(ALTlimmury.exp$mean.pr[ALTlimmury.exp$hour.rescale==0])) %>%
#   mutate(norm.sd.pr=sd.pr/mean(ALTlimmury.exp$mean.pr[ALTlimmury.exp$hour.rescale==0])) %>%
#   mutate(norm.se.pr=se.pr/mean(ALTlimmury.exp$mean.pr[ALTlimmury.exp$hour.rescale==0]))


# ALTlimad<-ALTlim
# ALTlimad$Perc<-ALTlimad$AliveSpores/ALTlimad$TotalSpores
# ggplot()+
#   geom_violin(data=ALTlimad, aes(x=hour.rescale, y=Perc, group=hour.rescale, fill=as.factor(hour.rescale)), alpha=.2, adjust=2.5, trim=F, bw=.007)+
#   geom_errorbar(data=ALTlimmury.norm, aes(x=hour.rescale, ymin=mean.pr-(se.pr), ymax=mean.pr+(se.pr)), width=.2)+
#   geom_line(data=ALTlimmury.norm, aes(x=hour.rescale, y=mean.pr, linetype=lambda),cex=1)+
#   scale_fill_manual(values=c("red", "orange", "yellow", "darkgreen", "blue"))+theme_bw()+theme(legend.position = "none")+
#   geom_point(data=ALTlimmury.norm.exp, aes(x=hour.rescale, y=norm.mean.pr, shape=lambda))
  
SOLlimmury <- SOLlim %>%
  group_by(hour.rescale, lambda) %>%
  summarise( 
    n=n(),
    mean.pr=mean(AliveSpores/TotalSpores),
    sd.pr=sd(AliveSpores/TotalSpores)
  ) %>%
  mutate (se.pr=sd.pr/sqrt(n)) 

SOLlimmury.norm <- SOLlimmury %>%
  group_by(hour.rescale) %>%
  mutate(norm.mean.pr=mean.pr/mean(SOLlimmury$mean.pr[SOLlimmury$hour.rescale==0])) %>%
  mutate(norm.sd.pr=sd.pr/mean(SOLlimmury$mean.pr[SOLlimmury$hour.rescale==0])) %>%
  mutate(norm.se.pr=se.pr/mean(SOLlimmury$mean.pr[SOLlimmury$hour.rescale==0]))

SOLlimad<-SOLlim
SOLlimad$Perc<-SOLlimad$AliveSpores/SOLlimad$TotalSpores
ggplot()+
  geom_violin(data=SOLlimad, aes(x=hour.rescale, y=Perc, group=hour.rescale, fill=as.factor(hour.rescale)), alpha=.2, adjust=2.5, trim=T, bw=.01)+
  geom_errorbar(data=SOLlimmury.norm, aes(x=hour.rescale, ymin=mean.pr-(se.pr), ymax=mean.pr+(se.pr)), width=.2)+
  geom_line(data=SOLlimmury.norm, aes(x=hour.rescale, y=mean.pr, linetype=lambda),cex=1)+
  scale_fill_manual(values=c("red", "orange", "yellow", "darkgreen", "blue"))+theme_bw()+theme(legend.position = "none")

#violin by base mean
SOLlimbase <- SOLlim %>%
  group_by(hour.rescale, lambda, base) %>%
  summarise( 
    n=n(),
    mean.pr=mean(AliveSpores/TotalSpores),
    sd.pr=sd(AliveSpores/TotalSpores)
  ) %>%
  mutate (se.pr=sd.pr/sqrt(n)) 

ggplot()+
  geom_violin(data=SOLlimbase, aes(x=hour.rescale, y=mean.pr, group=hour.rescale, fill=as.factor(hour.rescale)), alpha=.2, adjust=2, trim=F, bw=.01)+
  geom_errorbar(data=SOLlimmury.norm, aes(x=hour.rescale, ymin=mean.pr-(se.pr), ymax=mean.pr+(se.pr)), width=.2)+
  geom_line(data=SOLlimmury.norm, aes(x=hour.rescale, y=mean.pr, linetype=lambda),cex=1)+
  scale_fill_manual(values=c("red", "orange", "yellow", "darkgreen", "blue"))+theme_bw()+theme(legend.position = "none")+
  scale_linetype_manual(values=c("solid", "dotdash", "dotted"))+labs(x="Day", y="Prop. Germinated", title="SOL")


# together
grid.arrange(
  ggplot()+
    geom_violin(data=ALTlimad, aes(x=hour.rescale, y=Perc, group=hour.rescale, fill=as.factor(hour.rescale)), alpha=.2, adjust=2.5, trim=T, bw=.01)+
    geom_errorbar(data=ALTlimmury.norm, aes(x=hour.rescale, ymin=mean.pr-(se.pr), ymax=mean.pr+(se.pr)), width=.2)+
    geom_line(data=ALTlimmury.norm, aes(x=hour.rescale, y=mean.pr, linetype=lambda),cex=1)+
    scale_fill_manual(values=c("red", "orange", "yellow", "darkgreen", "blue"))+theme_bw()+theme(legend.position = "none")+
    scale_linetype_manual(values=c("solid", "dotdash", "dotted"))+labs(x="Day", y="Prop. Germinated", title="ALT"),
  
  ggplot()+
    geom_violin(data=SOLlimad, aes(x=hour.rescale, y=Perc, group=hour.rescale, fill=as.factor(hour.rescale)), alpha=.2, adjust=2.5, trim=T, bw=.01)+
    geom_errorbar(data=SOLlimmury.norm, aes(x=hour.rescale, ymin=mean.pr-(se.pr), ymax=mean.pr+(se.pr)), width=.2)+
    geom_line(data=SOLlimmury.norm, aes(x=hour.rescale, y=mean.pr, linetype=lambda),cex=1)+
    scale_fill_manual(values=c("red", "orange", "yellow", "darkgreen", "blue"))+theme_bw()+theme(legend.position = "none")+
    scale_linetype_manual(values=c("solid", "dotdash", "dotted"))+labs(x="Day", y="Prop. Germinated", title="SOL"),
  
  nrow=1, ncol=2
)

#THIS IS THE VERSION I LIKE
#violin by block
grid.arrange(
  ggplot()+
    geom_violin(data=ALTlimbase, aes(x=hour.rescale, y=mean.pr, group=hour.rescale, fill=as.factor(hour.rescale)), alpha=1, adjust=2, trim=F, bw=.01, color=NA)+
    geom_errorbar(data=ALTlimmury.norm, aes(x=hour.rescale, ymin=mean.pr-(se.pr), ymax=mean.pr+(se.pr)), width=.2)+
    geom_line(data=ALTlimmury.norm, aes(x=hour.rescale, y=mean.pr, linetype=lambda),cex=1)+
    scale_fill_manual(values=c("red", "orange", "yellow", "darkgreen", "blue"))+theme_bw()+theme(legend.position = "none")+
    scale_linetype_manual(values=c("solid", "dotdash", "dotted"))+labs(x="Day", y="Prop. Germinated", title="ALT"),
  
  ggplot()+
    geom_violin(data=SOLlimbase, aes(x=hour.rescale, y=mean.pr, group=hour.rescale, fill=as.factor(hour.rescale)), alpha=1, adjust=2, trim=F, bw=.01, color=NA)+
    geom_errorbar(data=SOLlimmury.norm, aes(x=hour.rescale, ymin=mean.pr-(se.pr), ymax=mean.pr+(se.pr)), width=.2)+
    geom_line(data=SOLlimmury.norm, aes(x=hour.rescale, y=mean.pr, linetype=lambda),cex=1)+
    scale_fill_manual(values=c("red", "orange", "yellow", "darkgreen", "blue"))+theme_bw()+theme(legend.position = "none")+
    scale_linetype_manual(values=c("solid", "dotdash", "dotted"))+labs(x="Day", y="Prop. Germinated", title="SOL"),
  
  
  nrow=1, ncol=2
)









##########################################################################################################
#
#                                   MODELING
#
##########################################################################################################
#produce summaries

library(dplyr)

ALTsum <- ALTlim %>%
  group_by(experiment, RH, temp, lambda, flux, BLOCK, hour.rescale, cart, base) %>%
  summarize(TotalSpores = mean(TotalSpores),
            AliveSpores = mean(AliveSpores))


SOLsum <- SOLlim %>%
  group_by(experiment, RH, temp, lambda, flux, BLOCK, hour.rescale, cart, base) %>%
  summarize(TotalSpores = mean(TotalSpores),
            AliveSpores = mean(AliveSpores))

##########################################################################################################
#######################           PDF SELECTION                               ############################
##########################################################################################################
#ALT
alt.mod.pois<-glmmTMB(AliveSpores~RH+hour.rescale+lambda+(1|experiment/cart:base)+offset(log(TotalSpores)),
              data=ALTsum, family=poisson)

alt.mod.nbin1<-glmmTMB(AliveSpores~RH+hour.rescale+lambda+(1|experiment/cart:base)+offset(log(TotalSpores)),
                      data=ALTsum, family=nbinom1)

alt.mod.nbin2<-glmmTMB(AliveSpores~RH+hour.rescale+lambda+(1|experiment/cart:base)+offset(log(TotalSpores)),
                      data=ALTsum, family=nbinom2)

AICtab(alt.mod.pois, alt.mod.nbin1, alt.mod.nbin2)
#nbin1 wins

# _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
#SOL
sol.mod.pois<-glmmTMB(AliveSpores~RH+hour.rescale+lambda+(1|experiment/cart:base)+offset(log(TotalSpores)),
                      data=SOLsum, family=poisson)

sol.mod.nbin1<-glmmTMB(AliveSpores~RH+hour.rescale+lambda+(1|experiment/cart:base)+offset(log(TotalSpores)),
                       data=SOLsum, family=nbinom1)

sol.mod.nbin2<-glmmTMB(AliveSpores~RH+hour.rescale+lambda+(1|experiment/cart:base)+offset(log(TotalSpores)),
                       data=SOLsum, family=nbinom2)

AICtab(sol.mod.pois, sol.mod.nbin1, sol.mod.nbin2)
#nbin2 wins


##########################################################################################################
#######################           PARAM SUBSETTING                            ############################
##########################################################################################################


alt.togon01<-glmmTMB(AliveSpores~RH+temp+hour.rescale+lambda+(1|experiment/cart:base)+offset(log(TotalSpores)),
                   data=ALTsum[ALTsum$hour.rescale !=0,], family=nbinom1)

alt.togon02<-glmmTMB(AliveSpores~RH+temp+hour.rescale+lambda+as.factor(flux)+(1|experiment/cart:base)+offset(log(TotalSpores)),
                   data=ALTsum[ALTsum$hour.rescale !=0,], family=nbinom1)

alt.togon03<-glmmTMB(AliveSpores~RH+temp+hour.rescale+(1|experiment/cart:base)+offset(log(TotalSpores)),
                   data=ALTsum[ALTsum$hour.rescale !=0,], family=nbinom1)



#ALT - RH

alt.mod1<-glmmTMB(AliveSpores~RH+hour.rescale+(1|experiment/cart:base)+offset(log(TotalSpores)),
                       data=ALTsum, family=nbinom1)

alt.mod2<-glmmTMB(AliveSpores~RH+hour.rescale+lambda+(1|experiment/cart:base)+offset(log(TotalSpores)),
                  data=ALTsum, family=nbinom1)


test<-glmmTMB(AliveSpores~RH+temp+hour.rescale+lambda+(1|experiment/cart:base)+offset(log(TotalSpores)),
                  data=ALTsum, family=nbinom1)
plot_model(test,show.values = TRUE, value.offset = .2, dot.size = 3, transform="exp", vline.color = "black")

test.no0<-glmmTMB(AliveSpores~RH+temp+hour.rescale+lambda+(1|experiment/cart:base)+offset(log(TotalSpores)),
              data=ALTsum[ALTsum$hour.rescale !=0,], family=nbinom1)
plot_model(test.no0,show.values = TRUE, value.offset = .2, dot.size = 3, transform="exp", vline.color = "black")

# alt.mod3<-glmmTMB(AliveSpores~RH+hour.rescale+lambda+as.factor(flux)+(1|experiment/cart:base)+offset(log(TotalSpores)),
#                   data=ALTsum, family=nbinom1)

# alt.mod4<-glmmTMB(AliveSpores~RH+hour.rescale+lambda+base+(1|experiment/cart:base)+offset(log(TotalSpores)),
#                   data=ALTsum, family=nbinom1)

# alt.mod5<-glmmTMB(AliveSpores~RH+hour.rescale+lambda+base+as.factor(flux)+(1|experiment/cart:base)+offset(log(TotalSpores)),
#                    data=ALTsum, family=nbinom1)

alt.mod6<-glmmTMB(AliveSpores~RH:lambda+hour.rescale+lambda+(1|experiment/cart:base)+offset(log(TotalSpores)),
                  data=ALTsum, family=nbinom1)

alt.mod2no0<-glmmTMB(AliveSpores~RH+hour.rescale+lambda+(1|experiment/cart:base)+offset(log(TotalSpores)),
                  data=ALTsum[ALTsum$hour.rescale != 0,], family=nbinom1)

alt.mod2.96<-glmmTMB(AliveSpores~RH+lambda+(1|experiment/cart:base)+offset(log(TotalSpores)),
                     data=ALTsum[ALTsum$hour.rescale == 4,], family=nbinom1)


AICtab(alt.mod1, alt.mod2,alt.mod6)
plot_model(alt.mod2,show.values = TRUE, value.offset = .2, dot.size = 3, transform="exp", vline.color = "black")
plot_model(alt.mod6,order.terms = c(1,2,3,5,6,4,8,9,7,11,12,10), show.values = TRUE, value.offset = .2, dot.size = 3, transform="exp", vline.color = "black")

plot_model(alt.mod2no0,show.values = TRUE, value.offset = .2, dot.size = 3, transform="exp", vline.color = "black")
plot_model(alt.mod2.96,show.values = TRUE, value.offset = .2, dot.size = 3, transform="exp", vline.color = "black")

# _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _

#ALT - T

alt.Tmod1<-glmmTMB(AliveSpores~temp+hour.rescale+(1|experiment/cart:base)+offset(log(TotalSpores)),
                  data=ALTsum, family=nbinom1)

alt.Tmod2<-glmmTMB(AliveSpores~temp+hour.rescale+lambda+(1|experiment/cart:base)+offset(log(TotalSpores)),
                  data=ALTsum, family=nbinom1)

# alt.Tmod3<-glmmTMB(AliveSpores~temp+hour.rescale+lambda+as.factor(flux)+(1|experiment/cart:base)+offset(log(TotalSpores)),
#                   data=ALTsum, family=nbinom1)

# alt.Tmod4<-glmmTMB(AliveSpores~temp+hour.rescale+lambda+base+(1|experiment/cart:base)+offset(log(TotalSpores)),
#                    data=ALTsum, family=nbinom1)

# alt.Tmod5<-glmmTMB(AliveSpores~temp:lambda+hour.rescale+lambda+(1|experiment/cart:base)+offset(log(TotalSpores)),
#                    data=ALTsum, family=nbinom1)

alt.Tmod2.no0<-glmmTMB(AliveSpores~temp+hour.rescale+lambda+(1|experiment/cart:base)+offset(log(TotalSpores)),
                   data=ALTsum[ALTsum$hour.rescale != 0,], family=nbinom1)

alt.Tmod2.96<-glmmTMB(AliveSpores~temp+lambda+(1|experiment/cart:base)+offset(log(TotalSpores)),
                       data=ALTsum[ALTsum$hour.rescale == 4,], family=nbinom1)


AICtab(alt.Tmod1, alt.Tmod2)
plot_model(alt.Tmod2,show.values = TRUE, value.offset = .2, dot.size = 3, transform="exp", vline.color = "black")

plot_model(alt.Tmod2.no0,show.values = TRUE, value.offset = .2, dot.size = 3, transform="exp", vline.color = "black")
plot_model(alt.Tmod2.96,show.values = TRUE, value.offset = .2, dot.size = 3, transform="exp", vline.color = "black")


# plot_model(alt.Tmod5,show.values = TRUE, value.offset = .2, dot.size = 3, transform="exp", vline.color = "black")

# _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
#ALT - T sub 90
ALTsum90<-ALTsum[ALTsum$RH == 90,]

alt.T90mod1<-glmmTMB(AliveSpores~temp+hour.rescale+(1|experiment/cart:base)+offset(log(TotalSpores)),
                   data=ALTsum90, family=nbinom1)

alt.T90mod2<-glmmTMB(AliveSpores~temp+hour.rescale+lambda+(1|experiment/cart:base)+offset(log(TotalSpores)),
                   data=ALTsum90, family=nbinom1)

# alt.T90mod3<-glmmTMB(AliveSpores~temp+hour.rescale+lambda+as.factor(flux)+(1|experiment/cart:base)+offset(log(TotalSpores)),
#                   data=ALTsum90, family=nbinom1)
# 
# alt.T90mod4<-glmmTMB(AliveSpores~temp+hour.rescale+lambda+base+(1|experiment/cart:base)+offset(log(TotalSpores)),
#                    data=ALTsum90, family=nbinom1)

# alt.T90mod5<-glmmTMB(AliveSpores~temp:lambda+hour.rescale+lambda+(1|experiment/cart:base)+offset(log(TotalSpores)),
#                      data=ALTsum90, family=nbinom1)

alt.T90mod2.no0<-glmmTMB(AliveSpores~temp+hour.rescale+lambda+(1|experiment/cart:base)+offset(log(TotalSpores)),
                     data=ALTsum90[ALTsum90$hour.rescale !=0,], family=nbinom1)

alt.T90mod2.96<-glmmTMB(AliveSpores~temp+lambda+(1|experiment/cart:base)+offset(log(TotalSpores)),
                         data=ALTsum90[ALTsum90$hour.rescale ==4,], family=nbinom1)

AICtab(alt.Tmod1, alt.Tmod2)
plot_model(alt.T90mod2,show.values = TRUE, value.offset = .2, dot.size = 3, transform="exp", vline.color = "black")

plot_model(alt.T90mod2.no0,show.values = TRUE, value.offset = .2, dot.size = 3, transform="exp", vline.color = "black")
plot_model(alt.T90mod2.96,show.values = TRUE, value.offset = .2, dot.size = 3, transform="exp", vline.color = "black")





exp(predict(alt.T90mod2, newdata = data.frame(lambda = "UVA",  temp = 10, hour.rescale = 0, TotalSpores = 1, experiment = "Exp10", cart="Cart1", base="b2")))
exp(predict(alt.T90mod2, newdata = data.frame(lambda = "UVA",  temp = 10, hour.rescale = 1, TotalSpores = 1, experiment = "Exp10", cart="Cart1", base="b2")))
exp(predict(alt.T90mod2, newdata = data.frame(lambda = "UVA",  temp = 10, hour.rescale = 2, TotalSpores = 1, experiment = "Exp10", cart="Cart1", base="b2")))
exp(predict(alt.T90mod2, newdata = data.frame(lambda = "UVA",  temp = 10, hour.rescale = 3, TotalSpores = 1, experiment = "Exp10", cart="Cart1", base="b2")))
exp(predict(alt.T90mod2, newdata = data.frame(lambda = "UVA",  temp = 10, hour.rescale = 4, TotalSpores = 1, experiment = "Exp10", cart="Cart1", base="b2")))

exp(predict(alt.T90mod2, newdata = data.frame(lambda = "UVA",  temp = 15, hour.rescale = 0, TotalSpores = 1, experiment = "Exp9", cart="Cart1", base="b2")))
exp(predict(alt.T90mod2, newdata = data.frame(lambda = "UVA",  temp = 15, hour.rescale = 1, TotalSpores = 1, experiment = "Exp9", cart="Cart1", base="b2")))
exp(predict(alt.T90mod2, newdata = data.frame(lambda = "UVA",  temp = 15, hour.rescale = 2, TotalSpores = 1, experiment = "Exp9", cart="Cart1", base="b2")))
exp(predict(alt.T90mod2, newdata = data.frame(lambda = "UVA",  temp = 15, hour.rescale = 3, TotalSpores = 1, experiment = "Exp9", cart="Cart1", base="b2")))
exp(predict(alt.T90mod2, newdata = data.frame(lambda = "UVA",  temp = 15, hour.rescale = 4, TotalSpores = 1, experiment = "Exp9", cart="Cart1", base="b2")))

exp(predict(alt.T90mod2, newdata = data.frame(lambda = "UVA",  temp = 20, hour.rescale = 0, TotalSpores = 1, experiment = "Exp8", cart="Cart1", base="b2")))
exp(predict(alt.T90mod2, newdata = data.frame(lambda = "UVA",  temp = 20, hour.rescale = 1, TotalSpores = 1, experiment = "Exp8", cart="Cart1", base="b2")))
exp(predict(alt.T90mod2, newdata = data.frame(lambda = "UVA",  temp = 20, hour.rescale = 2, TotalSpores = 1, experiment = "Exp8", cart="Cart1", base="b2")))
exp(predict(alt.T90mod2, newdata = data.frame(lambda = "UVA",  temp = 20, hour.rescale = 3, TotalSpores = 1, experiment = "Exp8", cart="Cart1", base="b2")))
exp(predict(alt.T90mod2, newdata = data.frame(lambda = "UVA",  temp = 20, hour.rescale = 4, TotalSpores = 1, experiment = "Exp8", cart="Cart1", base="b2")))

exp(predict(alt.T90mod2, newdata = data.frame(lambda = "UVA",  temp = 25, hour.rescale = 0, TotalSpores = 1, experiment = "Exp7", cart="Cart1", base="b2")))
exp(predict(alt.T90mod2, newdata = data.frame(lambda = "UVA",  temp = 25, hour.rescale = 1, TotalSpores = 1, experiment = "Exp7", cart="Cart1", base="b2")))
exp(predict(alt.T90mod2, newdata = data.frame(lambda = "UVA",  temp = 25, hour.rescale = 2, TotalSpores = 1, experiment = "Exp7", cart="Cart1", base="b2")))
exp(predict(alt.T90mod2, newdata = data.frame(lambda = "UVA",  temp = 25, hour.rescale = 3, TotalSpores = 1, experiment = "Exp7", cart="Cart1", base="b2")))
exp(predict(alt.T90mod2, newdata = data.frame(lambda = "UVA",  temp = 25, hour.rescale = 4, TotalSpores = 1, experiment = "Exp7", cart="Cart1", base="b2")))

# colnames(preds)<-"pred"
# 
# SODs<-rep(c(0,1,2,3,4), 4)
# temperature<-rep(c(10,15,20,25), each=5)
# 
# pred.mods<-cbind(SODs, temperature, preds)
# 
# ggplot(pred.mods, aes(x=SODs, y=pred, color=as.factor(temperature)))+geom_line()



# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

# SOL - RH

sol.togon01<-glmmTMB(AliveSpores~RH+temp+hour.rescale+lambda+(1|experiment/cart:base)+offset(log(TotalSpores)),
                     data=SOLsum[SOLsum$hour.rescale !=0,], family=nbinom2)

sol.togon02<-glmmTMB(AliveSpores~RH+temp+hour.rescale+(1|experiment/cart:base)+offset(log(TotalSpores)),
                     data=SOLsum[SOLsum$hour.rescale !=0,], family=nbinom2)
  
sol.togon03<-glmmTMB(AliveSpores~RH+temp+hour.rescale+as.factor(flux)+(1|experiment/cart:base)+offset(log(TotalSpores)),
                     data=SOLsum[SOLsum$hour.rescale !=0,], family=nbinom2)


  
sol.mod1<-glmmTMB(AliveSpores~RH+hour.rescale+(1|experiment/cart:base)+offset(log(TotalSpores)),
                  data=SOLsum, family=nbinom2)

sol.mod2<-glmmTMB(AliveSpores~RH+hour.rescale+lambda+(1|experiment/cart:base)+offset(log(TotalSpores)),
                  data=SOLsum, family=nbinom2)

# sol.mod3<-glmmTMB(AliveSpores~RH+hour.rescale+lambda+as.factor(flux)+(1|experiment/cart:base)+offset(log(TotalSpores)),
#                   data=SOLsum, family=nbinom2)

# sol.mod4<-glmmTMB(AliveSpores~RH+hour.rescale+lambda+base+(1|experiment/cart:base)+offset(log(TotalSpores)),
#                   data=SOLsum, family=nbinom2)

# sol.mod5<-glmmTMB(AliveSpores~RH+hour.rescale+lambda+base+as.factor(flux)+(1|experiment/cart:base)+offset(log(TotalSpores)),
#                     data=SOLsum, family=nbinom2)

sol.mod6<-glmmTMB(AliveSpores~RH:lambda+hour.rescale+lambda+(1|experiment/cart:base)+offset(log(TotalSpores)),
                  data=SOLsum, family=nbinom2)

sol.mod2.no0<-glmmTMB(AliveSpores~RH+hour.rescale+lambda+(1|experiment/cart:base)+offset(log(TotalSpores)),
                  data=SOLsum[SOLsum$hour.rescale != 0,], family=nbinom2)

sol.mod2.96<-glmmTMB(AliveSpores~RH+lambda+(1|experiment/cart:base)+offset(log(TotalSpores)),
                      data=SOLsum[SOLsum$hour.rescale == 4,], family=nbinom2)


soltest<-glmmTMB(AliveSpores~RH+temp+hour.rescale+lambda+(1|experiment/cart:base)+offset(log(TotalSpores)),
                  data=SOLsum, family=nbinom2)
plot_model(soltest,show.values = TRUE, value.offset = .2, dot.size = 3, transform="exp", vline.color = "black")

soltest.no0<-glmmTMB(AliveSpores~RH+temp+hour.rescale+lambda+(1|experiment/cart:base)+offset(log(TotalSpores)),
                 data=SOLsum[SOLsum$hour.rescale !=0,], family=nbinom2)
plot_model(soltest.no0,show.values = TRUE, value.offset = .2, dot.size = 3, transform="exp", vline.color = "black")


AICtab(sol.mod1,sol.mod2 )
plot_model(sol.mod2,show.values = TRUE, value.offset = .2, dot.size = 3, transform="exp", vline.color = "black")
plot_model(sol.mod6,show.values = TRUE, value.offset = .2, dot.size = 3, transform="exp", vline.color = "black")

plot_model(sol.mod2.no0,show.values = TRUE, value.offset = .2, dot.size = 3, transform="exp", vline.color = "black")
plot_model(sol.mod2.96,show.values = TRUE, value.offset = .2, dot.size = 3, transform="exp", vline.color = "black")


predicted <- predict(sol.mod2, newdata = data.frame(lambda = "UVA",  RH = 90, hour.rescale = 0, TotalSpores = 1, experiment = "Exp8", cart="Cart1", base="b4"))
predicted
exp(predicted)

# _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _

#SOL - T

sol.Tmod1<-glmmTMB(AliveSpores~temp+hour.rescale+(1|experiment/cart:base)+offset(log(TotalSpores)),
                   data=SOLsum, family=nbinom2)

sol.Tmod2<-glmmTMB(AliveSpores~temp+hour.rescale+lambda+(1|experiment/cart:base)+offset(log(TotalSpores)),
                   data=SOLsum, family=nbinom2)

sol.Tmod3<-glmmTMB(AliveSpores~temp+hour.rescale+lambda+as.factor(flux)+(1|experiment/cart:base)+offset(log(TotalSpores)),
                   data=SOLsum, family=nbinom2)

# sol.Tmod4<-glmmTMB(AliveSpores~temp+hour.rescale+lambda+base+(1|experiment/cart:base)+offset(log(TotalSpores)),
#                     data=SOLsum, family=nbinom2)

sol.Tmod5<-glmmTMB(AliveSpores~temp:lambda+hour.rescale+lambda+(1|experiment/cart:base)+offset(log(TotalSpores)),
                   data=SOLsum, family=nbinom2)

sol.Tmod2.no0<-glmmTMB(AliveSpores~temp+hour.rescale+lambda+(1|experiment/cart:base)+offset(log(TotalSpores)),
                   data=SOLsum[SOLsum$hour.rescale != 0,], family=nbinom2)

sol.Tmod2.n96<-glmmTMB(AliveSpores~temp+lambda+(1|experiment/cart:base)+offset(log(TotalSpores)),
                       data=SOLsum[SOLsum$hour.rescale == 4,], family=nbinom2)

AICtab(sol.Tmod1, sol.Tmod2, sol.Tmod3)
plot_model(sol.Tmod2,show.values = TRUE, value.offset = .2, dot.size = 3, transform="exp", vline.color = "black")
plot_model(sol.Tmod5,show.values = TRUE, value.offset = .2, dot.size = 3, transform="exp", vline.color = "black")

plot_model(sol.Tmod2.no0,show.values = TRUE, value.offset = .2, dot.size = 3, transform="exp", vline.color = "black")
plot_model(sol.Tmod2.n96,show.values = TRUE, value.offset = .2, dot.size = 3, transform="exp", vline.color = "black")

# _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _

#SOL - T sub 90
SOLsum90<-SOLsum[SOLsum$RH == "90",]

sol.T90mod1<-glmmTMB(AliveSpores~temp+hour.rescale+(1|experiment/cart:base)+offset(log(TotalSpores)),
                   data=SOLsum90, family=nbinom2)

sol.T90mod2<-glmmTMB(AliveSpores~temp+hour.rescale+lambda+(1|experiment/cart:base)+offset(log(TotalSpores)),
                   data=SOLsum90, family=nbinom2)

# sol.T90mod3<-glmmTMB(AliveSpores~temp+hour.rescale+lambda+as.factor(flux)+(1|experiment/cart:base)+offset(log(TotalSpores)),
#                   data=SOLsum90, family=nbinom1)
# 
# sol.T90mod4<-glmmTMB(AliveSpores~temp+hour.rescale+lambda+base+(1|experiment/cart:base)+offset(log(TotalSpores)),
#                    data=SOLsum90, family=nbinom1)

sol.T90mod5<-glmmTMB(AliveSpores~temp:lambda+hour.rescale+lambda+(1|experiment/cart:base)+offset(log(TotalSpores)),
                     data=SOLsum90, family=nbinom2)

sol.T90mod2.no0<-glmmTMB(AliveSpores~temp+hour.rescale+lambda+(1|experiment/cart:base)+offset(log(TotalSpores)),
                     data=SOLsum90[SOLsum90$hour.rescale != 0,], family=nbinom2)

sol.T90mod2.n96<-glmmTMB(AliveSpores~temp+lambda+(1|experiment/cart:base)+offset(log(TotalSpores)),
                         data=SOLsum90[SOLsum90$hour.rescale == 4,], family=nbinom2)

AICtab(sol.T90mod1, sol.T90mod2)
plot_model(sol.T90mod2,show.values = TRUE, value.offset = .2, dot.size = 3, transform="exp", vline.color = "black")
plot_model(sol.T90mod5,show.values = TRUE, value.offset = .2, dot.size = 3, transform="exp", vline.color = "black")

plot_model(sol.T90mod2.no0,show.values = TRUE, value.offset = .2, dot.size = 3, transform="exp", vline.color = "black")
plot_model(sol.T90mod2.n96,show.values = TRUE, value.offset = .2, dot.size = 3, transform="exp", vline.color = "black")


exp(predict(sol.T90mod2, newdata = data.frame(lambda = "UVA",  temp = 10, hour.rescale = 0, TotalSpores = 1, experiment = "Exp10", cart="Cart1", base="b2")))
exp(predict(sol.T90mod2, newdata = data.frame(lambda = "UVA",  temp = 10, hour.rescale = 1, TotalSpores = 1, experiment = "Exp10", cart="Cart1", base="b2")))
exp(predict(sol.T90mod2, newdata = data.frame(lambda = "UVA",  temp = 10, hour.rescale = 2, TotalSpores = 1, experiment = "Exp10", cart="Cart1", base="b2")))
exp(predict(sol.T90mod2, newdata = data.frame(lambda = "UVA",  temp = 10, hour.rescale = 3, TotalSpores = 1, experiment = "Exp10", cart="Cart1", base="b2")))
exp(predict(sol.T90mod2, newdata = data.frame(lambda = "UVA",  temp = 10, hour.rescale = 4, TotalSpores = 1, experiment = "Exp10", cart="Cart1", base="b2")))

exp(predict(sol.T90mod2, newdata = data.frame(lambda = "UVA",  temp = 15, hour.rescale = 0, TotalSpores = 1, experiment = "Exp9", cart="Cart1", base="b2")))
exp(predict(sol.T90mod2, newdata = data.frame(lambda = "UVA",  temp = 15, hour.rescale = 1, TotalSpores = 1, experiment = "Exp9", cart="Cart1", base="b2")))
exp(predict(sol.T90mod2, newdata = data.frame(lambda = "UVA",  temp = 15, hour.rescale = 2, TotalSpores = 1, experiment = "Exp9", cart="Cart1", base="b2")))
exp(predict(sol.T90mod2, newdata = data.frame(lambda = "UVA",  temp = 15, hour.rescale = 3, TotalSpores = 1, experiment = "Exp9", cart="Cart1", base="b2")))
exp(predict(sol.T90mod2, newdata = data.frame(lambda = "UVA",  temp = 15, hour.rescale = 4, TotalSpores = 1, experiment = "Exp9", cart="Cart1", base="b2")))

exp(predict(sol.T90mod2, newdata = data.frame(lambda = "UVA",  temp = 20, hour.rescale = 0, TotalSpores = 1, experiment = "Exp8", cart="Cart1", base="b2")))
exp(predict(sol.T90mod2, newdata = data.frame(lambda = "UVA",  temp = 20, hour.rescale = 1, TotalSpores = 1, experiment = "Exp8", cart="Cart1", base="b2")))
exp(predict(sol.T90mod2, newdata = data.frame(lambda = "UVA",  temp = 20, hour.rescale = 2, TotalSpores = 1, experiment = "Exp8", cart="Cart1", base="b2")))
exp(predict(sol.T90mod2, newdata = data.frame(lambda = "UVA",  temp = 20, hour.rescale = 3, TotalSpores = 1, experiment = "Exp8", cart="Cart1", base="b2")))
exp(predict(sol.T90mod2, newdata = data.frame(lambda = "UVA",  temp = 20, hour.rescale = 4, TotalSpores = 1, experiment = "Exp8", cart="Cart1", base="b2")))

exp(predict(sol.T90mod2, newdata = data.frame(lambda = "UVA",  temp = 25, hour.rescale = 0, TotalSpores = 1, experiment = "Exp7", cart="Cart1", base="b2")))
exp(predict(sol.T90mod2, newdata = data.frame(lambda = "UVA",  temp = 25, hour.rescale = 1, TotalSpores = 1, experiment = "Exp7", cart="Cart1", base="b2")))
exp(predict(sol.T90mod2, newdata = data.frame(lambda = "UVA",  temp = 25, hour.rescale = 2, TotalSpores = 1, experiment = "Exp7", cart="Cart1", base="b2")))
exp(predict(sol.T90mod2, newdata = data.frame(lambda = "UVA",  temp = 25, hour.rescale = 3, TotalSpores = 1, experiment = "Exp7", cart="Cart1", base="b2")))
exp(predict(sol.T90mod2, newdata = data.frame(lambda = "UVA",  temp = 25, hour.rescale = 4, TotalSpores = 1, experiment = "Exp7", cart="Cart1", base="b2")))

####################################################################################
# grid.arrange best models RH & T
library(gridExtra)
library(sjPlot)
library(ggplot2)
grid.arrange(
plot_model(alt.mod2,show.values = TRUE, value.offset = .2, dot.size = 3, transform="exp", vline.color = "black")+theme_bw()+ggtitle("ALT:RH"),
plot_model(alt.T90mod2,show.values = TRUE, value.offset = .2, dot.size = 3, transform="exp", vline.color = "black")+theme_bw()+ggtitle("ALT:T90"),
plot_model(sol.mod2,show.values = TRUE, value.offset = .2, dot.size = 3, transform="exp", vline.color = "black")+theme_bw()+ggtitle("SOL:RH"),
plot_model(sol.T90mod2,show.values = TRUE, value.offset = .2, dot.size = 3, transform="exp", vline.color = "black")+theme_bw()+ggtitle("SOL:T90"),
nrow=2, ncol=2)

grid.arrange(
  plot_model(alt.togon01,show.values = TRUE, value.offset = .4, dot.size = 3, transform="exp", vline.color = "black")+theme_bw()+ggtitle("ALT:24-96"),
  plot_model(sol.togon01,show.values = TRUE, value.offset = .4, dot.size = 3, transform="exp", vline.color = "black")+theme_bw()+ggtitle("SOL:24-96"),
  plot_model(alt.T90mod2.no0,show.values = TRUE, value.offset = .4, dot.size = 3, transform="exp", vline.color = "black")+theme_bw()+ggtitle("ALT:T90; 24-96"),
  plot_model(sol.T90mod2.no0,show.values = TRUE, value.offset = .4, dot.size = 3, transform="exp", vline.color = "black")+theme_bw()+ggtitle("SOL:T90; 24-96"),
nrow=2, ncol=2)
############################################################################################################
#
#                               POST_HOC
#
############################################################################################################
library(splines)
library(mgcv)
library(glmmTMB)
library(dplyr)
library(ggplot2)
library(multcomp)
glht_glmmTMB <- function (model, ..., component="cond") {
  glht(model, ...,
       coef. = function(x) fixef(x)[[component]],
       vcov. = function(x) vcov(x)[[component]],
       df = NULL)
}

modelparm.glmmTMB <- function (model, coef. = function(x) fixef(x)[[component]],
                               vcov. = function(x) vcov(x)[[component]],
                               df = NULL, component="cond", ...) {
                      multcomp:::modelparm.default(model, coef. = coef., vcov. = vcov.,
                               df = df, ...)
}

# _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
#ALT
glht_ALT.RH <- glht_glmmTMB(alt.mod2, linfct = mcp(RH = "Tukey"))
summary(glht_ALT.RH)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)   
# 60 - 50 == 0 -0.069174   0.097907  -0.707  0.89180   
# 75 - 50 == 0 -0.075575   0.072338  -1.045  0.71719   
# 90 - 50 == 0  0.149546   0.068430   2.185  0.12325   
# 75 - 60 == 0 -0.006401   0.092792  -0.069  0.99988   
# 90 - 60 == 0  0.218720   0.089773   2.436  0.06791 . 
# 90 - 75 == 0  0.225121   0.060900   3.697  0.00127 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Adjusted p values reported -- single-step method)

glht_ALT.T90 <- glht_glmmTMB(alt.T90mod2, linfct = mcp(temp = "Tukey"))
summary(glht_ALT.T90)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)    
# 15 - 10 == 0  0.07301    0.03420   2.135   0.1419    
# 20 - 10 == 0  0.08667    0.03682   2.354   0.0859 .  
# 25 - 10 == 0 -0.22657    0.03624  -6.252   <0.001 ***
# 20 - 15 == 0  0.01366    0.03571   0.383   0.9809    
# 25 - 15 == 0 -0.29958    0.03642  -8.226   <0.001 ***
# 25 - 20 == 0 -0.31324    0.03884  -8.066   <0.001 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Adjusted p values reported -- single-step method)

# * * * * * * no0

glht_ALT.RH.no0 <- glht_glmmTMB(alt.mod2no0, linfct = mcp(RH = "Tukey"))
summary(glht_ALT.RH.no0)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)   
# 60 - 50 == 0  -0.3577     0.1880  -1.903   0.2209   
# 75 - 50 == 0  -0.2077     0.1397  -1.487   0.4379   
# 90 - 50 == 0   0.1620     0.1324   1.224   0.6047   
# 75 - 60 == 0   0.1500     0.1774   0.846   0.8286   
# 90 - 60 == 0   0.5197     0.1718   3.025   0.0127 * 
# 90 - 75 == 0   0.3697     0.1170   3.160   0.0082 **


glht_ALT.T90.no0 <- glht_glmmTMB(alt.T90mod2.no0, linfct = mcp(temp = "Tukey"))
summary(glht_ALT.T90.no0)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)    
# 15 - 10 == 0  0.18788    0.04611   4.075   <0.001 ***
# 20 - 10 == 0  0.09412    0.04738   1.987    0.193    
# 25 - 10 == 0 -0.24829    0.04682  -5.303   <0.001 ***
# 20 - 15 == 0 -0.09376    0.04683  -2.002    0.187    
# 25 - 15 == 0 -0.43617    0.04636  -9.408   <0.001 ***
# 25 - 20 == 0 -0.34241    0.04761  -7.192   <0.001 ***



# _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
#SOL

glht_SOL.RH <- glht_glmmTMB(sol.mod2, linfct = mcp(RH = "Tukey"))
summary(glht_SOL.RH)

# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)    
# 60 - 50 == 0 -0.05970    0.10835  -0.551    0.945    
# 75 - 50 == 0 -0.16744    0.08064  -2.076    0.156    
# 90 - 50 == 0  0.39114    0.07646   5.116   <0.001 ***
# 75 - 60 == 0 -0.10774    0.10223  -1.054    0.712    
# 90 - 60 == 0  0.45084    0.09898   4.555   <0.001 ***
# 90 - 75 == 0  0.55858    0.06753   8.272   <0.001 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Adjusted p values reported -- single-step method)

glht_SOL.T90 <- glht_glmmTMB(sol.T90mod2, linfct = mcp(temp = "Tukey"))
summary(glht_SOL.T90)

# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)    
# 15 - 10 == 0  0.06030    0.03891   1.550    0.407    
# 20 - 10 == 0  0.24222    0.03707   6.534   <0.001 ***
# 25 - 10 == 0  0.14568    0.03724   3.912   <0.001 ***
# 20 - 15 == 0  0.18192    0.03732   4.874   <0.001 ***
# 25 - 15 == 0  0.08538    0.03742   2.282    0.102    
# 25 - 20 == 0 -0.09653    0.03556  -2.714    0.033 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Adjusted p values reported -- single-step method)

# * * * * * * no0

glht_SOL.RH.no0 <- glht_glmmTMB(sol.mod2.no0, linfct = mcp(RH = "Tukey"))
summary(glht_SOL.RH.no0)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)    
# 60 - 50 == 0 -0.13994    0.09456  -1.480    0.442    
# 75 - 50 == 0 -0.11588    0.07017  -1.652    0.342    
# 90 - 50 == 0  0.53530    0.06645   8.055   <1e-04 ***
# 75 - 60 == 0  0.02405    0.08930   0.269    0.993    
# 90 - 60 == 0  0.67524    0.08642   7.814   <1e-04 ***
# 90 - 75 == 0  0.65119    0.05874  11.087   <1e-04 ***

glht_SOL.T90.no0 <- glht_glmmTMB(sol.T90mod2.no0, linfct = mcp(temp = "Tukey"))
summary(glht_SOL.T90.no0)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)    
# 15 - 10 == 0  0.12779    0.06864   1.862    0.245    
# 20 - 10 == 0  0.27585    0.06840   4.033   <0.001 ***
# 25 - 10 == 0  0.15045    0.06838   2.200    0.123    
# 20 - 15 == 0  0.14806    0.06785   2.182    0.128    
# 25 - 15 == 0  0.02266    0.06783   0.334    0.987    
# 25 - 20 == 0 -0.12541    0.06758  -1.856    0.248    


# * * * * * * * * * * * * * * * * * * * * * * * * * *

multalttogo.T <- glht_glmmTMB(alt.togon01, linfct = mcp(temp = "Tukey"))
summary(multalttogo.T)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)   
# 15 - 10 == 0  0.01996    0.08823   0.226  0.99581   
# 20 - 10 == 0 -0.16494    0.08856  -1.862  0.23842   
# 25 - 10 == 0 -0.39412    0.11921  -3.306  0.00503 **
# 20 - 15 == 0 -0.18490    0.07554  -2.448  0.06613 . 
# 25 - 15 == 0 -0.41407    0.11613  -3.566  0.00196 **
# 25 - 20 == 0 -0.22918    0.11635  -1.970  0.19434  
multalttogo.RH <- glht_glmmTMB(alt.togon01, linfct = mcp(RH = "Tukey"))
summary(multalttogo.RH)
#Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)    
# 60 - 50 == 0 -0.42994    0.13930  -3.086   0.0103 *  
# 75 - 50 == 0 -0.23080    0.08855  -2.606   0.0425 *  
# 90 - 50 == 0  0.22411    0.08823   2.540   0.0506 .  
# 75 - 60 == 0  0.19914    0.12007   1.659   0.3338    
# 90 - 60 == 0  0.65405    0.12023   5.440   <0.001 ***
# 90 - 75 == 0  0.45491    0.07556   6.020   <0.001 ***
multalttogo.lambda <- glht_glmmTMB(alt.togon01, linfct = mcp(lambda = "Tukey"))
summary(multalttogo.lambda)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)    
# UVA - b21 == 0 -0.12105    0.05935   -2.04    0.095 .  
# UVB - b21 == 0 -1.80213    0.05998  -30.05   <1e-04 ***
# UVB - UVA == 0 -1.68107    0.02683  -62.65   <1e-04 ***

# _______________

multsoltogo.T <- glht_glmmTMB(sol.togon01, linfct = mcp(temp = "Tukey"))
summary(multsoltogo.T)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)
# 15 - 10 == 0  0.129105   0.060591   2.131    0.139
# 20 - 10 == 0  0.097507   0.060750   1.605    0.369
# 25 - 10 == 0  0.094246   0.081555   1.156    0.649
# 20 - 15 == 0 -0.031598   0.051693  -0.611    0.927
# 25 - 15 == 0 -0.034859   0.079363  -0.439    0.971
# 25 - 20 == 0 -0.003261   0.079382  -0.041    1.000
multsoltogo.RH <- glht_glmmTMB(sol.togon01, linfct = mcp(RH = "Tukey"))
summary(multsoltogo.RH)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)    
# 60 - 50 == 0 -0.02663    0.09544  -0.279    0.992    
# 75 - 50 == 0 -0.07780    0.06068  -1.282    0.562    
# 90 - 50 == 0  0.56834    0.06048   9.398   <1e-04 ***
# 75 - 60 == 0 -0.05117    0.08232  -0.622    0.922    
# 90 - 60 == 0  0.59497    0.08234   7.226   <1e-04 ***
# 90 - 75 == 0  0.64614    0.05175  12.487   <1e-04 ***
multsoltogo.lambda <- glht_glmmTMB(sol.togon01, linfct = mcp(lambda = "Tukey"))
summary(multsoltogo.lambda)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)    
# UVA - b21 == 0 -0.18607    0.08076  -2.304   0.0502 .  
# UVB - b21 == 0 -0.99040    0.08084 -12.251   <1e-04 ***
# UVB - UVA == 0 -0.80433    0.03470 -23.182   <1e-04 ***  

# best 90 sub, no0 _ _ _  _ _ _ _ _  _ _ _ _ ____ _ _ _
mult90.alt <- glht_glmmTMB(alt.T90mod2.no0, linfct = mcp(temp = "Tukey"))
summary(mult90.alt)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)    
# 15 - 10 == 0  0.18788    0.04611   4.075   <0.001 ***
# 20 - 10 == 0  0.09412    0.04738   1.987    0.193    
# 25 - 10 == 0 -0.24829    0.04682  -5.303   <0.001 ***
# 20 - 15 == 0 -0.09376    0.04683  -2.002    0.187    
# 25 - 15 == 0 -0.43617    0.04636  -9.408   <0.001 ***
# 25 - 20 == 0 -0.34241    0.04761  -7.192   <0.001 ***
mult90lam.alt <- glht_glmmTMB(alt.T90mod2.no0, linfct = mcp(lambda = "Tukey"))
summary(mult90lam.alt)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)    
# UVA - b21 == 0 -0.02910    0.05365  -0.542    0.844    
# UVB - b21 == 0 -1.03506    0.05579 -18.553   <1e-05 ***
# UVB - UVA == 0 -1.00596    0.02764 -36.398   <1e-05 ***

#    #    #    #    #    #    #
mult90.sol <- glht_glmmTMB(sol.T90mod2.no0, linfct = mcp(temp = "Tukey"))
summary(mult90.sol)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)    
# 15 - 10 == 0  0.13434    0.07165   1.875    0.239    
# 20 - 10 == 0  0.28740    0.07154   4.017   <0.001 ***
# 25 - 10 == 0  0.15756    0.07155   2.202    0.123    
# 20 - 15 == 0  0.15305    0.07129   2.147    0.138    
# 25 - 15 == 0  0.02321    0.07129   0.326    0.988    
# 25 - 20 == 0 -0.12984    0.07118  -1.824    0.262  

mult90lam.sol <- glht_glmmTMB(sol.T90mod2.no0, linfct = mcp(lambda = "Tukey"))
summary(mult90lam.sol)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)    
# UVA - b21 == 0 -0.10883    0.12078  -0.901    0.625    
# UVB - b21 == 0 -1.09296    0.12091  -9.040   <1e-05 ***
# UVB - UVA == 0 -0.98414    0.05176 -19.012   <1e-05 ***

########################################################################################################
########################################################################################################
########################################################################################################
##############################        TRY CROSS 1                    ###################################
########################################################################################################
########################################################################################################
########################################################################################################


#need to change RE to cart:base only re:Ané

#ALT
ALTsum.cr1<-ALTsum[ALTsum$experiment %in% c("Exp4", "Exp5", "Exp6", "Exp8", "Exp9", "Exp10"),]

alt.crmod1<-glmmTMB(AliveSpores~(RH*temp)+hour.rescale+(1|cart:base)+offset(log(TotalSpores)),
                  data=ALTsum.cr1, family=nbinom1)

alt.crmod2<-glmmTMB(AliveSpores~(RH*temp)+hour.rescale+lambda+(1|cart:base)+offset(log(TotalSpores)),
                  data=ALTsum.cr1, family=nbinom1)

alt.crmod2.no0<-glmmTMB(AliveSpores~(RH*temp)+hour.rescale+lambda+(1|cart:base)+offset(log(TotalSpores)),
                    data=ALTsum.cr1[ALTsum.cr1$hour.rescale !=0,], family=nbinom1)

AICtab(alt.crmod1, alt.crmod2)
plot_model(alt.crmod2,show.values = TRUE, value.offset = .2, dot.size = 3, transform="exp", vline.color = "black")


#SOL
SOLsum.cr1<-SOLsum[SOLsum$experiment %in% c("Exp4", "Exp5", "Exp6", "Exp8", "Exp9", "Exp10"),]

sol.crmod1<-glmmTMB(AliveSpores~(RH*temp)+hour.rescale+(1|cart:base)+offset(log(TotalSpores)),
                    data=SOLsum.cr1, family=nbinom2)

sol.crmod2<-glmmTMB(AliveSpores~(RH*temp)+hour.rescale+lambda+(1|cart:base)+offset(log(TotalSpores)),
                    data=SOLsum.cr1, family=nbinom2)

sol.crmod2.no0<-glmmTMB(AliveSpores~(RH*temp)+hour.rescale+lambda+(1|cart:base)+offset(log(TotalSpores)),
                    data=SOLsum.cr1[SOLsum.cr1$hour.rescale !=0,], family=nbinom2)

AICtab(sol.crmod1, sol.crmod2)
plot_model(sol.crmod2,show.values = TRUE, value.offset = .2, dot.size = 3, transform="exp", vline.color = "black")

grid.arrange(
  plot_model(alt.crmod2.no0,show.values = TRUE, value.offset = .4, dot.size = 3, transform="exp", vline.color = "black")+theme_bw()+ggtitle("CrALT:24-96"),
  plot_model(sol.crmod2.no0,show.values = TRUE, value.offset = .4, dot.size = 3, transform="exp", vline.color = "black")+theme_bw()+ggtitle("CrSOL:24-96"),
  nrow=1, ncol=2)

################################################################################
#Piece-wise regression
#Y = b0+b1x + b2(x-xk)xk


ALTsump<-ALTsum

dummy<-rep(0, length(ALTsump$AliveSpores))
dummy[ALTsump$hour.rescale > 1] = 1

ALTsump$Xdif=ALTsump$hour.rescale-1
ALTsump$DM<-dummy

ALTsump$X<-ALTsump$Xdif*ALTsump$DM

#RH
alt.mod2PW<-glmmTMB(AliveSpores~RH+hour.rescale   + X      +lambda+(1|experiment/cart:base)+offset(log(TotalSpores)),
                  data=ALTsump, family=nbinom1)
summary(alt.mod2PW)
plot_model(alt.mod2PW,show.values = TRUE, value.offset = .2, dot.size = 3, transform="exp", vline.color = "black")

# library(lme4)
# glmer.test<-glmer.nb(AliveSpores~RH+hour.rescale +lambda+(1|experiment/cart:base)+offset(log(TotalSpores)),
#                      data=ALTsum, nAGQ=0,
#                      control=glmerControl(optimizer = "nloptwrap"))
# summary(glmer.test)
# 
# glmer.test2<-glmer.nb(AliveSpores~RH+hour.rescale+ X +lambda+(1|experiment/cart:base)+offset(log(TotalSpores)),
#                      data=ALTsump, nAGQ=0,
#                      control=glmerControl(optimizer = "nloptwrap"))
# summary(glmer.test2)
# 
# estimates_glmer.test2<-as.data.frame(summary(glmer.test2)$coef)
# estimates_glmer.test2$exp_trans<-exp(estimates_glmer.test2$Estimate)
# plot_model(glmer.test2,show.values = TRUE, value.offset = .2, dot.size = 3, transform="exp", vline.color = "black")
# 
# 
# AICtab(alt.mod2PW, alt.mod2) # so basically it slightly improves the models


#T
alt.Tmod2PW<-glmmTMB(AliveSpores~temp+hour.rescale   + X      +lambda+(1|experiment/cart:base)+offset(log(TotalSpores)),
                    data=ALTsump, family=nbinom1)
summary(alt.Tmod2PW)
plot_model(alt.Tmod2PW,show.values = TRUE, value.offset = .2, dot.size = 3, transform="exp", vline.color = "black")

#sub90
alt.sub90mod2PW<-glmmTMB(AliveSpores~temp+hour.rescale   + X      +lambda+(1|experiment/cart:base)+offset(log(TotalSpores)),
                     data=ALTsump[ALTsump$RH == 90,], family=nbinom1)

summary(alt.sub90mod2PW)

plot_model(alt.sub90mod2PW,show.values = TRUE, value.offset = .2, dot.size = 3, transform="exp", vline.color = "black")

glht_ALT.sub90PW <- glht_glmmTMB(alt.sub90mod2PW, linfct = mcp(temp = "Tukey"))
summary(glht_ALT.sub90PW)

# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)    
# 15 - 10 == 0  0.04932    0.03413   1.445    0.471    
# 20 - 10 == 0  0.06507    0.03663   1.776    0.284    
# 25 - 10 == 0 -0.22618    0.03590  -6.301   <1e-04 ***
# 20 - 15 == 0  0.01574    0.03545   0.444    0.971    
# 25 - 15 == 0 -0.27551    0.03616  -7.620   <1e-04 ***
# 25 - 20 == 0 -0.29125    0.03850  -7.565   <1e-04 ***

# ** * *  * * * * * *** *** *  ** * * ***   ** ** * ** *  **** * ** * * **  *** 

SOLsump<-SOLsum

dummys<-rep(0, length(SOLsump$AliveSpores))
dummys[SOLsump$hour.rescale > 1] = 1

SOLsump$Xdif=SOLsump$hour.rescale-1
SOLsump$DM<-dummys

SOLsump$X<-SOLsump$Xdif*SOLsump$DM

#RH
sol.mod2PW<-glmmTMB(AliveSpores~RH+hour.rescale   + X      +lambda+(1|experiment/cart:base)+offset(log(TotalSpores)),
                    data=SOLsump, family=nbinom2)
summary(sol.mod2PW)
plot_model(sol.mod2PW,show.values = TRUE, value.offset = .2, dot.size = 3, transform="exp", vline.color = "black")


-0.176031
0.008299
-0.131923
0.245930
-0.86453

-0.065923
-0.464598

og<-(exp(-0.176031)*3)+(exp(0.008299)*3)+(exp(-0.131923)*3)+(exp(0.245930)*3)+(exp(-0.86453)*3)+(exp(-0.065923)*3)+(exp(-0.464598)*3)

PW<-exp(0.790555)*((3-1)*1)


#T
sol.Tmod2PW<-glmmTMB(AliveSpores~temp+hour.rescale   + X      +lambda+(1|experiment/cart:base)+offset(log(TotalSpores)),
                    data=SOLsump, family=nbinom2)
summary(sol.Tmod2PW)
plot_model(sol.Tmod2PW,show.values = TRUE, value.offset = .2, dot.size = 3, transform="exp", vline.color = "black")

#90 sub
sol.sub90mod2PW<-glmmTMB(AliveSpores~temp+hour.rescale   + X      +lambda+(1|experiment/cart:base)+offset(log(TotalSpores)),
                     data=SOLsump[SOLsump$RH ==90,], family=nbinom2)
sol.sub90mod2PW2<-glmmTMB(AliveSpores~temp:lambda+hour.rescale   + X      +lambda+(1|experiment/cart:base)+offset(log(TotalSpores)),
                         data=SOLsump[SOLsump$RH ==90,], family=nbinom2)

summary(sol.sub90mod2PW)
plot_model(sol.sub90mod2PW,show.values = TRUE, value.offset = .2, dot.size = 3, transform="exp", vline.color = "black")
plot_model(sol.sub90mod2PW2,show.values = TRUE, value.offset = .2, dot.size = 3, transform="exp", vline.color = "black")



glht_SOL.sub90PW <- glht_glmmTMB(sol.sub90mod2PW, linfct = mcp(temp = "Tukey"))
summary(glht_SOL.sub90PW)
# 15 - 10 == 0  0.10784    0.04870   2.214  0.11927    
# 20 - 10 == 0  0.28122    0.04856   5.792  < 0.001 ***
# 25 - 10 == 0  0.18215    0.04856   3.751  0.00102 ** 
# 20 - 15 == 0  0.17338    0.04843   3.580  0.00191 ** 
# 25 - 15 == 0  0.07430    0.04844   1.534  0.41704    
# 25 - 20 == 0 -0.09908    0.04829  -2.052  0.16928    

################# Compile

grid.arrange(
  
  plot_model(alt.mod2PW,show.values = TRUE, value.offset = .2, dot.size = 3, transform="exp", vline.color = "black")+theme_bw()+ggtitle("PW: ALT RH"),
  plot_model(sol.mod2PW,show.values = TRUE, value.offset = .2, dot.size = 3, transform="exp", vline.color = "black")+theme_bw()+ggtitle("PW: SOL RH"),
  
  plot_model(alt.sub90mod2PW,show.values = TRUE, value.offset = .2, dot.size = 3, transform="exp", vline.color = "black")+theme_bw()+ggtitle("PW: ALT T90"),
  plot_model(sol.sub90mod2PW,show.values = TRUE, value.offset = .2, dot.size = 3, transform="exp", vline.color = "black")+theme_bw()+ggtitle("PW: SOL T90"),
  
  nrow=2, ncol=2
)



estimates_alt.mod2PW<-as.data.frame(summary(alt.mod2PW)$coef$cond)
estimates_alt.mod2PW$exp_trans<-exp(estimates_alt.mod2PW$Estimate)
estimates_alt.mod2PW




grid.arrange(
  plot_model(alt.mod2,show.values = TRUE, value.offset = .2, dot.size = 3, transform="exp", vline.color = "black")+theme_bw()+ggtitle("ALT:RH"),
  plot_model(sol.mod2,show.values = TRUE, value.offset = .2, dot.size = 3, transform="exp", vline.color = "black")+theme_bw()+ggtitle("SOL:RH"),
  plot_model(alt.T90mod2,show.values = TRUE, value.offset = .2, dot.size = 3, transform="exp", vline.color = "black")+theme_bw()+ggtitle("ALT:T90"),
  plot_model(sol.T90mod2,show.values = TRUE, value.offset = .2, dot.size = 3, transform="exp", vline.color = "black")+theme_bw()+ggtitle("SOL:T90"),
  nrow=2, ncol=2)




############################################################################################
############################################################################################
#                                   Krskall-Wallis; post hoc Dunn test
############################################################################################
############################################################################################
#TukeyHSD(ALTsumUVArez, which = "temp")
library(multcomp)
library(FSA)
library(dunn.test)

#By lambda
ALTsumLAMBDArez<-aov(AliveSpores ~ lambda + hour.rescale, data = ALTsum)
summary(glht(ALTsumLAMBDArez, linfct = mcp(lambda = "Tukey"))) #takes unbalanced design into consideration

hist(log(ALTsum[ALTsum$lambda == "UVA",]$AliveSpores))
hist(log(ALTsum[ALTsum$lambda == "UVB",]$AliveSpores))
hist(log(ALTsum[ALTsum$lambda == "b21",]$AliveSpores))

ALTsumno0<-ALTsum[ALTsum$hour.rescale != 0,]
hist((ALTsumno0[ALTsumno0$lambda == "UVA",]$AliveSpores)/(ALTsumno0[ALTsumno0$lambda == "UVA",]$TotalSpores))
hist((ALTsumno0[ALTsumno0$lambda == "UVB",]$AliveSpores)/(ALTsumno0[ALTsumno0$lambda == "UVB",]$TotalSpores))
hist((ALTsumno0[ALTsumno0$lambda == "b21",]$AliveSpores)/(ALTsumno0[ALTsumno0$lambda == "b21",]$TotalSpores))

#box plot made with ALTlimTex and SOLlimTex, so use these
ALTlimTex.no<-ALTlimTex[ALTlimTex$hour.rescale !=0,]
hist((ALTlimTex.no[ALTlimTex.no$lambda == "UVA",]$mean.T))
hist((ALTlimTex.no[ALTlimTex.no$lambda == "UVB",]$mean.T))
hist((ALTlimTex.no[ALTlimTex.no$lambda == "b21",]$mean.T))

# "bh" the FDR is controlled using the Benjamini-Hochberg adjustment (1995), 
# a step-down procedure appropriate to independent tests or tests that are 
# positively dependent. p-values are ordered from largest to smallest,
# and adjusted p-values = max[1, pm/(m+1-i)], where i indexes the ordering.
# All tests after and including the first to be rejected are also rejected.

#two-sided tells me if group 1 is different from group 2
#one sided tells me if group 1 is higher/lower than group 2
#because two-sided, sig p-values is alpha/2, where alpha is my p-value cutoff

#BH p.adjust prob fine: http://www.biostathandbook.com/multiplecomparisons.html

################### ALT

# **limTEX are made on line ~126

#BY lambda
kruskal.test(mean.T ~ lambda,
             data = ALTlimTex.no)
dunn.test(ALTlimTex.no$mean.T, ALTlimTex.no$lambda, method="bh", list=TRUE)

#BY T
kruskal.test(mean.T ~ temp,
             data = ALTlimTex.no[ALTlimTex.no$lambda == "UVA",])
table<-dunn.test(ALTlimTex.no[ALTlimTex.no$lambda == "UVA",]$mean.T, ALTlimTex.no[ALTlimTex.no$lambda == "UVA",]$temp, method="bh", list=TRUE)

kruskal.test(mean.T ~ temp,
             data = ALTlimTex.no[ALTlimTex.no$lambda == "UVB",])
dunn.test(ALTlimTex.no[ALTlimTex.no$lambda == "UVB",]$mean.T, ALTlimTex.no[ALTlimTex.no$lambda == "UVB",]$temp, method="bh", list=TRUE)

kruskal.test(mean.T ~ temp,
             data = ALTlimTex.no[ALTlimTex.no$lambda == "b21",])
dunn.test(ALTlimTex.no[ALTlimTex.no$lambda == "b21",]$mean.T, ALTlimTex.no[ALTlimTex.no$lambda == "b21",]$temp, method="bh", list=TRUE)

#BY RH
kruskal.test(mean.T ~ RH,
             data = ALTlimTex.no[ALTlimTex.no$lambda == "UVA",])
dunn.test(ALTlimTex.no[ALTlimTex.no$lambda == "UVA",]$mean.T, ALTlimTex.no[ALTlimTex.no$lambda == "UVA",]$RH, method="bh", list=TRUE)

kruskal.test(mean.T ~ RH,
             data = ALTlimTex.no[ALTlimTex.no$lambda == "UVB",])
dunn.test(ALTlimTex.no[ALTlimTex.no$lambda == "UVB",]$mean.T, ALTlimTex.no[ALTlimTex.no$lambda == "UVB",]$RH, method="bh", list=TRUE)

kruskal.test(mean.T ~ RH,
             data = ALTlimTex.no[ALTlimTex.no$lambda == "b21",])
dunn.test(ALTlimTex.no[ALTlimTex.no$lambda == "b21",]$mean.T, ALTlimTex.no[ALTlimTex.no$lambda == "b21",]$RH, method="bh", list=TRUE)


#BY experiment
kruskal.test(mean.T ~ experiment,
             data = ALTlimTex.no[ALTlimTex.no$lambda == "UVA",])
dunn.test(ALTlimTex.no[ALTlimTex.no$lambda == "UVA",]$mean.T, ALTlimTex.no[ALTlimTex.no$lambda == "UVA",]$experiment, method="bh", list=TRUE)

kruskal.test(mean.T ~ experiment,
             data = ALTlimTex.no[ALTlimTex.no$lambda == "UVB",])
dunn.test(ALTlimTex.no[ALTlimTex.no$lambda == "UVB",]$mean.T, ALTlimTex.no[ALTlimTex.no$lambda == "UVB",]$experiment, method="bh", list=TRUE)

kruskal.test(mean.T ~ experiment,
             data = ALTlimTex.no[ALTlimTex.no$lambda == "b21",])
dunn.test(ALTlimTex.no[ALTlimTex.no$lambda == "b21",]$mean.T, ALTlimTex.no[ALTlimTex.no$lambda == "b21",]$experiment, method="bh", list=TRUE)

#BY experiment just 90

# Day 3
kruskal.test(mean.T ~ experiment,
             data = ALTlimTexbase[ALTlimTexbase$lambda %in% c("UVA") & ALTlimTexbase$hour.rescale %in% c(3) & ALTlimTexbase$RH %in% c(90),])
dunn.test(ALTlimTexbase[ALTlimTexbase$lambda %in% c("UVA") & ALTlimTexbase$hour.rescale %in% c(3) & ALTlimTexbase$RH %in% c(90),]$mean.T, 
          ALTlimTexbase[ALTlimTexbase$lambda %in% c("UVA") & ALTlimTexbase$hour.rescale %in% c(3) & ALTlimTexbase$RH %in% c(90),]$experiment, 
          method="bh", list=TRUE)

#Day 4
kruskal.test(mean.T ~ experiment,
             data = ALTlimTexbase[ALTlimTexbase$lambda %in% c("UVA") & ALTlimTexbase$hour.rescale %in% c(4) & ALTlimTexbase$RH %in% c(90),])
dunn.test(ALTlimTexbase[ALTlimTexbase$lambda %in% c("UVA") & ALTlimTexbase$hour.rescale %in% c(4) & ALTlimTexbase$RH %in% c(90),]$mean.T, 
          ALTlimTexbase[ALTlimTexbase$lambda %in% c("UVA") & ALTlimTexbase$hour.rescale %in% c(4) & ALTlimTexbase$RH %in% c(90),]$experiment, 
          method="bh", list=TRUE)


#################SOL
SOLlimTex.no<-SOLlimTex[SOLlimTex$hour.rescale !=0,]
hist((SOLlimTex.no[SOLlimTex.no$lambda == "UVA",]$mean.T))
hist((SOLlimTex.no[SOLlimTex.no$lambda == "UVB",]$mean.T))
hist((SOLlimTex.no[SOLlimTex.no$lambda == "b21",]$mean.T))

#BY lambda
kruskal.test(mean.T ~ lambda,
             data = SOLlimTex.no)
dunn.test(SOLlimTex.no$mean.T, SOLlimTex.no$lambda, method="bh", list=TRUE)

#BY T
kruskal.test(mean.T ~ temp,
             data = SOLlimTex.no[SOLlimTex.no$lambda == "UVA",])
dunn.test(SOLlimTex.no[SOLlimTex.no$lambda == "UVA",]$mean.T, SOLlimTex.no[SOLlimTex.no$lambda == "UVA",]$temp, method="bh", list=TRUE)

kruskal.test(mean.T ~ temp,
             data = SOLlimTex.no[SOLlimTex.no$lambda == "UVB",])
dunn.test(SOLlimTex.no[SOLlimTex.no$lambda == "UVB",]$mean.T, SOLlimTex.no[SOLlimTex.no$lambda == "UVB",]$temp, method="bh", list=TRUE)

kruskal.test(mean.T ~ temp,
             data = SOLlimTex.no[SOLlimTex.no$lambda == "b21",])
dunn.test(SOLlimTex.no[SOLlimTex.no$lambda == "b21",]$mean.T, SOLlimTex.no[SOLlimTex.no$lambda == "b21",]$temp, method="bh", list=TRUE)

#BY RH
kruskal.test(mean.T ~ RH,
             data = SOLlimTex.no[SOLlimTex.no$lambda == "UVA",])
dunn.test(SOLlimTex.no[SOLlimTex.no$lambda == "UVA",]$mean.T, SOLlimTex.no[SOLlimTex.no$lambda == "UVA",]$RH, method="bh", list=TRUE)

kruskal.test(mean.T ~ RH,
             data = SOLlimTex.no[SOLlimTex.no$lambda == "UVB",])
dunn.test(SOLlimTex.no[SOLlimTex.no$lambda == "UVB",]$mean.T, SOLlimTex.no[SOLlimTex.no$lambda == "UVB",]$RH, method="bh", list=TRUE)

kruskal.test(mean.T ~ RH,
             data = SOLlimTex.no[SOLlimTex.no$lambda == "b21",])
dunn.test(SOLlimTex.no[SOLlimTex.no$lambda == "b21",]$mean.T, SOLlimTex.no[SOLlimTex.no$lambda == "b21",]$RH, method="bh", list=TRUE)


#BY experiment
kruskal.test(mean.T ~ experiment,
             data = SOLlimTex.no[SOLlimTex.no$lambda == "UVA",])
dunn.test(SOLlimTex.no[SOLlimTex.no$lambda == "UVA",]$mean.T, SOLlimTex.no[SOLlimTex.no$lambda == "UVA",]$experiment, method="bh", list=TRUE)

kruskal.test(mean.T ~ experiment,
             data = SOLlimTex.no[SOLlimTex.no$lambda == "UVB",])
dunn.test(SOLlimTex.no[SOLlimTex.no$lambda == "UVB",]$mean.T, SOLlimTex.no[SOLlimTex.no$lambda == "UVB",]$experiment, method="bh", list=TRUE)

kruskal.test(mean.T ~ experiment,
             data = SOLlimTex.no[SOLlimTex.no$lambda == "b21",])
dunn.test(SOLlimTex.no[SOLlimTex.no$lambda == "b21",]$mean.T, SOLlimTex.no[SOLlimTex.no$lambda == "b21",]$experiment, method="bh", list=TRUE)

#BY experiment just 90

# Day 3
kruskal.test(mean.T ~ experiment,
             data = SOLlimTexbase[SOLlimTexbase$lambda %in% c("UVA") & SOLlimTexbase$hour.rescale %in% c(4) & SOLlimTexbase$RH %in% c(90) & SOLlimTexbase$base %in% c("b1", "b2", "b3", "b4", "b5", "b6"),])
dunn.test(SOLlimTexbase[SOLlimTexbase$lambda %in% c("UVA") & SOLlimTexbase$hour.rescale %in% c(4) & SOLlimTexbase$RH %in% c(90),]$mean.T, 
          SOLlimTexbase[SOLlimTexbase$lambda %in% c("UVA") & SOLlimTexbase$hour.rescale %in% c(4) & SOLlimTexbase$RH %in% c(90),]$experiment, 
          method="bh", list=TRUE)
# Day 4
kruskal.test(mean.T ~ experiment,
             data = SOLlimTexbase[SOLlimTexbase$lambda %in% c("UVA") & SOLlimTexbase$hour.rescale %in% c(4) & SOLlimTexbase$RH %in% c(90) & SOLlimTexbase$base %in% c("b1", "b2", "b3", "b4", "b5", "b6"),])
dunn.test(SOLlimTexbase[SOLlimTexbase$lambda %in% c("UVA") & SOLlimTexbase$hour.rescale %in% c(4) & SOLlimTexbase$RH %in% c(90) & SOLlimTexbase$base %in% c("b1", "b2", "b3", "b4", "b5", "b6"),]$mean.T, 
          SOLlimTexbase[SOLlimTexbase$lambda %in% c("UVA") & SOLlimTexbase$hour.rescale %in% c(4) & SOLlimTexbase$RH %in% c(90) & SOLlimTexbase$base %in% c("b1", "b2", "b3", "b4", "b5", "b6"),]$experiment, 
          method="bh", list=TRUE)
# dunn.test(SOLlimTexbase[SOLlimTexbase$lambda %in% c("UVA") & SOLlimTexbase$hour.rescale %in% c(4) & SOLlimTexbase$RH %in% c(90) ,]$mean.T, 
#           SOLlimTexbase[SOLlimTexbase$lambda %in% c("UVA") & SOLlimTexbase$hour.rescale %in% c(4) & SOLlimTexbase$RH %in% c(90) ,]$experiment, 
#           method="bh", list=TRUE)



# test<-pairwise.wilcox.test(SOLlimTex.no[SOLlimTex.no$lambda == "b21",]$mean.T, SOLlimTex.no[SOLlimTex.no$lambda == "b21",]$experiment, p.adjust.method = "BH")
# test$p.value
# 
# library(corrplot)
# palette = colorRampPalette(c("green", "yellow", "red", "blue")) (40)
# corrplot(as.matrix(cha$P.adjusted), col=palette, pch.cex=2, method = "number")
# 
# 
# #HOW TO MAKE A LIST INTO A SYM MATRIX FINALLLLLLY
# raw<-as.data.frame(cbind(cha$comparisons, cha$P.adjusted))
# library(stringr)
# foo <- data.frame(do.call('rbind', strsplit(as.character(raw[,1]),'-',fixed=TRUE)))
# foo$padj<-raw[,2]
# colnames(foo)<-c("var1", "var2", "padj")
# 
# fuck <- reshape(foo, direction="wide", idvar="var2", timevar="var1")
# rownames(fuck)<-fuck[,1]
# fuck<-fuck[,-1]
# fuckmat<-as.matrix(fuck)
# fuckdist<-as.dist(fuckmat)
# 
# fuckdf<-as.data.frame(fuckmat)
# colnames(fuckdf)<-c("Exp1.1 ", "Exp10 ",  "Exp2.1 ", "Exp3 "  , "Exp4 "  , "Exp5 "  , "Exp6 "  , "Exp7 "  , "Exp8 "  )
# 
# 
# corrplot(as.matrix(fuckdist), is.corr = F, method = "square", col=palette, order="alphabet", tl.pos = "d")
# corrplot(as.matrix(fuckdist), type = "lower", order = "alphabet", tl.col = "black", 
#          tl.srt = 45, method="number", col=palette, tl.pos = "d", p.mat = as.matrix(fuckdist),
#          sig.level = 0.025, insig = "blank", cl.pos = "l")
# 
# library(ggcorrplot)
# 
# ggcorrplot(as.matrix(fuckdist),method = "square", hc.order = TRUE, ggtheme = ggplot2::theme_bw,
#            colors = c("#375e97", "#ffbb00", "#3f681c"), lab = TRUE, p.mat = as.matrix(fuckdist))
# 
# ggcorrplot(as.matrix(fuckdist), sig.level=0.025, lab_size = 4 , p.mat = NULL, 
#            insig = c("pch", "blank"), pch = 1, pch.col = "black", pch.cex =1,
#            tl.cex = 14, colors = c("yellow", "orange", "red"), lab = TRUE) +
#   theme(axis.text.x = element_text(margin=margin(-2,0,0,0)),  # Order: top, right, bottom, left
#         axis.text.y = element_text(margin=margin(0,-2,0,0))) +
#   geom_vline(xintercept=1:ncol(mtcars)-0.5, colour="white", size=2) +
#   geom_hline(yintercept=1:ncol(mtcars)-0.5, colour="white", size=2)
# 
# 
#   
# 
# ggplot(melt(as.matrix(fuckdist)), aes(Var1, Var2, fill=value)) +
#   geom_tile(height=0.8, width=0.8) +
#   scale_fill_gradient2(low="yellow", mid="orange", high="red") +
#   theme_bw() +
#   coord_equal() +
#   labs(x="",y="",fill="Adj. P-value") +
#   theme(axis.text.x=element_text(size=10, angle=90, vjust=0, hjust=0, 
#                                  margin=margin(-3,0,0,0)),
#         axis.text.y=element_text(size=10, margin=margin(0,-3,0,0)),
#         panel.grid.major=element_blank()) 
# 
# 
#   scale_x_discrete(limits=c("Exp1.1", "Exp2.1", "Exp3", "Exp4", "Exp5", "Exp6", "Exp7", 'Exp8', "Exp9", "Exp10"))
#   scale_y_discrete(limits=c("Exp1.1", "Exp2.1", "Exp3", "Exp4", "Exp5", "Exp6", "Exp7", 'Exp8', "Exp9", "Exp10"))
# 
#   level_order <-c("Exp1.1", "Exp2.1", "Exp3", "Exp4", "Exp5", "Exp6", "Exp7", 'Exp8', "Exp9", "Exp10")
# 


# table = cbind.data.frame(table$comparisons,table$Z,table$P.adjusted)
# table[order(table$table$`P.adjusted`),]
# 
# library(ggpubr)
# my_comparisons <- list(c("10", "15"), c("10", "20"), c("10", "25"), c("15", "20"), c("15", "25"), c("20", "25"))
# 
# ggboxplot(ALTlimTex.no[ALTlimTex.no$lambda == "UVA",], x = "temp", y = "mean.T",
#           fill = "RH", palette = "jco", outlier.shape = NA)+ 
#   stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
#   stat_compare_means(label.y = 1)     # Add global p-value
# 
# position = position_dodge(preserve = "single")





# ANOVA - but dont do because data is not normal
#By lambda
ALTsumLAMBDArez<-aov(AliveSpores ~ lambda + hour.rescale, data = ALTsum)
summary(glht(ALTsumLAMBDArez, linfct = mcp(lambda = "Tukey"))) #takes unbalanced design into consideration


#BY temp
ALTsumUVArez<-aov(AliveSpores ~ temp + hour.rescale, data = ALTsum[ALTsum$lambda == "UVA",])
summary(glht(ALTsumUVArez, linfct = mcp(temp = "Tukey"))) #takes unbalanced design into consideration

ALTsumUVBrez<-aov(AliveSpores ~ temp + hour.rescale, data = ALTsum[ALTsum$lambda == "UVB",])
summary(glht(ALTsumUVBrez, linfct = mcp(temp = "Tukey"))) #takes unbalanced design into consideration

ALTsumDARKrez<-aov(AliveSpores ~ temp + hour.rescale, data = ALTsum[ALTsum$lambda == "b21",])
summary(glht(ALTsumDARKrez, linfct = mcp(temp = "Tukey"))) #takes unbalanced design into consideration


#BY RH
ALTsumUVArezRH<-aov(AliveSpores ~ RH + hour.rescale, data = ALTsum[ALTsum$lambda == "UVA",])
summary(glht(ALTsumUVArezRH, linfct = mcp(RH = "Tukey"))) #takes unbalanced design into consideration

ALTsumUVBrezRH<-aov(AliveSpores ~ RH + hour.rescale, data = ALTsum[ALTsum$lambda == "UVB",])
summary(glht(ALTsumUVBrezRH, linfct = mcp(RH = "Tukey"))) #takes unbalanced design into consideration

ALTsumDARKrezRH<-aov(AliveSpores ~ RH + hour.rescale, data = ALTsum[ALTsum$lambda == "b21",])
summary(glht(ALTsumDARKrezRH, linfct = mcp(RH = "Tukey"))) #takes unbalanced design into consideration

#BY Experiment
ALTsumUVArezEX<-aov(AliveSpores ~ experiment + hour.rescale, data = ALTsum[ALTsum$lambda == "UVA",])
summary(glht(ALTsumUVArezEX, linfct = mcp(experiment = "Tukey"))) #takes unbalanced design into consideration

ALTsumUVBrezEX<-aov(AliveSpores ~ experiment + hour.rescale, data = ALTsum[ALTsum$lambda == "UVB",])
summary(glht(ALTsumUVBrezEX, linfct = mcp(experiment = "Tukey"))) #takes unbalanced design into consideration

ALTsumDARKrezEX<-aov(AliveSpores ~ experiment + hour.rescale, data = ALTsum[ALTsum$lambda == "b21",])
summary(glht(ALTsumDARKrezEX, linfct = mcp(experiment = "Tukey"))) #takes unbalanced design into consideration



plot(ALTsumDARKrez, 2)
shapiro.test(x = residuals(object = ALTsumDARKrez) )





#residuals vs. fitted
plot(ALTsumUVArez, 1)
#qq plot of residuals
plot(ALTsumUVArez, 2)

#test homogeneity of variance
library(car)
leveneTest(AliveSpores ~ temp, data = ALTsum[ALTsum$lambda == "UVA",])
#Shapiro Wilks test for normality, p>.5 means normal
shapiro.test(x = residuals(object = ALTsumUVArez) )









################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
# MANUAL VERIFICATION

library(MASS) # to access Animals data sets
library(scales) # to access break formatting functions
library(ggplot2)
library(gridExtra)

setwd("~/Dropbox/Madison/Biotron")

manalt<-read.csv("TOCROP/ALT_manual/MANUAL_ALT.csv", head=T, stringsAsFactors = F)
ggplot(manalt, aes(x=(MAN_alive/MAN_all), y=(Count...Alive.Spores/Count...All.Spores)))+geom_point()+theme_bw()

fit<-lm((Count...Alive.Spores/Count...All.Spores)~(MAN_alive/MAN_all), data=manalt)
ggplot(manalt, aes(x=(MAN_alive/MAN_dead), y=(Count...Alive.Spores/Count...Dead.Spores)))+
  geom_abline(intercept = 0, color="darkgrey", linetype=2)+
  geom_point(color="#6666ff", cex=2)+theme_bw()+
  #geom_smooth(method="lm", col="darkgrey")+
  labs(x="Manual Count - Alive:Dead", y="Automated Count - Alive:Dead", title="A. alternata",
    caption = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 1),
                     "Intercept =",signif(fit$coef[[1]],1 ),
                     " Slope =",signif(fit$coef[[2]], 1),
                     " P =",signif(summary(fit)$coef[2,4], 1)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())

mansol<-read.csv("TOCROP/SOL_manual/MANUAL_SOL.csv", head=T, stringsAsFactors = F)
ggplot(mansol, aes(x=(MAN_alive/MAN_all), y=(Count...Alive.Spores/Count...Total.Spores)))+geom_point()+theme_bw()

fit<-lm((Count...Alive.Spores/Count...All.Spores)~(MAN_alive/MAN_all), data=mansol)
ggplot(mansol, aes(x=(MAN_alive/MAN_dead), y=(Count...Alive.Spores/Count...Dead.Spores)))+
  geom_abline(intercept = 0, color="darkgrey", linetype=2)+
  geom_point(color="#ff0000", cex=2)+theme_bw()+
  #geom_smooth(method="lm", col="darkgrey")+
  labs(x="Manual Count - Alive:Dead", y="Automated Count - Alive:Dead", title="A. solani",
       caption = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 1),
                       "Intercept =",signif(fit$coef[[1]],1 ),
                       " Slope =",signif(fit$coef[[2]], 1),
                       " P =",signif(summary(fit)$coef[2,4], 1)))+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+ annotation_logticks()  +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.background = element_blank())


grid.arrange(
  
  ggplot(manalt, aes(x=(MAN_alive/MAN_dead), y=(Count...Alive.Spores/Count...Dead.Spores)))+
    geom_abline(intercept = 0, color="darkgrey", linetype=2)+
    geom_point(color="#6666ff", cex=2)+theme_bw()+
    #geom_smooth(method="lm", col="darkgrey")+
    labs(x="Manual Count - Alive:Dead", y="Automated Count - Alive:Dead", title="A. alternata",
         caption = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 1),
                         "Intercept =",signif(fit$coef[[1]],1 ),
                         " Slope =",signif(fit$coef[[2]], 1),
                         " P =",signif(summary(fit)$coef[2,4], 1)))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank()),
  
  ggplot(mansol, aes(x=(MAN_alive/MAN_dead), y=(Count...Alive.Spores/Count...Dead.Spores)))+
    geom_abline(intercept = 0, color="darkgrey", linetype=2)+
    geom_point(color="#ff0000", cex=2)+theme_bw()+
    #geom_smooth(method="lm", col="darkgrey")+
    labs(x="Manual Count - Alive:Dead", y="Automated Count - Alive:Dead", title="A. solani",
         caption = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 1),
                         "Intercept =",signif(fit$coef[[1]],1 ),
                         " Slope =",signif(fit$coef[[2]], 1),
                         " P =",signif(summary(fit)$coef[2,4], 1)))+
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)))+ annotation_logticks()  +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank()),
  
  nrow=1, ncol=2
)






###################################
# Ask TEDWARD, but I think the fact that they have different start prop alive doesnt matter, because the modeling is comparing slopes
# 
# ALTlimmury <- ALTlim %>%
#   group_by(hour.rescale, lambda) %>%
#   summarise( 
#     n=n(),
#     mean.pr=mean(AliveSpores/TotalSpores),
#     sd.pr=sd(AliveSpores/TotalSpores)
#   ) %>%
#   mutate (se.pr=sd.pr/sqrt(n)) 
# 
# ALTsum2 <- ALTlim %>%
#   group_by(experiment, RH, temp, lambda, flux, BLOCK, hour.rescale, cart, base) %>%
#   summarize(TotalSpores = mean(TotalSpores),
#             AliveSpores = mean(AliveSpores))
# 
# ALTsum2.norm <- ALTsum2 %>%
#   group_by(experiment, RH, temp, lambda, flux, BLOCK, hour.rescale, cart, base) %>%
#   mutate(TotalSpores.norm=TotalSpores/mean(ALTsum2$TotalSpores[ALTsum2$hour.rescale==0])) %>%
#   mutate(AliveSpores.norm=AliveSpores/mean(ALTsum2$AliveSpores[ALTsum2$hour.rescale==0]))
# 
# 
# alt.mod2.norm<-glmmTMB(AliveSpores.norm~RH+hour.rescale+lambda+(1|experiment/cart:base)+offset(log(TotalSpores.norm)),
#                   data=ALTsum2.norm, family=poisson)
# 
# hist(log(ALTsum2.norm$TotalSpores.norm))
# 
# alt.mod2no0.norm<-glmmTMB(AliveSpores.norm~RH+hour.rescale+lambda+(1|experiment/cart:base)+offset(log(TotalSpores.norm)),
#                      data=ALTsum2.norm[ALTsum2.norm$hour.rescale != 0,], family=poisson)
# 
# plot_model(alt.mod2.norm, show.values = TRUE, value.offset = .2, dot.size = 3, transform="exp", vline.color = "black")
# plot_model(alt.mod2no0.norm, show.values = TRUE, value.offset = .2, dot.size = 3, transform="exp", vline.color = "black")
# 
# 
# ggplot(ALTsum2.norm[ALTsum2.norm$lambda == "UVA",], aes(x=temp, y=AliveSpores.norm/TotalSpores.norm, group=temp))+geom_boxplot()
# 
# 
# ALTsum2.norm.sum <- ALTsum2.norm %>%
#   group_by(hour.rescale, lambda, experiment) %>%
#   summarise( 
#     n=n(),
#     mean=mean(AliveSpores.norm/TotalSpores.norm),
#     sd=sd(AliveSpores.norm/TotalSpores.norm)
#   ) %>%
#   mutate (se=sd/sqrt(n)) 
# 
# 
# ggplot(ALTsum2.norm.sum[ALTsum2.norm.sum$lambda == "UVA",], aes(x=hour.rescale, y=mean, color=experiment))+geom_line()



###################### Tedward trying to figure outpiecewise regression

# SOLsum90$period <- SOLsum90$hour.rescale <=1
# 
# sol.T90mod2<-glmmTMB(AliveSpores~temp+hour.rescale+lambda+(1|experiment/cart:base)+offset(log(TotalSpores)),
#                      data=SOLsum90, family=nbinom1)
# 
# 
# sol.T90mod2_piece<-glmmTMB(AliveSpores~temp+ hour.rescale*period +lambda+(1|experiment/cart:base)+offset(log(TotalSpores)),
#                            data=SOLsum90, family=nbinom1)
# 
# 
# predframe = data.frame(lambda = "UVA",  temp = 10, period = c(T,T,F,F,F), hour.rescale = seq(0,4,1), TotalSpores = 250, experiment = "Exp10", cart="Cart1", base="b2")
# 
# p1 <- predict(sol.T90mod2, predframe)
# p2 <- predict(sol.T90mod2_piece,predframe)
# 
# plot(seq(0,4,1),exp(p1))
# plot(seq(0,4,1),exp(p2))



