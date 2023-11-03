setwd("~/Documents/2020_AphalloidesPopgen/")
library(adegenet)
library(dplyr)
library(poppr)
library(ggplot2)
library(sp)

##AFLP

coords_aflp<-read.csv("data/ALL_COORDS_geo.csv", head=T, stringsAsFactors = F)


hier<-read.genalex("data/aflp/AFLP_genalex_AMANITABASED.csv")

pops.aflp<-c(
  replicate(13, "roch_1"),
  replicate(7, "roch_2"),
  replicate(2, "roch_3"),
  replicate(1, "Austria"),
  replicate(37, "CESAC06"),
  replicate(25, "CESAC02"),
  replicate(19, "Drake1"),
  replicate(13, "Drake2"),
  replicate(8, "Drake3"),
  replicate(6, "Drake4"),
  replicate(1, "Spain_2006"),
  replicate(7, "HD1"),
  replicate(7, "HD2"),
  replicate(4, "HD3"),
  replicate(44, "JL"),
  replicate(1, "Arfons_2007"),
  replicate(1, "Monterrey_Pisto_2006"),
  replicate(14, "S"),
  replicate(11, "str")
)
strata(hier)$Pop<-pops.aflp

# fix pops
coords_aflp$pop<-ifelse(coords_aflp$pop == "roch",
                        paste("roch",gsub("\\_.*","",coords_aflp$pop.t1),sep="_"),
                        coords_aflp$pop)

# fix coordinates for Serbia
coords_aflp$x<-ifelse(coords_aflp$pop=="S",coords_aflp$x/100,coords_aflp$x)
coords_aflp$y<-ifelse(coords_aflp$pop=="S",coords_aflp$y/100,coords_aflp$y)

#get genets
hier.gc<-as.genclone(hier)
mlg<-data.frame(ID=indNames(hier.gc),mlg=hier.gc@mlg[])

### calculate distances

coords_aflp$mlg<-mlg$mlg[match(coords_aflp$FINAlab, mlg$ID,)]

coords_aflp<-coords_aflp[coords_aflp$FINAlab%in%mlg$ID,]

coords_aflpchull<-coords_aflp %>% 
  group_by(pop,mlg) %>%
  summarise(meanX=mean(x),
            meanY=mean(y),
            mindiameter=ifelse(n()>1,max(dist(matrix(c(x,y),ncol=2))),NA),
            chull=ifelse(n()>2,Polygon(matrix(c(x,y),ncol=2)[chull(x,y),])@area,NA))%>%
  as.data.frame()


popnametranslater<-data.frame(
  oldname=c("S","CESAC02","CESAC06",
            "JL","str","roch_1","roch_2","roch_3",
            "Drake1","Drake2","Drake3","Drake4",
            "HD1","HD2","HD3"),
  goodname=c("Serbia","CESAC 2002","CESAC 2006","Jake's Landing","Round Valley",
             "Rochester 1","Rochester 2","Rochester 3",
             "Drake 1","Drake 2","Drake 3","Drake 4",
             "Heart's Desire 1","Heart's Desire 2","Heart's Desire 3"),
  color=c("#0E575B","#1CADB6","#A4DEE2","#001399","#99A1D6",
          "#fc9105","#cc8f02","#948103",
          "#D61A3C","#77001D","gray","#943541",
          "#E26313","#D55209","#E27429"))

coords_aflpchull$pop<-popnametranslater$goodname[match(coords_aflpchull$pop,popnametranslater$oldname)]
cols<-popnametranslater$color
names(cols)<-popnametranslater$goodname
  
ggplot(coords_aflpchull)+
  geom_boxplot(aes(x=pop,y=mindiameter))+theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(coords_aflpchull)+
  geom_boxplot(aes(x=pop,y=chull))+theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

dir.create("results/chull")
write.csv(coords_aflpchull,file = "results/chull/chullaflp.csv",row.names = F)

##SNP

coords_snp<-read.csv("data/GENETS_coords_THISONE.csv", head=T, stringsAsFactors = F) #just so that I had the right labels
mlg<-read.csv("data/Golan_mlg_loose.csv")
colnames(mlg)[2]<-"mlg"
coords_snp<-merge(coords_snp,mlg,by.y="sample",by.x="Sample_Name")
coords_snp$mlg<-as.character(coords_snp$mlg)


### calculate distances

coords_snpchull<-coords_snp %>% 
  group_by(pop,mlg) %>%
  summarise(meanX=mean(x),
            meanY=mean(y),
            mindiameter=ifelse(n()>1&sum(!is.na(x))>1,max(dist(matrix(c(x,y),ncol=2))),NA),
            chull=ifelse(n()>2&sum(!is.na(x))>2,Polygon(matrix(c(x,y),ncol=2)[chull(x,y),])@area,NA))%>%
  as.data.frame()


ggplot(coords_snpchull)+
  geom_boxplot(aes(x=pop,y=mindiameter))+theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(coords_snpchull)+
  geom_boxplot(aes(x=pop,y=chull))+theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


write.csv(coords_snpchull,file = "results/chull/chullsnp.csv",row.names = F)

coords_aflpchull$method<-"AFLP"
coords_snpchull$method<-"SNP"

chullmerge<-rbind(coords_aflpchull,coords_snpchull)

ggplot(chullmerge)+
  geom_boxplot(aes(x=pop,y=mindiameter,color=method))+theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggplot(chullmerge)+
  geom_boxplot(aes(x=pop,y=chull,color=method))+theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



wantedpopslevels<-c("Mira","CESAC 2002","CESAC 2006","Serbia",
                    "Round Valley","Rochester 1","Rochester 2","Rochester 3","Jake's Landing",
                    "Heart's Desire 1","Heart's Desire 2","Heart's Desire 3",
                    "Drake 1","Drake 2","D2.04","D2.14","D2.15",
                    "Drake 3","D3.04","D3.14","D3.15","Drake 4")

chullmergeforfinalplot<-chullmerge[chullmerge$pop%in%wantedpopslevels,]

chullmergeforfinalplot$pop<-factor(chullmergeforfinalplot$pop,levels=wantedpopslevels)

pdf("Figs/genetsize.pdf")

tempforplot<-chullmergeforfinalplot%>%group_by(pop)%>%
  summarise(maxmindiam=max(mindiameter,na.rm=T),
            maxchull=ifelse(all(is.na(chull)),NA,round(max(chull,na.rm = T),digits = 2)))

ggplot()+
  geom_boxplot(data=chullmergeforfinalplot,aes(x=pop,y=mindiameter,color=method))+
  theme_bw()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_text(data=tempforplot,aes(x=pop,y=maxmindiam,label=maxchull),
               hjust=0.5,vjust = -0.9)+
  lims(y=c(0,12))
dev.off()


ggplot(chullmergeforfinalplot)+
  geom_boxplot(aes(x=pop,y=chull,color=method))+theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



pdf("Figs/genetsize.pdf")

outhull<-chullmerge%>%group_by(pop)%>%
  summarise(maxmindiam=max(mindiameter,na.rm=T),
            maxchull=ifelse(all(is.na(chull)),NA,round(max(chull,na.rm = T),digits = 2)))
write.csv(outhull,"results/chull/maxchullperpop.csv")
