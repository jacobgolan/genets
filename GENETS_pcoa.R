################################ ################################ ################################ ################################ 
################################ ################################ ################################ ################################ 
################################ ################################ ################################ ################################ 
################################ ################################ ################################ ################################ 
################################ ################################ ################################ ################################ 
#                                                               PCoA
################################ ################################ ################################ ################################ 
################################ ################################ ################################ ################################ 
################################ ################################ ################################ ################################ 
################################ ################################ ################################ ################################ 
################################ ################################ ################################ ################################ 



setwd("~/Dropbox/finalizing_genets")
library(adegenet)
library(poppr)
library(ggplot2)

library(devtools)
#BiocManager::install(c("SNPRelate", "qvalue"))
#install_github("green-striped-gecko/dartR")
library(dartR)


maf006<-read.PLINK("HESS_fix/mx0.5/vcf.MMDP60.maf0.006.minQ30.recode.raw", 
                                       map.file = "HESS_fix/mx0.5/vcf.MMDP60.maf0.006.minQ30.recode.map", quiet = FALSE,
                                       parallel = require("parallel"), n.cores = NULL)

PCoA.snp<-ape::pcoa(as.matrix(bitwise.dist(maf006, percent = F, missing_match = T, scale_missing = T)), correction="none", rn=NULL)
#biplot(PCoA.snp, cex=.02)
1.695416e-01*100
9.748778e-02*100

PCoA.snp2<-PCoA.snp
PCoA.snp2$vectors<-apply(PCoA.snp2$vectors,2,scale,center=T,scale=T)
rownames(PCoA.snp2$vectors)<-rownames(PCoA.snp$vectors)
victors<-PCoA.snp2$vectors
victors<- data.frame(names = row.names(victors), victors) #convert rown names to column
setwd("~/Dropbox/GENETS_finishing")
coords<-read.csv("GENETS_coords_THISONE.csv", head=T, stringsAsFactors = F) #just so that I had the right labels
PCoA.axis.labelled<-merge(victors, coords, by.x="names", by.y="Sample_Name")

forms<-c('turd',	'turd',	'turd',	'turd',	'tri',	'turd',	'sta',	'turd',	'cir',	'cir',	'cir',	'cir',	'cir',	'cir',	'cir',	
         'cir',	'cir',	'cir',	'cir',	'cir',	'cir',	'cir',	'cir',	'cir',	'cir',	'cir',	'tri',	'squ',	'squ',	'squ',
         'squ',	'squ',	'squ',	'squ',	'squ',	'squ',	'squ',	'squ',	'squ',	'squ',	'squ',	'squ',	'squ',	'squ',	'squ',
         'squ',	'squ',	'squ',	'squ',	'squ',	'squ',	'squ',	'squ',	'squ',	'squ',	'squ',	'squ',	'squ',	'squ',	'squ',
         'squ',	'tri',	'tri',	'tri',	'tri',	'tri',	'tri',	'tri',	'tri',	'tri',	'tri',	'tri',	'rho',	'rho',	'rho',
         'rho',	'rho',	'rho',	'rho',	'rho',	'rho',	'rho',	'rho',	'rho',	'rho',	'rho')

PCoA.axis.labelled$forms<-forms

ggplot(PCoA.axis.labelled, aes(x=Axis.1, y=Axis.2, size=3, color=factor(PCoA.axis.labelled$pop), shape=as.factor(PCoA.axis.labelled$forms)))+geom_point()+
  theme_classic(base_size = 16)+labs(x="PCoA Axis 1", y="PCoA Axis 2", title="PCoA: MAF0.006, Euclidean Distance (genetic)", fill="Population")+ 
  theme(legend.title=element_blank())+
  scale_shape_manual(values = c("tri"=17, "cir"=19, "squ"=15, "rho"=18, "sta"=8, "turd"= 4))+
  scale_color_manual(values=c("#cc33ff", "#663300", "#ff0000", "#b30000", "#ff6666", "#ff9900", "#cc7a00", "#cc9900", "#9966ff", "#660066", "#6600ff"))+
  guides(colour = guide_legend(override.aes = list(size=4)),shape = guide_legend(override.aes = list(size=4)))
#scale_color_discrete(breaks=c("CA_old","D2.04","D2.14","D2.15","D3.04","D3.14","D3.15","Agraria","Mira","Vilarinho", "old_EU"))

# No legend
ggplot(PCoA.axis.labelled, aes(x=Axis.1, y=Axis.2, size=3, color=factor(PCoA.axis.labelled$pop), shape=as.factor(PCoA.axis.labelled$forms)))+geom_point()+
  theme_classic(base_size = 16)+labs(x="PCoA Axis 1", y="PCoA Axis 2", title="PCoA: MAF0.006, Euclidean Distance (genetic)", fill="Population")+ 
  theme(legend.position="none")+
  scale_shape_manual(values = c("tri"=17, "cir"=19, "squ"=15, "rho"=18, "sta"=8, "turd"= 4))+
  scale_color_manual(values=c("#cc33ff", "#663300", "#ff0000", "#b30000", "#ff6666", "#ff9900", "#cc7a00", "#cc9900", "#9966ff", "#660066", "#6600ff"))+
  guides(colour = guide_legend(override.aes = list(size=4)),shape = guide_legend(override.aes = list(size=4)))










####### Ellipses
PCoA.axis.labelled$region<-PCoA.axis.labelled$name

library(plyr)
PCoA.axis.labelled$region<-revalue(PCoA.axis.labelled$region, c(
  '10004'="EU",	'10007'="EU",	'10016'="EU",	'10018'="EU",	'10019'="EU",	'10169'="EU",	'10170'="D2",	'10171'="EU",	'10221'="D2",	'10222'="D2",	'10223'="D2",	
  '10224'="D2",	'10225'="D2",	'10226'="D2",	'10227'="D2",	'10228'="D2",	'10229'="D2",	'10230'="D2",	'10231'="D2",	'10232'="D2",	'10233'="D2",	'10237'="D3",	
  '10238'="D3",	'10239'="D3",	'10240'="D3",	'10241'="D3",	'10277'="EU",	'10280'="D2",	'10281'="D2",	'10282'="D2",	'10283'="D2",	'10287'="D2",	'10288'="D2",	
  '10292'="D2",	'10293'="D2",	'10294'="D2",	'10295'="D2",	'10298'="D2",	'10299'="D2",	'10300'="D2",	'10301'="D2",	'10303'="D2",	'10304'="D2",	'10306'="D2",	
  '10309'="D2",	'10326'="D2",	'10327'="D2",	'10328'="D2",	'10329'="D2",	'10330'="D2",	'10331'="D2",	'10334'="D2",	'10347'="D3",	'10348'="D3",	'10349'="D3",	
  '10350'="D3",	'10354'="D3",	'10355'="D3",	'10356'="D3",	'10380'="D3",	'10384'="D3",	'10502'="EU",	'10503'="EU",	'10504'="EU",	'10505'="EU",	'10506'="EU",	
  '10508'="EU",	'10509'="EU",	'10510'="EU",	'10511'="EU",	'10512'="EU",	'10513'="EU",	'10707'="D2",	'10708'="D2",	'10709'="D2",	'10710'="D2",	'10711'="D2",	
  '10712'="D2",	'10713'="D2",	'10715'="D2",	'10716'="D2",	'10717'="D2",	'10718'="D2",	'10719'="D3",	'10720'="D3",	'10721'="D3"
))



ggplot(data=NULL)+
  geom_point(data=PCoA.axis.labelled, aes(x=Axis.1, y=Axis.2, cex=5,color=factor(PCoA.axis.labelled$pop), shape=as.factor(PCoA.axis.labelled$forms)))+
  stat_ellipse(data=PCoA.axis.labelled, aes(x=Axis.1, y=Axis.2, group=PCoA.axis.labelled$region), linetype="dotted")+
  scale_shape_manual(values = c("tri"=17, "cir"=19, "squ"=15, "rho"=18, "sta"=8, "turd"= 4))+
  scale_color_manual(values=c("#cc33ff", "#663300", "#ff0000", "#b30000", "#ff6666", "#ff9900", "#cc7a00", "#cc9900", "#9966ff", "#660066", "#6600ff"))+
  labs(x="PCoA Axis 1", y="PCoA Axis 2", title="", fill="")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1, fill="white"), legend.position="none")

# tres amigos (3 ellipses)
PCoA.axis.labelled2<-PCoA.axis.labelled
PCoA.axis.labelled2$region<-PCoA.axis.labelled2$name

PCoA.axis.labelled2$region<-revalue(PCoA.axis.labelled2$region, c(
  '10004'="EU",	'10007'="EU",	'10016'="EU",	'10018'="EU",	'10019'="EU",	'10169'="EU",	'10170'="D2",	'10171'="EU",	'10221'="D2",	'10222'="D2",	'10223'="D2",	
  '10224'="D2",	'10225'="D2",	'10226'="D2",'10227'="D2",	'10228'="D2",	'10229'="D2",	'10230'="D2",	'10231'="D2",	'10232'="D2",	'10233'="D2",	'10237'="D3",	
  '10238'="D3",	'10239'="D3",	'10240'="D3",	'10241'="D3",	'10277'="EU",	'10280'="D2s",'10281'="D2s",'10282'="D2s",'10283'="D2s",'10287'="D2",	'10288'="D2",	
  '10292'="D2",	'10293'="D2",	'10294'="D2",	'10295'="D2",	'10298'="D2",	'10299'="D2",	'10300'="D2",	'10301'="D2",	'10303'="D2",	'10304'="D2s",'10306'="D2s",	
  '10309'="D2",	'10326'="D2s",'10327'="D2s",'10328'="D2s",'10329'="D2s",'10330'="D2s",'10331'="D2s",'10334'="D2",'10347'="D3",'10348'="D3",'10349'="D3",	
  '10350'="D3",	'10354'="D3",	'10355'="D3",	'10356'="D3",	'10380'="D3",	'10384'="D3",	'10502'="EU",	'10503'="EU",	'10504'="EU",	'10505'="EU",	'10506'="EU",	
  '10508'="EU",	'10509'="EU",	'10510'="EU",	'10511'="EU",	'10512'="EU",	'10513'="EU",	'10707'="D2",	'10708'="D2",	'10709'="D2",	'10710'="D2s",'10711'="D2s",	
  '10712'="D2s",'10713'="D2",	'10715'="D2s",'10716'="D2",	'10717'="D2",	'10718'="D2",	'10719'="D3",	'10720'="D3",	'10721'="D3"
))



ggplot(data=NULL)+
  geom_point(data=PCoA.axis.labelled2, aes(x=Axis.1, y=Axis.2, cex=5,color=factor(PCoA.axis.labelled2$pop), shape=as.factor(PCoA.axis.labelled2$forms)))+
  stat_ellipse(data=PCoA.axis.labelled2, aes(x=Axis.1, y=Axis.2, group=PCoA.axis.labelled2$region), linetype="dotted")+
  scale_shape_manual(values = c("tri"=17, "cir"=19, "squ"=15, "rho"=18, "sta"=8, "turd"= 4))+
  scale_color_manual(values=c("#cc33ff", "#663300", "#ff0000", "#b30000", "#ff6666", "#ff9900", "#cc7a00", "#cc9900", "#9966ff", "#660066", "#6600ff"))+
  labs(x="PCoA Axis 1", y="PCoA Axis 2", title="", fill="")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1, fill="white"), legend.position="none")

#AXIS 3
# ggplot(PCoA.axis.labelled, aes(x=Axis.1, y=Axis.3, size=3, color=factor(PCoA.axis.labelled$pop), shape=as.factor(PCoA.axis.labelled$forms)))+geom_point()+
#   theme_classic(base_size = 16)+labs(x="PCoA Axis 1", y="PCoA Axis 3", title="PCoA: MAF0.006, Euclidean Distance (genetic)", fill="Population")+ 
#   theme(legend.title=element_blank())+
#   scale_shape_manual(values = c("tri"=17, "cir"=19, "squ"=15, "rho"=18, "sta"=8))+
#   scale_color_manual(values=c("#cc33ff", "#663300", "#ff0000", "#b30000", "#ff6666", "#ff9900", "#cc7a00", "#cc9900", "#9966ff", "#660066", "#6600ff"))+
#   guides(colour = guide_legend(override.aes = list(size=4)),shape = guide_legend(override.aes = list(size=4)))
# #scale_color_discrete(breaks=c("CA_old","D2.04","D2.14","D2.15","D3.04","D3.14","D3.15","Agraria","Mira","Vilarinho", "old_EU"))

# points labelled
# ggplot(PCoA.axis.labelled, aes(x=Axis.1, y=Axis.2, size=3, color=factor(PCoA.axis.labelled$pop), shape=as.factor(PCoA.axis.labelled$forms)))+geom_point()+
#   theme_classic(base_size = 16)+labs(x="PCoA Axis 1", y="PCoA Axis 2", title="PCoA: MAF0.006, Euclidean Distance (genetic)", fill="Population")+ 
#   theme(legend.title=element_blank())+
#   scale_shape_manual(values = c("tri"=17, "cir"=19, "squ"=15, "rho"=18, "sta"=8))+
#   scale_color_manual(values=c("#cc33ff", "#663300", "#ff0000", "#b30000", "#ff6666", "#ff9900", "#cc7a00", "#cc9900", "#9966ff", "#660066", "#6600ff"))+
#   guides(colour = guide_legend(override.aes = list(size=4)),shape = guide_legend(override.aes = list(size=4)))+
#   geom_text_repel(aes(label = PCoA.axis.labelled$names))


#line to centroid
PCoA.axis.labelled3<-PCoA.axis.labelled
PCoA.axis.labelled3$region<-PCoA.axis.labelled3$name



PCoA.axis.labelled3$region<-revalue(PCoA.axis.labelled3$region, c(
  '10004'="EU",	'10007'="EU",	'10016'="EU",	'10018'="EU",	'10019'="EU",	'10169'="EU",	'10170'="TRH",	'10171'="EU",	'10221'="D2.04",	'10222'="D2.04",	'10223'="D2.04",	
  '10224'="D2.04",	'10225'="D2.04",	'10226'="D2s.04",'10227'="D2.04",	'10228'="D2.04",	'10229'="D2.04",	'10230'="D2.04",	'10231'="D2.04",	'10232'="D2.04",	'10233'="D2.04",	'10237'="D3",	
  '10238'="D3",	'10239'="D3",	'10240'="D3",	'10241'="D3",	'10277'="EU",	'10280'="D2s.14",'10281'="D2s.14",'10282'="D2s.14",'10283'="D2s.14",'10287'="D2.14",	'10288'="D2.14",	
  '10292'="D2.14",	'10293'="D2.14",	'10294'="D2.14",	'10295'="D2.14",	'10298'="D2.14",	'10299'="D2.14",	'10300'="D2.14",	'10301'="D2.14",	'10303'="D2.14",	'10304'="D2s.14",'10306'="D2s.14",	
  '10309'="D2.14",	'10326'="D2s.14",'10327'="D2s.14",'10328'="D2s.14",'10329'="D2s.14",'10330'="D2s.14",'10331'="D2s.14",'10334'="D2",'10347'="D3",'10348'="D3",'10349'="D3",	
  '10350'="D3",	'10354'="D3",	'10355'="D3",	'10356'="D3",	'10380'="D3",	'10384'="D3",	'10502'="EU",	'10503'="EU",	'10504'="EU",	'10505'="EU",	'10506'="EU",	
  '10508'="EU",	'10509'="EU",	'10510'="EU",	'10511'="EU",	'10512'="EU",	'10513'="EU",	'10707'="D2.15",	'10708'="D2.15",	'10709'="D2.15",	'10710'="D2s.15",'10711'="D2s.15",	
  '10712'="D2s.15",'10713'="D2.15",	'10715'="D2s.15",'10716'="D2.15",	'10717'="D2.15",	'10718'="D2.15",	'10719'="D3",	'10720'="D3",	'10721'="D3"
))

gg <- merge(PCoA.axis.labelled3,aggregate(cbind(mean.x=PCoA.axis.labelled3$Axis.1,mean.y=PCoA.axis.labelled3$Axis.2)~region,PCoA.axis.labelled3,mean),by="region")


ggplot(gg, aes(Axis.1,Axis.2,color=factor(region)))+geom_point(size=3)+
  geom_point(aes(x=mean.x,y=mean.y),size=5)+
  geom_segment(aes(x=mean.x, y=mean.y, xend=Axis.1, yend=Axis.2))


ggplot(gg, aes(Axis.1,Axis.2,shape=factor(forms),color=factor(pop)))+geom_point(size=3)+
  geom_point(aes(x=mean.x,y=mean.y),size=0)+
  geom_segment(aes(x=mean.x, y=mean.y, xend=Axis.1, yend=Axis.2), size=.5, alpha=.2)+
  scale_shape_manual(values = c("tri"=17, "cir"=19, "squ"=15, "rho"=18, "sta"=8, "turd"= 4))+
  scale_color_manual(values=c("#cc33ff", "#663300", "#ff0000", "#b30000", "#ff6666", "#ff9900", "#cc7a00", "#cc9900", "#9966ff", "#660066", "#6600ff"))+
  labs(x="PCoA Axis 1", y="PCoA Axis 2", title="", fill="")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1, fill="white"), legend.position="none")


# _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ 
#AFLP

setwd("~/Dropbox/Genet Size/Jacob Analyses")
hier <- read.genalex("AFLP_genalex_AMANITABASED.csv")
# mean(bitwise.dist(hier, percent = F, missing_match = T, scale_missing = T))
# mean(diss.dist(hier, percent = F))

PCoA.aflp<-ape::pcoa(as.matrix(diss.dist(hier, percent = F)))

#biplot(PCoA.aflp)
2.742490e-01*100
1.627170e-01*100

PCoA.aflp2<-PCoA.aflp
PCoA.aflp2$vectors<-apply(PCoA.aflp2$vectors,2,scale,center=T,scale=T)
rownames(PCoA.aflp2$vectors)<-rownames(PCoA.aflp$vectors)
victors<-PCoA.aflp2$vectors
victors<- data.frame(names = row.names(victors), victors) #convert rown names to column

#victors.aflp<-PCoA.aflp$vectors
#rownames(victors)[86]<-"10004"

#victors.aflp<- data.frame(names = row.names(victors.aflp), victors.aflp) #convert rown names to column

pops.aflp<-c(
  replicate(13, "Rochester_NY_2007_1"),
  replicate(7, "Rochester_NY_2007_2"),
  replicate(2, "Rochester_NY_2007_3"),
  replicate(1, "Austria"),
  replicate(37, "CESAC_2006"),
  replicate(25, "CESAC_2002"),
  replicate(19, "Drake_1_2004"),
  replicate(13, "Drake_2_2004"),
  replicate(8, "Drake_3_2004"),
  replicate(6, "Drake_4_2004"),
  replicate(1, "Spain_2006"),
  replicate(7, "Hearts_Desire_1_2004"),
  replicate(7, "Hearts_Desire_2_2004"),
  replicate(4, "Hearts_Desire_3_2004"),
  replicate(44, "Jakes_Landing_2006"),
  replicate(1, "Arfons_2007"),
  replicate(1, "Monterrey_Pisto_2006"),
  replicate(14, "Serbia_2007"),
  replicate(11, "New_Jersey_str_2006")
)

PCoA.axis.labelled.aflp<-cbind(victors, pops.aflp)

library(ggrepel)

# ggplot(PCoA.axis.labelled.aflp, aes(x=Axis.1, y=Axis.2))+geom_point()+
#   geom_label_repel(
#     aes(Axis.1, Axis.2, fill = factor(pops.aflp), label = names),
#     fontface = 'bold', color = 'white',
#     box.padding = 0.035, point.padding = 0.5, size = 1.25, segment.alpha = .6,
#     segment.color = 'grey50'
#   ) +
#   theme_classic(base_size = 16)+labs(x="PCoA Axis 1", y="PCoA Axis 2", title="PCoA: AFLP, Euclidean Distance (genetic)", fill="Population")

# ggplot(PCoA.axis.labelled.aflp, aes(x=Axis.1, y=Axis.2, color=factor(PCoA.axis.labelled$pops.aflp)))+geom_point()+
#   geom_text_repel(aes(label = PCoA.axis.labelled$names), segment.size = .25, size=3,segment.alpha = .4)+
#   theme_classic(base_size = 16)+labs(x="PCoA Axis 1", y="PCoA Axis 2", title="PCoA: AFLP, Euclidean Distance (genetic)", fill="Population")+ theme(legend.title=element_blank())


ggplot(PCoA.axis.labelled.aflp, aes(x=Axis.1, y=Axis.2, size=3, color=factor(PCoA.axis.labelled.aflp$pops.aflp)))+geom_point()+
  theme_classic(base_size = 16)+labs(x="PCoA Axis 1", y="PCoA Axis 2", title="PCoA: AFLP, Euclidean Distance (genetic)", fill="Population")+ 
  theme(legend.title=element_blank())+
  scale_shape_manual(values = 19)+
  scale_color_manual(values=c("#6600cc", "#d966ff", "#9900cc", "#993366", "#ff6699", "#ff0000", "#ff9900", "#cc5200", "#ffff00", "#999900","#e6b800","#0000ff",
                              "#ffcccc", "#8080ff", "#0099ff", "#006bb3", "#b3e0ff", "#df9fbf", "#4d0066"))+
  guides(colour = guide_legend(override.aes = list(size=4)))

# no legend
ggplot(PCoA.axis.labelled.aflp, aes(x=Axis.1, y=Axis.2, size=3, color=factor(PCoA.axis.labelled.aflp$pops.aflp)))+geom_point()+
  theme_classic(base_size = 16)+labs(x="PCoA Axis 1", y="PCoA Axis 2", title="PCoA: AFLP, Euclidean Distance (genetic)", fill="Population")+ 
  theme(legend.position="none")+
  scale_shape_manual(values = 19)+
  scale_color_manual(values=c("#6600cc", "#d966ff", "#9900cc", "#993366", "#ff6699", "#ff0000", "#ff9900", "#cc5200", "#ffff00", "#999900","#e6b800","#0000ff",
                              "#ffcccc", "#8080ff", "#0099ff", "#006bb3", "#b3e0ff", "#df9fbf", "#4d0066"))+
  guides(colour = guide_legend(override.aes = list(size=4)))



###################################
#with ellipse
PCoA.axis.labelled.aflp$region<-PCoA.axis.labelled.aflp$pops.aflp
# PCoA.axis.labelled.aflp$region[PCoA.axis.labelled.aflp$region %in% c("Arfons_2007","CESAC_2002",
#                                                                      "Austria", "CESAC_2006", "Spain_2006", "Serbia_2007")]<-"EU"

library(plyr)
PCoA.axis.labelled.aflp$region<-revalue(PCoA.axis.labelled.aflp$region, c("Rochester_NY_2007_1"  = "EC",
                                                                          "Rochester_NY_2007_2" =  "EC",
                                                                          "Rochester_NY_2007_3" = "EC",
                                                                          "Austria" = "EU",            
                                                                          "CESAC_2006" = "EU",           
                                                                          "CESAC_2002" = "EU",          
                                                                          "Drake_1_2004" = "CA",        
                                                                          "Drake_2_2004" = "CA",         
                                                                          "Drake_3_2004"  = "CA",       
                                                                          "Drake_4_2004"  = "CA",      
                                                                          "Spain_2006"  = "EU",           
                                                                          "Hearts_Desire_1_2004" = "CA",
                                                                          "Hearts_Desire_2_2004" = "CA", 
                                                                          "Hearts_Desire_3_2004" = "CA",
                                                                          "Jakes_Landing_2006"  = "EC", 
                                                                          "Arfons_2007" = "EU",          
                                                                          "Monterrey_Pisto_2006" = "CA", 
                                                                          "Serbia_2007" = "EU",          
                                                                          "New_Jersey_str_2006" = "EC" ))

ggplot(data=NULL)+
  geom_point(data=PCoA.axis.labelled.aflp, aes(x=Axis.1, y=Axis.2, size=3, color=factor(PCoA.axis.labelled.aflp$pops.aflp)))+
  stat_ellipse(data=PCoA.axis.labelled.aflp, aes(x=Axis.1, y=Axis.2, shape = PCoA.axis.labelled.aflp$region), linetype="dotted", cex=.5, inherit.aes = F)+
  scale_color_manual(values=c("#6600cc", "#d966ff", "#9900cc", "#993366", "#ff6699", "#ff0000", "#ff9900", "#cc5200", "#ffff00", "#999900","#e6b800","#0000ff",
                              "#ffcccc", "#8080ff", "#0099ff", "#006bb3", "#b3e0ff", "#df9fbf", "#4d0066"))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1, fill="white"))+
  labs(x="PCoA Axis 1", y="PCoA Axis 2", title="", fill="Population")+theme(legend.position="none")


# Axis 3
# ggplot(PCoA.axis.labelled.aflp, aes(x=Axis.1, y=Axis.3, size=3, color=factor(PCoA.axis.labelled.aflp$pops.aflp)))+geom_point()+
#   theme_classic(base_size = 16)+labs(x="PCoA Axis 1", y="PCoA Axis 3", title="PCoA: AFLP, Euclidean Distance (genetic)", fill="Population")+ 
#   theme(legend.title=element_blank())+
#   scale_shape_manual(values = 19)+
#   scale_color_manual(values=c("#6600cc", "#d966ff", "#9900cc", "#993366", "#ff6699", "#ff0000", "#ff9900", "#cc5200", "#ffff00", "#999900","#e6b800","#0000ff",
#                               "#ffcccc", "#8080ff", "#0099ff", "#006bb3", "#b3e0ff", "#df9fbf", "#4d0066"))+
#   guides(colour = guide_legend(override.aes = list(size=4)))