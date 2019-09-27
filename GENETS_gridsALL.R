#####################################################################################################
#####################################################################################################
#####################################################################################################
#####################################################################################################
#############################           Pop Grids           #########################################
#####################################################################################################
#####################################################################################################
#####################################################################################################
#####################################################################################################
#####################################################################################################

setwd("~/Dropbox/finalizing_genets")
library(phangorn)
library(hierfstat)
library(poppr)
library(adegenet)
library(RColorBrewer)
library(gplots)
library(vcfR)
library(magrittr)
library(poppr)
library(ggtree)



d3.15<-as.genind(as.matrix(read.PLINK("HESS_fix/pops/vcf.MMDP60.maf0.006.minQ30.recode.D3.15.HWE.raw", 
                map.file = "HESS_fix/pops/vcf.MMDP60.maf0.006.minQ30.recode.D3.15.HWE.map", quiet = FALSE,
                parallel = require("parallel"), n.cores = NULL)))
d3.14<-as.genind(as.matrix(read.PLINK("HESS_fix/pops/vcf.MMDP60.maf0.006.minQ30.recode.D3.14.HWE.raw", 
                  map.file = "HESS_fix/pops/vcf.MMDP60.maf0.006.minQ30.recode.D3.14.HWE.map", quiet = FALSE,
                  parallel = require("parallel"), n.cores = NULL)))
d3.04<-as.genind(as.matrix(read.PLINK("HESS_fix/pops/vcf.MMDP60.maf0.006.minQ30.recode.D3.04.HWE.raw", 
                  map.file = "HESS_fix/pops/vcf.MMDP60.maf0.006.minQ30.recode.D3.04.HWE.map", quiet = FALSE,
                  parallel = require("parallel"), n.cores = NULL)))

d2.15<-as.genind(as.matrix(read.PLINK("HESS_fix/pops/vcf.MMDP60.maf0.006.minQ30.recode.D2.15.HWE.raw", 
                  map.file = "HESS_fix/pops/vcf.MMDP60.maf0.006.minQ30.recode.D2.15.HWE.map", quiet = FALSE,
                  parallel = require("parallel"), n.cores = NULL)))
d2.14<-as.genind(as.matrix(read.PLINK("HESS_fix/pops/vcf.MMDP60.maf0.006.minQ30.recode.D2.14.HWE.raw", 
                  map.file = "HESS_fix/pops/vcf.MMDP60.maf0.006.minQ30.recode.D2.14.HWE.map", quiet = FALSE,
                  parallel = require("parallel"), n.cores = NULL)))
d2.04<-as.genind(as.matrix(read.PLINK("HESS_fix/pops/vcf.MMDP60.maf0.006.minQ30.recode.D2.04.HWE.raw", 
                  map.file = "HESS_fix/pops/vcf.MMDP60.maf0.006.minQ30.recode.D2.04.HWE.map", quiet = FALSE,
                  parallel = require("parallel"), n.cores = NULL)))

agr<-as.genind(as.matrix(read.PLINK("HESS_fix/pops/vcf.MMDP60.maf0.006.minQ30.recode.agraria.HWE.raw", 
                  map.file = "HESS_fix/pops/vcf.MMDP60.maf0.006.minQ30.recode.agraria.HWE.map", quiet = FALSE,
                  parallel = require("parallel"), n.cores = NULL)))
vil<-as.genind(as.matrix(read.PLINK("HESS_fix/pops/vcf.MMDP60.maf0.006.minQ30.recode.vilarinho.HWE.raw", 
                  map.file = "HESS_fix/pops/vcf.MMDP60.maf0.006.minQ30.recode.vilarinho.HWE.map", quiet = FALSE,
                  parallel = require("parallel"), n.cores = NULL)))
mira<-as.genind(as.matrix(read.PLINK("HESS_fix/pops/vcf.MMDP60.maf0.006.minQ30.recode.mira.HWE.raw", 
                  map.file = "HESS_fix/pops/vcf.MMDP60.maf0.006.minQ30.recode.mira.HWE.map", quiet = FALSE,
                  parallel = require("parallel"), n.cores = NULL)))


#AFLP
setwd("~/Dropbox/Genet Size/Jacob Analyses")

hier <- read.genalex("AFLP_genalex_AMANITABASED.csv")
hier.gi<-as.genind(as.matrix(hier))
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
pop(hier.gi)<-pops.aflp

roch1<-popsub(hier.gi, sublist = "Rochester_NY_2007_1")
roch2<-popsub(hier.gi, sublist = "Rochester_NY_2007_2")
roch3<-popsub(hier.gi, sublist = "Rochester_NY_2007_3")

ce02<-popsub(hier.gi, sublist = "CESAC_2002")
ce06<-popsub(hier.gi, sublist = "CESAC_2006")

d1<-popsub(hier.gi, sublist = "Drake_1_2004")
d2<-popsub(hier.gi, sublist = "Drake_2_2004")
d3<-popsub(hier.gi, sublist = "Drake_3_2004")
d4<-popsub(hier.gi, sublist = "Drake_4_2004")

hd1<-popsub(hier.gi, sublist = "Hearts_Desire_1_2004")
hd2<-popsub(hier.gi, sublist = "Hearts_Desire_2_2004")
hd3<-popsub(hier.gi, sublist = "Hearts_Desire_3_2004")

jl<-popsub(hier.gi, sublist = "Jakes_Landing_2006")
str<-popsub(hier.gi, sublist = "New_Jersey_str_2006")

S<-popsub(hier.gi, sublist = "Serbia_2007")

plot(upgma(dist(S)))


###########################################################
# coords
setwd("~/Dropbox/GENETS_finishing")
coordies<-read.csv("ALL_COORDS_geo_AMANITABASED.csv", head=T, stringsAsFactors = F)

roch1.c<-coordies[coordies$FINAlab %in% rownames(roch1$tab),]
roch2.c<-coordies[coordies$FINAlab %in% rownames(roch2$tab),]
roch3.c<-coordies[coordies$FINAlab %in% rownames(roch3$tab),]

ce02.c<-coordies[coordies$FINAlab %in% rownames(ce02$tab),]
ce06.c<-coordies[coordies$FINAlab %in% rownames(ce06$tab),]

d1.c<-coordies[coordies$FINAlab %in% rownames(d1$tab),]
d2.c<-coordies[coordies$FINAlab %in% c(rownames(d2$tab),rownames(d2.04$tab)),]
d3.c<-coordies[coordies$FINAlab %in% c(rownames(d3$tab),rownames(d3.04$tab)),]
d4.c<-coordies[coordies$FINAlab %in% rownames(d4$tab),]

hd1.c<-coordies[coordies$FINAlab %in% rownames(hd1$tab),]
hd2.c<-coordies[coordies$FINAlab %in% rownames(hd2$tab),]
hd3.c<-coordies[coordies$FINAlab %in% rownames(hd3$tab),]

jl.c<-coordies[coordies$FINAlab %in% rownames(jl$tab),]
jl.c$genet.1 <- sub("^$", "sing", jl.c$genet.1)
str.c<-coordies[coordies$FINAlab %in% rownames(str$tab),]

S.c<-coordies[coordies$FINAlab %in% rownames(S$tab),]
S.c[,5:6]<-S.c[,5:6]/100


mira.c<-coordies[coordies$FINAlab %in% rownames(mira$tab),]
agr.c<-coordies[coordies$FINAlab %in% rownames(agr$tab),]


d2.14.c<-coordies[coordies$FINAlab %in% c(rownames(d2.14$tab)),]
d2.15.c<-coordies[coordies$FINAlab %in% c(rownames(d2.15$tab)),]

d3.14.c<-coordies[coordies$FINAlab %in% c(rownames(d3.14$tab)),]
d3.15.c<-coordies[coordies$FINAlab %in% c(rownames(d3.15$tab)),]


#################################
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(grid)
library(ggtree)
############ Now I can plot, haliluya
mytheme <- theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
                 panel.grid.major = element_line(colour="white",size = rel(0.5)),
                 plot.background = element_rect(fill = "white"),
                 panel.background = element_rect(fill = "white"),
                 legend.position="none")

ggplot(d2.c, aes(x=x, y=y, label=ID, color=genet, shape=shape))+
  geom_point(size=3)+mytheme+labs(x="Distance (m)", y="Distance (m)")+
  geom_text_repel(aes(label=ID),hjust=0, vjust=0)+
  scale_shape_manual(values=c("ui"=15, "gen"=19, "tri"=17))+
  scale_color_manual(values=c("Unknown"="black", "genet"="darkgreen"))+coord_fixed(ratio = 1)

ggplot(d2.c, aes(x=x, y=y, label=ID, color=genet, shape=shape))+
  geom_point(size=3)+mytheme+labs(x="Distance (m)", y="Distance (m)")+
  geom_text_repel(aes(label=ID),hjust=0, vjust=0)+
  scale_shape_manual(values=c("ui"=15, "gen"=19, "tri"=17))+
  scale_color_manual(values=c("Unknown"="black", "genet"="darkgreen"))+coord_fixed(ratio = 1)+
  scale_x_continuous(limits=c(-15,23))+scale_y_continuous(limits=c(-8, 22))

d2.g.pre<-ggplot(d2.c, aes(x=x, y=y, label=ID, color=genet, shape=shape))+
  geom_point(size=1)+mytheme+labs(x="Distance (m)", y="Distance (m)")+
  geom_text_repel(aes(label=ID),hjust=0, vjust=0, size=2, segment.size = .1, segment.alpha = 1)+
  scale_shape_manual(values=c("ui"=15, "gen"=19, "tri"=17))+
  scale_color_manual(values=c("Unknown"="black", "genet"="darkgreen"))+ggtitle("Drake 2 2004")

d2snp.g.pre<-ggplot(d2.c, aes(x=x, y=y, label=ID, color=genet, shape=shape))+
  geom_point(size=1)+mytheme+labs(x="Distance (m)", y="Distance (m)")+
  geom_text_repel(aes(label=ID),hjust=0, vjust=0, size=1.5, segment.size = .05, segment.alpha = 1)+
  scale_shape_manual(values=c("ui"=15, "gen"=19, "tri"=17))+
  scale_color_manual(values=c("Unknown"="black", "genet"="darkgreen"))+ggtitle("Drake 2 2004")




# d2.p<-ggtree(d2.p.pre, color="#9400D3")+geom_tiplab(hjust=.5, angle=0)+scale_x_reverse()+coord_flip()+geom_tippoint(shape=18,size=3, col="black")
# d2.04.p <- ggtree(d2.04.p.pre, color="#FFA500")+geom_tiplab(hjust=.5, angle=-45)+scale_x_reverse()+geom_tippoint(shape=18,size=3, col="black")
# d2.g<-list(ggplotGrob(d2.g.pre),ggplotGrob(d2.p), ggplotGrob(d2.04.p))
# lay <- rbind(c(2,2,2, 2, NA),c(1,1,1,1,3),c(1,1,1,1,3), c(1,1,1,1,3), c(1,1,1,1,3))
# gaga.d2.14<- grid.arrange(grobs = d2.g, layout_matrix = lay)

#d1
d1.g.pre<-ggplot(d1.c, aes(x=x, y=y, label=ID, color=genet.1, shape=genet.1))+
  geom_point(size=1)+mytheme+labs(x="Distance (m)", y="Distance (m)")+
  geom_text_repel(aes(label=ID),hjust=0, vjust=0,size=1.5, segment.size = .05, segment.alpha = 1)+
  scale_shape_manual(values=c("ui"=15, "d1A"=19, "d1B"=19, "d1C"=19, "d1D"=19, "d1E"=19, "d1F"=19))+
  scale_color_manual(values=c("d1A"="dodgerblue",
                              "d1B"="firebrick",
                              "d1C"="forestgreen",
                              "d1D"="darkslateblue",
                              "d1E"="darkorange2",
                              "d1F"="darkmagenta",
                              "ui"="black"))+ggtitle("Drake 1 2004")


#d4
d4.c$genet.1 <- sub("^$", "ui", d4.c$genet.1)

d4.g.pre<-ggplot(d4.c, aes(x=x, y=y, label=ID, color=genet.1, shape=genet.1))+
  geom_point(size=1)+mytheme+labs(x="Distance (m)", y="Distance (m)")+
  geom_text_repel(aes(label=ID),hjust=0, vjust=0, size=1.5, segment.size = .05, segment.alpha = 1)+
  scale_shape_manual(values=c("ui"=15, "d4.A"=19))+
  scale_color_manual(values=c("d4.A"="dodgerblue",
                              "ui"="black"))+ggtitle("Drake 4 2004")




# d3


ggplot(d3.c, aes(x=x, y=y, label=ID, color=genet, shape=shape))+
  geom_point(size=3)+mytheme+labs(x="Distance (m)", y="Distance (m)")+
  geom_text_repel(aes(label=ID),hjust=0, vjust=0)+
  scale_shape_manual(values=c("ui"=15, "gen"=19, "tri"=17))+
  scale_color_manual(values=c("Unknown"="black", "genet"="darkgreen"))+coord_fixed(ratio=1)




d3.g.pre<-ggplot(d3.c, aes(x=x, y=y, label=ID, color=genet, shape=shape))+
  geom_point(size=1)+mytheme+labs(x="Distance (m)", y="Distance (m)")+
  geom_text_repel(aes(label=ID),hjust=0, vjust=0, size=1.5, segment.size = .05, segment.alpha = 1)+
  scale_shape_manual(values=c("ui"=15, "gen"=19, "tri"=17))+
  scale_color_manual(values=c("Unknown"="black", "genet"="darkgreen"))+ggtitle("Drake 3 2004")


d3snp.g.pre<-ggplot(d3.c[2:nrow(d3.c),], aes(x=x, y=y, label=ID, color=genet, shape=shape))+
  geom_point(size=2)+mytheme+labs(x="Distance (m)", y="Distance (m)")+
  geom_text_repel(aes(label=ID),hjust=0, vjust=0, size=1.5, segment.size = .05, segment.alpha = 1)+
  scale_shape_manual(values=c("ui"=15, "gen"=15, "tri"=15))+
  scale_color_manual(values=c("Unknown"="black", "genet"="black"))+ggtitle("Drake 3 2004")


d3.p.pre<-upgma(dist(d3))
d3.04.p.pre<-upgma(dist(d3.04))

# d3.p<-ggtree(d3.p.pre, color="#9400D3")+geom_tiplab(hjust=.5, angle=0)+scale_x_reverse()+coord_flip()+geom_tippoint(shape=18,size=3, col="black")
# d3.04.p <- ggtree(d3.04.p.pre, color="#FFA500")+geom_tiplab(hjust=.5, angle=-45)+scale_x_reverse()+geom_tippoint(shape=18,size=3, col="black")
# d3.g<-list(ggplotGrob(d3.g.pre),ggplotGrob(d3.p), ggplotGrob(d3.04.p))
# lay <- rbind(c(2,2,2, 2, NA),c(1,1,1,1,3),c(1,1,1,1,3), c(1,1,1,1,3), c(1,1,1,1,3))
# gaga.d3.04<- grid.arrange(grobs = d3.g, layout_matrix = lay)

# d2.14
ggplot(d2.14.c, aes(x=x, y=y, label=ID, fill="black"))+
  geom_point(size=3, shape=15)+mytheme+labs(x="Distance (m)", y="Distance (m)")+
  geom_text_repel(aes(label=ID),hjust=0, vjust=0)+coord_fixed(ratio = 1)

ggplot(d2.14.c, aes(x=x, y=y, label=ID, fill="black"))+
  geom_point(size=3, shape=15)+mytheme+labs(x="Distance (m)", y="Distance (m)")+
  geom_text_repel(aes(label=ID),hjust=0, vjust=0)+coord_fixed(ratio = 1)+
  scale_x_continuous(limits=c(-15,23))+scale_y_continuous(limits=c(-8,22))

d2.14.g.pre<-ggplot(d2.14.c, aes(x=x, y=y, label=ID, fill="black"))+
  geom_point(size=2, shape=15)+mytheme+labs(x="Distance (m)", y="Distance (m)")+
  geom_text_repel(aes(label=ID),hjust=0, vjust=0, size=1.5, segment.size = .05, segment.alpha = 1)+ggtitle("Drake 2 2014")

# d2.14.p.pre<-upgma(dist(d2.14))
# 
# r<-rectGrob(gp=gpar(fill="white", col="white"))
# d2.14.p <- ggtree(d2.14.p.pre, color="#FFA500")+geom_tiplab(hjust=.5, angle=-45, size=3)+scale_x_reverse()+geom_tippoint(shape=18,size=2, col="black")
# d2.14.g<-list(ggplotGrob(d2.14.g.pre),(r), ggplotGrob(d2.14.p))
# lay <- rbind(c(2,2,2, 2, NA),c(1,1,1,1,3),c(1,1,1,1,3), c(1,1,1,1,3), c(1,1,1,1,3))
# gaga.d2.14<- grid.arrange(grobs = d2.14.g, layout_matrix = lay)


#d3.14

ggplot(d3.14.c, aes(x=x, y=y, label=ID, fill="black"))+
  geom_point(size=3, shape=15)+mytheme+labs(x="Distance (m)", y="Distance (m)")+
  geom_text_repel(aes(label=ID),hjust=0, vjust=0)+coord_fixed(ratio=1)


d3.14.g.pre<-ggplot(d3.14.c, aes(x=x, y=y, label=ID, fill="black"))+
  geom_point(size=2, shape=15)+mytheme+labs(x="Distance (m)", y="Distance (m)")+
  geom_text_repel(aes(label=ID),hjust=0, vjust=0,size=1.5, segment.size = .05, segment.alpha = 1)+ggtitle("Drake 3 2014")

# d3.14.p.pre<-upgma(dist(d3.14))
# 
# r<-rectGrob(gp=gpar(fill="white", col="white"))
# d3.14.p <- ggtree(d3.14.p.pre, color="#FFA500")+geom_tiplab(hjust=.5, angle=-45, size=3)+scale_x_reverse()+geom_tippoint(shape=18,size=2, col="black")
# d3.14.g<-list(ggplotGrob(d3.14.g.pre),(r), ggplotGrob(d3.14.p))
# lay <- rbind(c(2,2,2, 2, NA),c(1,1,1,1,3),c(1,1,1,1,3), c(1,1,1,1,3), c(1,1,1,1,3))
# gaga.d3.14<- grid.arrange(grobs = d3.14.g, layout_matrix = lay)


# d2.15


ggplot(d2.15.c, aes(x=x, y=y, label=ID, fill="black"))+
  geom_point(size=3, shape=15)+mytheme+labs(x="Distance (m)", y="Distance (m)")+
  geom_text_repel(aes(label=ID),hjust=0, vjust=0)+coord_fixed(ratio = 1)

ggplot(d2.15.c, aes(x=x, y=y, label=ID, fill="black"))+
  geom_point(size=2, shape=15)+mytheme+labs(x="Distance (m)", y="Distance (m)")+
  geom_text_repel(aes(label=ID),hjust=0, vjust=0,size=1.5, segment.size = .05, segment.alpha = 1)+coord_fixed(ratio = 1)+
  scale_x_continuous(limits=c(-15,23))+scale_y_continuous(limits=c(-8,22))


d2.15.g.pre<-ggplot(d2.15.c, aes(x=x, y=y, label=ID, fill="black"))+
  geom_point(size=2, shape=15)+mytheme+labs(x="Distance (m)", y="Distance (m)")+
  geom_text_repel(aes(label=ID),hjust=0, vjust=0,size=1.5, segment.size = .05, segment.alpha = 1)+ggtitle("Drake 2 2015")

# d2.15.p.pre<-upgma(dist(d2.15))
# 
# r<-rectGrob(gp=gpar(fill="white", col="white"))
# d2.15.p <- ggtree(d2.15.p.pre, color="#FFA500")+geom_tiplab(hjust=.5, angle=-45, size=3)+scale_x_reverse()+geom_tippoint(shape=18,size=2, col="black")
# d2.15.g<-list(ggplotGrob(d2.15.g.pre),(r), ggplotGrob(d2.15.p))
# lay <- rbind(c(2,2,2, 2, NA),c(1,1,1,1,3),c(1,1,1,1,3), c(1,1,1,1,3), c(1,1,1,1,3))
# gaga.d2.15<- grid.arrange(grobs = d2.15.g, layout_matrix = lay)


#d3.15

ggplot(d3.15.c, aes(x=x, y=y, label=ID, fill="black"))+
  geom_point(size=3, shape=15)+mytheme+labs(x="Distance (m)", y="Distance (m)")+
  geom_text_repel(aes(label=ID),hjust=0, vjust=0)+coord_fixed(ratio=1)

d3.15.g.pre<-ggplot(d3.15.c, aes(x=x, y=y, label=ID, fill="black"))+
  geom_point(size=2, shape=15)+mytheme+labs(x="Distance (m)", y="Distance (m)")+
  geom_text_repel(aes(label=ID),hjust=0, vjust=0,size=1.5, segment.size = .05, segment.alpha = 1)+ggtitle("Drake 3 2015")

# d3.15.p.pre<-upgma(dist(d3.15))
# 
# r<-rectGrob(gp=gpar(fill="white", col="white"))
# d3.15.p <- ggtree(d3.15.p.pre, color="#FFA500")+geom_tiplab(hjust=.5, angle=-45, size=3)+scale_x_reverse()+geom_tippoint(shape=18,size=2, col="black")
# d3.15.g<-list(ggplotGrob(d3.15.g.pre),(r), ggplotGrob(d3.15.p))
# lay <- rbind(c(2,2,2, 2, NA),c(1,1,1,1,3),c(1,1,1,1,3), c(1,1,1,1,3), c(1,1,1,1,3))
# gaga.d3.15<- grid.arrange(grobs = d3.15.g, layout_matrix = lay)

#agr

ggplot(agr.c, aes(x=x, y=y, label=ID, fill="black"))+
  geom_point(size=3, shape=15)+mytheme+labs(x="Distance (m)", y="Distance (m)")+
  geom_text_repel(aes(label=ID),hjust=0, vjust=0)+coord_fixed(ratio = 1, ylim=c(-1,1))

agr.g.pre<-ggplot(agr.c, aes(x=x, y=y, label=ID, fill="black"))+
  geom_point(size=2, shape=15)+mytheme+labs(x="Distance (m)", y="Distance (m)")+
  geom_text_repel(aes(label=ID),hjust=0, vjust=0,size=1.5, segment.size = .05, segment.alpha = 1)+ggtitle("Agraria 2015")

# agr.p.pre<-upgma(dist(agr))
# 
# r<-rectGrob(gp=gpar(fill="white", col="white"))
# agr.p <- ggtree(agr.p.pre, color="#FFA500")+geom_tiplab(hjust=.5, angle=-45, size=3)+scale_x_reverse()+geom_tippoint(shape=18,size=2, col="black")
# agr.g<-list(ggplotGrob(agr.g.pre),(r), ggplotGrob(agr.p))
# lay <- rbind(c(2,2,2, 2, NA),c(1,1,1,1,3),c(1,1,1,1,3), c(1,1,1,1,3), c(1,1,1,1,3))
# gaga.agr<- grid.arrange(grobs = agr.g, layout_matrix = lay)

#mira

ggplot(mira.c, aes(x=x, y=y, label=ID, fill="black"))+
  geom_point(size=3, shape=15)+mytheme+labs(x="Distance (m)", y="Distance (m)")+
  geom_text_repel(aes(label=ID),hjust=0, vjust=0)+coord_fixed(ratio=1, ylim=c(-1,1))

mira.g.pre<-ggplot(mira.c, aes(x=x, y=y, label=ID, fill="black"))+
  geom_point(size=2, shape=15)+mytheme+labs(x="Distance (m)", y="Distance (m)")+
  geom_text_repel(aes(label=ID),hjust=0, vjust=0,size=1.5, segment.size = .05, segment.alpha = 1)+ggtitle("Mira 2015")

# mira.p.pre<-upgma(dist(mira))
# 
# r<-rectGrob(gp=gpar(fill="white", col="white"))
# mira.p <- ggtree(mira.p.pre, color="#FFA500")+geom_tiplab(hjust=.5, angle=-45, size=3)+scale_x_reverse()+geom_tippoint(shape=18,size=2, col="black")
# mira.g<-list(ggplotGrob(mira.g.pre),(r), ggplotGrob(mira.p))
# lay <- rbind(c(2,2,2, 2, NA),c(1,1,1,1,3),c(1,1,1,1,3), c(1,1,1,1,3), c(1,1,1,1,3))
# gaga.mira<- grid.arrange(grobs = mira.g, layout_matrix = lay)

#CESAC02

ggplot(ce02.c, aes(x=x, y=y, label=FINAlab, color=genet.1, shape=shape))+
  mytheme+labs(x="Distance (m)", y="Distance (m)")+
  geom_text_repel(aes(label=FINAlab),hjust=0, vjust=0, size=1.5)+
  geom_point(size=1)+
  scale_shape_manual(values=c("ui"=15, "gen"=19, "conf"=17))+
  scale_color_manual(values=c("F"="dodgerblue",
                              "D"="firebrick",
                              "Singleton"="black"))+coord_fixed(ratio=1)

ggplot(ce02.c, aes(x=x, y=y, label=FINAlab, color=genet.1, shape=shape))+
  mytheme+labs(x="Distance (m)", y="Distance (m)")+
  geom_text_repel(aes(label=FINAlab),hjust=0, vjust=0, size=2.5)+
  geom_point(size=3)+
  scale_shape_manual(values=c("ui"=15, "gen"=19, "conf"=17))+
  scale_color_manual(values=c("F"="dodgerblue",
                              "D"="firebrick",
                              "Singleton"="black"))+coord_fixed(ratio=1)+
  scale_x_continuous(limits=c(-10, 15)) +
  scale_y_continuous(limits=c(-10, 15))


ce02.g.pre<-ggplot(ce02.c, aes(x=x, y=y, label=FINAlab, color=genet.1, shape=shape))+
  mytheme+labs(x="Distance (m)", y="Distance (m)")+
  geom_text_repel(aes(label=FINAlab),hjust=0, vjust=0, size=1.5, segment.size = .05, segment.alpha = 1)+
  geom_point(size=1)+
  scale_shape_manual(values=c("ui"=15, "gen"=19, "conf"=17))+
  scale_color_manual(values=c("F"="dodgerblue",
                              "D"="firebrick",
                              "Singleton"="black"))+ggtitle("CESAC 2002")

# ce02.p.pre<-upgma(dist(ce02))
# 
# ce02.p<-ggtree(ce02.p.pre, color="#9400D3")+geom_tiplab(hjust=.5, angle=-45, size=2)+scale_x_reverse()+coord_flip()+geom_tippoint(shape=18,size=3, col="black")
# r<-rectGrob(gp=gpar(fill="white", col="white"))
# ce02.g<-list(ggplotGrob(ce02.g.pre), ggplotGrob(ce02.p), r)
# lay <- rbind(c(2,2,2, 2, NA),c(1,1,1,1,3),c(1,1,1,1,3), c(1,1,1,1,3), c(1,1,1,1,3))
# gaga.ce02<- grid.arrange(grobs = ce02.g, layout_matrix = lay)

#CESAC06

ggplot(ce06.c, aes(x=x, y=y, label=FINAlab, color=genet.1, shape=shape))+
  mytheme+labs(x="Distance (m)", y="Distance (m)")+
  geom_text_repel(aes(label=FINAlab),hjust=0, vjust=0, size=2)+
  geom_point(size=2)+
  scale_shape_manual(values=c("ui"=15, "gen"=19, "conf"=17))+
  scale_color_manual(values=c("A"="dodgerblue",
                              "B"="firebrick",
                              "C"="forestgreen",
                              "D"="darkslateblue",
                              "E"="darkorange2",
                              "Singleton"="black"))+coord_fixed(ratio=1)

ggplot(ce06.c, aes(x=x, y=y, label=FINAlab, color=genet.1, shape=shape))+
  mytheme+labs(x="Distance (m)", y="Distance (m)")+
  geom_text_repel(aes(label=FINAlab),hjust=0, vjust=0, size=2)+
  geom_point(size=2)+
  scale_shape_manual(values=c("ui"=15, "gen"=19, "conf"=17))+
  scale_color_manual(values=c("A"="dodgerblue",
                              "B"="firebrick",
                              "C"="forestgreen",
                              "D"="darkslateblue",
                              "E"="darkorange2",
                              "Singleton"="black"))+coord_fixed(ratio=1)+
  scale_x_continuous(limits=c(-10, 15)) +
  scale_y_continuous(limits=c(-10, 15))

ce06.g.pre<-ggplot(ce06.c, aes(x=x, y=y, label=FINAlab, color=genet.1, shape=shape))+
  mytheme+labs(x="Distance (m)", y="Distance (m)")+
  geom_text_repel(aes(label=FINAlab),hjust=0, vjust=0, size=1.5, segment.size = .05, segment.alpha = 1)+
  geom_point(size=1)+
  scale_shape_manual(values=c("ui"=15, "gen"=19, "conf"=17))+
  scale_color_manual(values=c("A"="dodgerblue",
                              "B"="firebrick",
                              "C"="forestgreen",
                              "D"="darkslateblue",
                              "E"="darkorange2",
                              "Singleton"="black"))+ggtitle("CESAC 2006")

# ce06.p.pre<-upgma(dist(ce06))
# 
# ce06.p<-ggtree(ce06.p.pre, color="#9400D3")+geom_tiplab(hjust=.5, angle=-45, size=2)+scale_x_reverse()+coord_flip()+geom_tippoint(shape=18,size=3, col="black")
# r<-rectGrob(gp=gpar(fill="white", col="white"))
# ce06.g<-list(ggplotGrob(ce06.g.pre), ggplotGrob(ce06.p), r)
# lay <- rbind(c(2,2,2, 2, NA),c(1,1,1,1,3),c(1,1,1,1,3), c(1,1,1,1,3), c(1,1,1,1,3))
# gaga.ce06<- grid.arrange(grobs = ce06.g, layout_matrix = lay)



#jl

ggplot(jl.c, aes(x=x, y=y, label=FINAlab, color=genet.1, shape=shape))+
  geom_point(size=1)+mytheme+labs(x="Distance (m)", y="Distance (m)")+
  geom_text_repel(aes(label=FINAlab),hjust=0, vjust=0, size=1.5, segment.size = .05, segment.alpha = 1)+
  scale_shape_manual(values=c("ui"=15, "gen"=19, "conf"=17))+
  scale_color_manual(values=c("A"="dodgerblue",
                              "B"="firebrick",
                              "C"="forestgreen",
                              "D"="gold3",
                              "E"="darkolivegreen",
                              "F"="darkorange2",
                              "G"="darkmagenta",
                              "H"="darkslateblue",
                              "I"="chartreuse3",
                              "sing"="black"))+coord_fixed(ratio=1)



jl.g.pre<-ggplot(jl.c, aes(x=x, y=y, label=FINAlab, color=genet.1, shape=shape))+
  geom_point(size=1)+mytheme+labs(x="Distance (m)", y="Distance (m)")+
  geom_text_repel(aes(label=FINAlab),hjust=0, vjust=0, size=1.5, segment.size = .05, segment.alpha = 1)+
  scale_shape_manual(values=c("ui"=15, "gen"=19, "conf"=17))+
  scale_color_manual(values=c("A"="dodgerblue",
                              "B"="firebrick",
                              "C"="forestgreen",
                              "D"="gold3",
                              "E"="darkolivegreen",
                              "F"="darkorange2",
                              "G"="darkmagenta",
                              "H"="darkslateblue",
                              "I"="chartreuse3",
                              "sing"="black"))+ggtitle("JL 2006")

# jl.p.pre<-upgma(dist(jl))
# 
# jl.p<-ggtree(jl.p.pre, color="#9400D3")+geom_tiplab(hjust=.5, angle=-45, size=2)+scale_x_reverse()+coord_flip()+geom_tippoint(shape=18,size=3, col="black")
# r<-rectGrob(gp=gpar(fill="white", col="white"))
# jl.g<-list(ggplotGrob(jl.g.pre), ggplotGrob(jl.p), r)
# lay <- rbind(c(2,2,2, 2, NA),c(1,1,1,1,3),c(1,1,1,1,3), c(1,1,1,1,3), c(1,1,1,1,3))
# gaga.jl<- grid.arrange(grobs = jl.g, layout_matrix = lay)


#str

ggplot(str.c, aes(x=x, y=y, label=FINAlab, color=genet.1, shape=genet.1))+
  geom_point(size=3)+mytheme+labs(x="Distance (m)", y="Distance (m)")+
  geom_text_repel(aes(label=FINAlab),hjust=0, vjust=0)+
  scale_shape_manual(values=c("A"=19, "Singleton"=15))+
  scale_color_manual(values=c("A"="dodgerblue",
                              "Singleton"="black"))+coord_fixed(ratio=1)


str.g.pre<-ggplot(str.c, aes(x=x, y=y, label=FINAlab, color=genet.1, shape=genet.1))+
  geom_point(size=1)+mytheme+labs(x="Distance (m)", y="Distance (m)")+
  geom_text_repel(aes(label=FINAlab),hjust=0, vjust=0, size=1.5, segment.size = .05, segment.alpha = 1)+
  scale_shape_manual(values=c("A"=19, "Singleton"=15))+
  scale_color_manual(values=c("A"="dodgerblue",
                              "Singleton"="black"))+ggtitle("str 2006")

# str.p.pre<-upgma(dist(str))
# 
# str.p<-ggtree(str.p.pre, color="#9400D3")+geom_tiplab(hjust=.5, angle=-45, size=3)+scale_x_reverse()+coord_flip()+geom_tippoint(shape=18,size=3, col="black")
# r<-rectGrob(gp=gpar(fill="white", col="white"))
# str.g<-list(ggplotGrob(str.g.pre), ggplotGrob(str.p), r)
# lay <- rbind(c(2,2,2, 2, NA),c(1,1,1,1,3),c(1,1,1,1,3), c(1,1,1,1,3), c(1,1,1,1,3))
# gaga.str<- grid.arrange(grobs = str.g, layout_matrix = lay)

#Serbia

ggplot(S.c, aes(x=x, y=y, label=FINAlab, color=genet.1, shape=genet.1))+
  geom_point(size=1)+mytheme+labs(x="Distance (m)", y="Distance (m)")+
  geom_text_repel(aes(label=FINAlab),hjust=0, vjust=0, size=1.5, segment.size = .05, segment.alpha = 1)+
  scale_shape_manual(values=c("SA"=19, "SB"=19, "SC"=19, "ui"=15))+
  scale_color_manual(values=c("SA"="dodgerblue",
                              "SB"="firebrick",
                              "SC"="forestgreen",
                              "ui"="black"))+ggtitle("S 2007")

S.g.pre<-ggplot(S.c, aes(x=x, y=y, label=FINAlab, color=genet.1, shape=genet.1))+
  geom_point(size=1)+mytheme+labs(x="Distance (m)", y="Distance (m)")+
  geom_text_repel(aes(label=FINAlab),hjust=0, vjust=0, size=1.5, segment.size = .05, segment.alpha = 1)+
  scale_shape_manual(values=c("SA"=19, "SB"=19, "SC"=19, "ui"=15))+
  scale_color_manual(values=c("SA"="dodgerblue",
                              "SB"="firebrick",
                              "SC"="forestgreen",
                              "ui"="black"))+ggtitle("S 2007")

# S.p.pre<-upgma(dist(S))
# 
# S.p<-ggtree(S.p.pre, color="#9400D3")+geom_tiplab(hjust=.5, angle=-45, size=3)+scale_x_reverse()+coord_flip()+geom_tippoint(shape=18,size=3, col="black")
# r<-rectGrob(gp=gpar(fill="white", col="white"))
# S.g<-list(ggplotGrob(S.g.pre), ggplotGrob(S.p), r)
# lay <- rbind(c(2,2,2, 2, NA),c(1,1,1,1,3),c(1,1,1,1,3), c(1,1,1,1,3), c(1,1,1,1,3))
# gaga.S<- grid.arrange(grobs = S.g, layout_matrix = lay)

#HD1

ggplot(hd1.c, aes(x=x, y=y, label=FINAlab, color=shape, shape=shape))+
  geom_point(size=3)+mytheme+labs(x="Distance (m)", y="Distance (m)")+
  geom_text_repel(aes(label=FINAlab),hjust=0, vjust=0)+
  scale_shape_manual(values=c("ui"=15))+
  scale_color_manual(values=c("ui"="black"))+coord_fixed(ratio = 1)


hd1.g.pre<-ggplot(hd1.c, aes(x=x, y=y, label=FINAlab, color=shape, shape=shape))+
  geom_point(size=1)+mytheme+labs(x="Distance (m)", y="Distance (m)")+
  geom_text_repel(aes(label=FINAlab),hjust=0, vjust=0, size=1.5, segment.size = .05, segment.alpha = 1)+
  scale_shape_manual(values=c("ui"=19))+
  scale_color_manual(values=c("ui"="black"))+ggtitle("HD1 2004")

# hd1.p.pre<-upgma(dist(hd1))
# 
# hd1.p<-ggtree(hd1.p.pre, color="#9400D3")+geom_tiplab(hjust=.5, angle=-45, size=3)+scale_x_reverse()+coord_flip()+geom_tippoint(shape=18,size=3, col="black")
# r<-rectGrob(gp=gpar(fill="white", col="white"))
# hd1.g<-list(ggplotGrob(hd1.g.pre), ggplotGrob(hd1.p), r)
# lay <- rbind(c(2,2,2, 2, NA),c(1,1,1,1,3),c(1,1,1,1,3), c(1,1,1,1,3), c(1,1,1,1,3))
# gaga.hd1<- grid.arrange(grobs = hd1.g, layout_matrix = lay)

#HD2

ggplot(hd2.c, aes(x=x, y=y, label=FINAlab, color=genet.1, shape=genet.1))+
  geom_point(size=3)+mytheme+labs(x="Distance (m)", y="Distance (m)")+
  geom_text_repel(aes(label=FINAlab),hjust=0, vjust=0)+
  scale_shape_manual(values=c("A"=19, "B" =19,"ui"=15))+
  scale_color_manual(values=c("A"="dodgerblue","B"="firebrick","ui"="black"))+coord_fixed(ratio = 1)

hd2.g.pre<-ggplot(hd2.c, aes(x=x, y=y, label=FINAlab, color=genet.1, shape=genet.1))+
  geom_point(size=1)+mytheme+labs(x="Distance (m)", y="Distance (m)")+
  geom_text_repel(aes(label=FINAlab),hjust=0, vjust=0, size=1.5, segment.size = .05, segment.alpha = 1)+
  scale_shape_manual(values=c("A"=19, "B" =19,"ui"=15))+
  scale_color_manual(values=c("A"="dodgerblue","B"="firebrick","ui"="black"))+ggtitle("HD2 2004")

# hd2.p.pre<-upgma(dist(hd2))
# 
# hd2.p<-ggtree(hd2.p.pre, color="#9400D3")+geom_tiplab(hjust=.5, angle=-45, size=3)+scale_x_reverse()+coord_flip()+geom_tippoint(shape=18,size=3, col="black")
# r<-rectGrob(gp=gpar(fill="white", col="white"))
# hd2.g<-list(ggplotGrob(hd2.g.pre), ggplotGrob(hd2.p), r)
# lay <- rbind(c(2,2,2, 2, NA),c(1,1,1,1,3),c(1,1,1,1,3), c(1,1,1,1,3), c(1,1,1,1,3))
# gaga.hd2<- grid.arrange(grobs = hd2.g, layout_matrix = lay)

#HD3

ggplot(hd3.c, aes(x=x, y=y, label=FINAlab, color=genet.1, shape=genet.1))+
  geom_point(size=3)+mytheme+labs(x="Distance (m)", y="Distance (m)")+
  geom_text_repel(aes(label=FINAlab),hjust=0, vjust=0)+
  scale_shape_manual(values=c("A"=19, "ui"=15))+
  scale_color_manual(values=c("A"="dodgerblue","ui"="black"))+coord_fixed(ratio=1)

hd3.g.pre<-ggplot(hd3.c, aes(x=x, y=y, label=FINAlab, color=genet.1, shape=genet.1))+
  geom_point(size=1)+mytheme+labs(x="Distance (m)", y="Distance (m)")+
  geom_text_repel(aes(label=FINAlab),hjust=0, vjust=0, size=1.5, segment.size = .05, segment.alpha = 1)+
  scale_shape_manual(values=c("A"=19, "ui"=15))+
  scale_color_manual(values=c("A"="dodgerblue","ui"="black"))+ggtitle("HD3 2004")

# hd3.p.pre<-upgma(dist(hd3))
# 
# hd3.p<-ggtree(hd3.p.pre, color="#9400D3")+geom_tiplab(hjust=.5, angle=-45, size=3)+scale_x_reverse()+coord_flip()+geom_tippoint(shape=18,size=3, col="black")
# r<-rectGrob(gp=gpar(fill="white", col="white"))
# hd3.g<-list(ggplotGrob(hd3.g.pre), ggplotGrob(hd3.p), r)
# lay <- rbind(c(2,2,2, 2, NA),c(1,1,1,1,3),c(1,1,1,1,3), c(1,1,1,1,3), c(1,1,1,1,3))
# gaga.hd3<- grid.arrange(grobs = hd3.g, layout_matrix = lay)


#roch1

ggplot(roch1.c, aes(x=x, y=y, label=FINAlab, color=genet.1, shape=shape))+
  geom_point(size=2)+mytheme+labs(x="Distance (m)", y="Distance (m)")+
  geom_text_repel(aes(label=FINAlab),hjust=0, vjust=0, size=3)+
  scale_shape_manual(values=c("gen"=19, "ui"=15))+
  scale_color_manual(values=c("A1"="dodgerblue","ui"="black"))+coord_fixed(ratio = 1)

roch1.g.pre<-ggplot(roch1.c, aes(x=x, y=y, label=FINAlab, color=genet.1, shape=shape))+
  geom_point(size=1)+mytheme+labs(x="Distance (m)", y="Distance (m)")+
  geom_text_repel(aes(label=FINAlab),hjust=0, vjust=0, size=1.5, segment.size = .05, segment.alpha = 1)+
  scale_shape_manual(values=c("gen"=19, "ui"=15))+
  scale_color_manual(values=c("A1"="dodgerblue","ui"="black"))+ggtitle("Roch1 2007")

# roch1.p.pre<-upgma(dist(roch1))
# 
# roch1.p<-ggtree(roch1.p.pre, color="#9400D3")+geom_tiplab(hjust=.5, angle=-45, size=3)+scale_x_reverse()+coord_flip()+geom_tippoint(shape=18,size=3, col="black")
# r<-rectGrob(gp=gpar(fill="white", col="white"))
# roch1.g<-list(ggplotGrob(roch1.g.pre), ggplotGrob(roch1.p), r)
# lay <- rbind(c(2,2,2, 2, NA),c(1,1,1,1,3),c(1,1,1,1,3), c(1,1,1,1,3), c(1,1,1,1,3))
# gaga.hd3<- grid.arrange(grobs = roch1.g, layout_matrix = lay)

#roch2

ggplot(roch2.c, aes(x=x, y=y, label=FINAlab, color=genet.1, shape=shape))+
  geom_point(size=3)+mytheme+labs(x="Distance (m)", y="Distance (m)")+
  geom_text_repel(aes(label=FINAlab),hjust=0, vjust=0)+
  scale_shape_manual(values=c("gen"=19, "ui"=15))+
  scale_color_manual(values=c("A2"="dodgerblue","ui"="black"))+coord_fixed(ratio = 1)

roch2.g.pre<-ggplot(roch2.c, aes(x=x, y=y, label=FINAlab, color=genet.1, shape=shape))+
  geom_point(size=1)+mytheme+labs(x="Distance (m)", y="Distance (m)")+
  geom_text_repel(aes(label=FINAlab),hjust=0, vjust=0, size=1.5, segment.size = .05, segment.alpha = 1)+
  scale_shape_manual(values=c("gen"=19, "ui"=15))+
  scale_color_manual(values=c("A2"="dodgerblue","ui"="black"))+ggtitle("Roch2 2007")

# roch2.p.pre<-upgma(dist(roch2))
# 
# roch2.p<-ggtree(roch2.p.pre, color="#9400D3")+geom_tiplab(hjust=.5, angle=-45, size=3)+scale_x_reverse()+coord_flip()+geom_tippoint(shape=18,size=3, col="black")
# r<-rectGrob(gp=gpar(fill="white", col="white"))
# roch2.g<-list(ggplotGrob(roch2.g.pre), ggplotGrob(roch2.p), r)
# lay <- rbind(c(2,2,2, 2, NA),c(1,1,1,1,3),c(1,1,1,1,3), c(1,1,1,1,3), c(1,1,1,1,3))
# gaga.hd3<- grid.arrange(grobs = roch2.g, layout_matrix = lay)

#roch3

ggplot(roch3.c, aes(x=x, y=y, label=FINAlab, color=genet.1, shape=shape))+
  geom_point(size=2)+mytheme+labs(x="Distance (m)", y="Distance (m)")+
  geom_text_repel(aes(label=FINAlab),hjust=0, vjust=0, size=3)+
  scale_shape_manual(values=c("gen"=19, "ui"=15))+
  scale_color_manual(values=c("A2"="dodgerblue","ui"="black"))+coord_fixed(ratio = 1)

roch3.g.pre<-ggplot(roch3.c, aes(x=x, y=y, label=FINAlab, color=genet.1, shape=shape))+
  geom_point(size=1)+mytheme+labs(x="Distance (m)", y="Distance (m)")+
  geom_text_repel(aes(label=FINAlab),hjust=0, vjust=0, size=1.5, segment.size = .05, segment.alpha = 1)+
  scale_shape_manual(values=c("gen"=19, "ui"=15))+
  scale_color_manual(values=c("A2"="dodgerblue","ui"="black"))+ggtitle("Roch3 2007")

# roch3.p.pre<-upgma(dist(roch3))
# 
# roch3.p<-ggtree(roch3.p.pre, color="#9400D3")+geom_tiplab(hjust=.5, angle=-45, size=3)+scale_x_reverse()+coord_flip()+geom_tippoint(shape=18,size=3, col="black")
# r<-rectGrob(gp=gpar(fill="white", col="white"))
# roch3.g<-list(ggplotGrob(roch3.g.pre), ggplotGrob(roch3.p), r)
# lay <- rbind(c(2,2,2, 2, NA),c(1,1,1,1,3),c(1,1,1,1,3), c(1,1,1,1,3), c(1,1,1,1,3))
# gaga.hd3<- grid.arrange(grobs = roch3.g, layout_matrix = lay)


######### COMBINE INTO PANELS

plotlist.aflp2<-list(S.g.pre, ce02.g.pre, ce06.g.pre, jl.g.pre, str.g.pre, r, roch1.g.pre, roch2.g.pre, roch3.g.pre, hd1.g.pre, hd2.g.pre, hd3.g.pre, d1.g.pre, d2.g.pre, r, d3.g.pre, d4.g.pre, r)


plotlist.aflp1<-list(S.g.pre, ce02.g.pre, ce06.g.pre, jl.g.pre, str.g.pre, r, roch1.g.pre, roch2.g.pre, roch3.g.pre)
grid.arrange(grobs=plotlist.aflp1, ncol=3, nrow=3)

plotlist.aflp2<-list(hd1.g.pre, hd2.g.pre, hd3.g.pre, d1.g.pre, d2.g.pre, r, d3.g.pre, d4.g.pre, r)
grid.arrange(grobs=plotlist.aflp2, ncol=3, nrow=3)


plotlist.snp<-list(mira.g.pre, agr.g.pre, d2snp.g.pre, d3snp.g.pre, d2.14.g.pre, d3.14.g.pre, d2.15.g.pre, d3.15.g.pre)
grid.arrange(grobs=plotlist.snp, nrow=4, ncol=2)
