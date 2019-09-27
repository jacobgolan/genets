##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
########################################                    Fst                             ##############################################
##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################

setwd("~/Dropbox/Genet Size/Jacob Analyses")

hier <- read.genalex("AFLP_genalex_4hier.csv")

hier.gi<-as.genind(as.matrix(hier))
hier.gl<-as.genlight(as.matrix(hier))


library(hierfstat)
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



test<-as.data.frame(cbind(pops.aflp, hier.gi$tab))
test[, -1] <- sapply(test[, -1], as.numeric)
efestee<-genet.dist(test, diploid=T, method="WC84")
efestee.dist<-as.dist(efestee)
library(lattice)
heatmap(as.matrix(efestee.dist),symm=TRUE)

new.palette=colorRampPalette(c("black","red","yellow","white"),space="rgb") 
levelplot(as.matrix(efestee.dist)[1:ncol(as.matrix(efestee.dist)),ncol(as.matrix(efestee.dist)):1],col.regions=new.palette(20),scales=list(x=list(rot=90)), panel=as.matrix(efestee))

library(gplots)
heatmap.2(as.matrix(efestee.dist),           # cell labeling
          cellnote=round(as.matrix(efestee.dist), digits =3),
          notecex=.5, notecol="black", trace="none", key=F, cexRow = .75, cexCol = .75, dendrogram="row")

cleanFST<-as.matrix(efestee.dist)
cleanFST<-cleanFST[-c(1,2,13,19), -c(1,2,13,19)]
heatmap.2(cleanFST,           # cell labeling
          cellnote=round(cleanFST, digits =3),
          notecex=.5, notecol="black", trace="none", key=F, cexRow = .75, cexCol = .75, dendrogram="row")

# clone correct ********************** THIS ONE for AFLP
hier.gc<-as.genclone(hier.gi)

strata(hier.gc) <- as.data.frame(pop(hier.gi)) #set population as stratum
hier.gc$strata
hier.gc.cc<-clonecorrect(hier.gc)


plot(nj(dist(hier.gc.cc)))

pop.cc<-as.data.frame(hier.gc.cc$pop)

test2<-as.data.frame(cbind(pop.cc, hier.gc.cc$tab))
test2[, -1] <- sapply(test2[, -1], as.numeric)
colnames(test2)[1]<-"pops.aflp.cc"
efestee.cc<-genet.dist(test2, diploid=T, method="WC84")
efestee.cc.dist<-as.dist(efestee.cc)

cleanFST.cc<-as.matrix(efestee.cc.dist)
cleanFST.cc<-cleanFST.cc[-c(3,4,11,16,17), -c(3,4,11,16,17)] # *********  I had to remove rochester 3 because there are onlt two ind
heatmap.2(cleanFST.cc,           # cell labeling
          cellnote=round(cleanFST.cc, digits =3),
          notecex=.5, notecol="black", trace="none", key=F, cexRow = .75, cexCol = .75, dendrogram="row")

### Reorder matrices to follow Anne's configuration
# Serbia, CESAC, Rochester, JL, str, HD, Drake 1-4


cleanFST.cc.reorder<-cleanFST.cc[c(13,4,3,1,2,12,14,9,10,11,5,6,7,8),c(13,4,3,1,2,12,14,9,10,11,5,6,7,8)]

cleanFST.cc.reorder.upper<-cleanFST.cc.reorder
cleanFST.cc.reorder.upper[upper.tri(cleanFST.cc.reorder.upper, diag=T)] <- NA

heatmap.2(cleanFST.cc.reorder.upper,           # cell labeling
          cellnote=round(cleanFST.cc.reorder.upper, digits =3),
          notecex=1, notecol="black", trace="none", key=F, cexRow = 1.25, cexCol = 1.25, dendrogram="none", symm=F ,Rowv=FALSE, Colv=FALSE, srtCol=45)

cleanFST.cc.reorder.dist<-as.dist(cleanFST.cc.reorder)
#plot(hclust(cleanFST.cc.reorder.dist))

hc <- upgma((cleanFST.cc.reorder.dist))
plot(hc)


popsub(hier.gc.cc, sublist = "Rochester_NY_2007_1")
popsub(hier.gc.cc, sublist = "Rochester_NY_2007_2")
popsub(hier.gc.cc, sublist = "Rochester_NY_2007_3")

popsub(hier.gc, sublist = "Rochester_NY_2007_1")
popsub(hier.gc, sublist = "Rochester_NY_2007_2")
popsub(hier.gc, sublist = "Rochester_NY_2007_3")


#Drake
efestee.mat<-as.matrix(efestee)
efestee.drakes<-efestee.mat[c(5,6,7,8),c(5,6,7,8)]
levelplot(efestee.drakes[1:ncol(efestee.drakes),ncol(efestee.drakes):1],col.regions=new.palette(20),scales=list(x=list(rot=90)))

#HD
efestee.HD<-efestee.mat[c(9,10,11),c(9,10,11)]
levelplot(efestee.HD[1:ncol(efestee.HD),ncol(efestee.HD):1],col.regions=new.palette(20),scales=list(x=list(rot=90)))

#east
efestee.east<-efestee.mat[c(12,14,15,16,17),c(12,14,15,16,17)]
levelplot(efestee.east[1:ncol(efestee.east),ncol(efestee.east):1],col.regions=new.palette(20),scales=list(x=list(rot=90)))


#EUROPE
efestee.EU<-efestee.mat[c(3,4,18),c(3,4,18)]
levelplot(efestee.EU[1:ncol(efestee.EU),ncol(efestee.EU):1],col.regions=new.palette(20),scales=list(x=list(rot=90)))



write.csv(as.matrix(efestee.dist), "AFLP_pairwiseFst.csv")


View(test)
test.sub.compare2SNP.drake<-test[c(106:117,120:124),]
efestee.d2.d3comSNP<-genet.dist(test.sub.compare2SNP.drake, diploid=T, method="WC84") #0.01875621

############

#SNP
setwd("~/Dropbox/finalizing_genets")

fst.mx.5.time<-read.csv("Fst.csv", head=T, stringsAsFactors = F)
rownames(fst.mx.5.time)<-fst.mx.5.time[,1]
fst.mx.5.time<-as.matrix(fst.mx.5.time[,-1])

gplots::heatmap.2(fst.mx.5.time,           # cell labeling
          cellnote=round(fst.mx.5.time, digits =3),
          notecex=1, notecol="black", trace="none", key=F, cexRow = .5, cexCol = .5, dendrogram="none", symm=F ,Rowv=FALSE, Colv=FALSE, srtCol=45)


fst.mx.5.time.mlt<-na.omit(reshape2::melt(fst.mx.5.time))
fst.mx.5.time.mlt$pair<-paste(fst.mx.5.time.mlt$Var1, fst.mx.5.time.mlt$Var2,sep="&")

fst.mx.5.time.mlt$pair<-plyr::revalue(fst.mx.5.time.mlt$pair, c(
  'agraria_2015&mira_2015'='PT',
  'agraria_2015&vilarinho_2015'='PT',
  'mira_2015&vilarinho_2015'='PT',
  'agraria_2015&Drake2_2004'='X_yr',
  'mira_2015&Drake2_2004'='X_yr',
  'vilarinho_2015&Drake2_2004'='X_yr',
  'agraria_2015&Drake2_2014'='X_yr',
  'mira_2015&Drake2_2014'='X_yr',
  'vilarinho_2015&Drake2_2014'='X_yr',
  'Drake2_2004&Drake2_2014'='US_yr',
  'agraria_2015&Drake2_2015'='X_yr',
  'mira_2015&Drake2_2015'='X',
  'vilarinho_2015&Drake2_2015'='X',
  'Drake2_2004&Drake2_2015'='US_yr',
  'Drake2_2014&Drake2_2015'='US_yr',
  'agraria_2015&Drake3_2004'='X_yr',
  'mira_2015&Drake3_2004'='X_yr',
  'vilarinho_2015&Drake3_2004'='X_yr',
  'Drake2_2004&Drake3_2004'='USX',
  'Drake2_2014&Drake3_2004'='USX_yr',
  'Drake2_2015&Drake3_2004'='USX_yr',
  'agraria_2015&Drake3_2014'='X_yr',
  'mira_2015&Drake3_2014'='X_yr',
  'vilarinho_2015&Drake3_2014'='X_yr',
  'Drake2_2004&Drake3_2014'='USX_yr',
  'Drake2_2014&Drake3_2014'='USX',
  'Drake2_2015&Drake3_2014'='USX_yr',
  'Drake3_2004&Drake3_2014'='US_yr',
  'agraria_2015&Drake3_2015'='X',
  'mira_2015&Drake3_2015'='X',
  'vilarinho_2015&Drake3_2015'='X',
  'Drake2_2004&Drake3_2015'='USX_yr',
  'Drake2_2014&Drake3_2015'='USX_yr',
  'Drake2_2015&Drake3_2015'='USX',
  'Drake3_2004&Drake3_2015'='US_yr',
  'Drake3_2014&Drake3_2015'='US_yr'

))


library(magrittr)
fst.mx.5.time.mlt[,c("value", "pair")] %>% 
  dplyr::group_by(pair) %>% 
  dplyr::summarise(
    n = n(),
    mean=mean(value),
    median=median(value),
    sd=sd(value),
    se=sd/sqrt(n),
    max=max(value),
    min=min(value)) 















fst.mx.5.subpops<-read.csv("Fst_w_subpops.csv", head=T, stringsAsFactors = F)
rownames(fst.mx.5.subpops)<-fst.mx.5.subpops[,1]
fst.mx.5.subpops<-as.matrix(fst.mx.5.subpops[,-1])
fst.mx.5.subpops<-fst.mx.5.subpops[-7,-7]

fst.mx.5.subpops.mlt<-na.omit(reshape2::melt(fst.mx.5.subpops))
fst.mx.5.subpops.mlt$pair<-paste(fst.mx.5.subpops.mlt$Var1, fst.mx.5.subpops.mlt$Var2,sep="&")
#View(unique(fst.mx.5.subpops.mlt$pair))

fst.mx.5.subpops.mlt$pair<-plyr::revalue(fst.mx.5.subpops.mlt$pair, c(
'agraria_2015&mira_2015'='PT',
'agraria_2015&vilarinho_2015'='PT',
'mira_2015&vilarinho_2015'='PT',
'agraria_2015&Drake2_2004'='X_yr',
'mira_2015&Drake2_2004'='X_yr',
'vilarinho_2015&Drake2_2004'='X_yr',
'agraria_2015&Drake2_2014'='X_yr',
'mira_2015&Drake2_2014'='X_yr',
'vilarinho_2015&Drake2_2014'='X_yr',
'Drake2_2004&Drake2_2014'='CUS_yr',
'agraria_2015&Drake2_2015'='X',
'mira_2015&Drake2_2015'='X',
'vilarinho_2015&Drake2_2015'='X',
'Drake2_2004&Drake2_2015'='CUS_yr',
'Drake2_2014&Drake2_2015'='CUS_yr',
'agraria_2015&Drake2s_2004'='X_yr',
'mira_2015&Drake2s_2004'='X_yr',
'vilarinho_2015&Drake2s_2004'='X_yr',
'Drake2_2004&Drake2s_2004'='CUSX',
'Drake2_2014&Drake2s_2004'='CUSX_yr',
'Drake2_2015&Drake2s_2004'='CUSX_yr',
'agraria_2015&Drake2s_2014'='X_yr',
'mira_2015&Drake2s_2014'='X_yr',
'vilarinho_2015&Drake2s_2014'='X_yr',
'Drake2_2004&Drake2s_2014'='CUSX_yr',
'Drake2_2014&Drake2s_2014'='CUSX',
'Drake2_2015&Drake2s_2014'='CUSX_yr',
'Drake2s_2004&Drake2s_2014'='CUS_yr',
'agraria_2015&Drake2s_2015'='X',
'mira_2015&Drake2s_2015'='X',
'vilarinho_2015&Drake2s_2015'='X',
'Drake2_2004&Drake2s_2015'='CUSX_yr',
'Drake2_2014&Drake2s_2015'='CUSX_yr',
'Drake2_2015&Drake2s_2015'='CUSX',
'Drake2s_2004&Drake2s_2015'='CUS_yr',
'Drake2s_2014&Drake2s_2015'='CUS_yr',
'agraria_2015&Drake3_2004'='X_yr',
'mira_2015&Drake3_2004'='X_yr',
'vilarinho_2015&Drake3_2004'='X_yr',
'Drake2_2004&Drake3_2004'='CUSX',
'Drake2_2014&Drake3_2004'='CUSX_yr',
'Drake2_2015&Drake3_2004'='CUSX_yr',
'Drake2s_2004&Drake3_2004'='CUSX',
'Drake2s_2014&Drake3_2004'='CUSX_yr',
'Drake2s_2015&Drake3_2004'='CUSX_yr',
'agraria_2015&Drake3_2014'='X_yr',
'mira_2015&Drake3_2014'='X_yr',
'vilarinho_2015&Drake3_2014'='X_yr',
'Drake2_2004&Drake3_2014'='CUSX_yr',
'Drake2_2014&Drake3_2014'='CUSX',
'Drake2_2015&Drake3_2014'='CUSX_yr',
'Drake2s_2004&Drake3_2014'='CUSX_yr',
'Drake2s_2014&Drake3_2014'='CUSX',
'Drake2s_2015&Drake3_2014'='CUSX_yr',
'Drake3_2004&Drake3_2014'='CUS_yr',
'agraria_2015&Drake3_2015'='X',
'mira_2015&Drake3_2015'='X',
'vilarinho_2015&Drake3_2015'='X',
'Drake2_2004&Drake3_2015'='CUSX_yr',
'Drake2_2014&Drake3_2015'='CUSX_yr',
'Drake2_2015&Drake3_2015'='CUSX',
'Drake2s_2004&Drake3_2015'='CUSX_yr',
'Drake2s_2014&Drake3_2015'='CUSX_yr',
'Drake2s_2015&Drake3_2015'='CUSX',
'Drake3_2004&Drake3_2015'='CUS_yr',
'Drake3_2014&Drake3_2015'='CUS_yr'))



fst.mx.5.subpops.mlt$subpair<-paste(fst.mx.5.subpops.mlt$Var1, fst.mx.5.subpops.mlt$Var2,sep="&")
fst.mx.5.subpops.mlt$subpair<-plyr::revalue(fst.mx.5.subpops.mlt$subpair, c(
  
'agraria_2015&mira_2015'='PT',
'agraria_2015&vilarinho_2015'='PT',
'mira_2015&vilarinho_2015'='PT',
'agraria_2015&Drake2_2004'='X_yr',
'mira_2015&Drake2_2004'='X_yr',
'vilarinho_2015&Drake2_2004'='X_yr',
'agraria_2015&Drake2_2014'='X_yr',
'mira_2015&Drake2_2014'='X_yr',
'vilarinho_2015&Drake2_2014'='X_yr',
'Drake2_2004&Drake2_2014'='CUS_yr',
'agraria_2015&Drake2_2015'='X',
'mira_2015&Drake2_2015'='X',
'vilarinho_2015&Drake2_2015'='X',
'Drake2_2004&Drake2_2015'='CUS_yr',
'Drake2_2014&Drake2_2015'='CUS_yr',
'agraria_2015&Drake2s_2004'='X_yr',
'mira_2015&Drake2s_2004'='X_yr',
'vilarinho_2015&Drake2s_2004'='X_yr',
'Drake2_2004&Drake2s_2004'='CUSX',
'Drake2_2014&Drake2s_2004'='cc_yr',
'Drake2_2015&Drake2s_2004'='cc_yr',
'agraria_2015&Drake2s_2014'='X_yr',
'mira_2015&Drake2s_2014'='X_yr',
'vilarinho_2015&Drake2s_2014'='X_yr',
'Drake2_2004&Drake2s_2014'='cc_yr',
'Drake2_2014&Drake2s_2014'='cc',
'Drake2_2015&Drake2s_2014'='cc_yr',
'Drake2s_2004&Drake2s_2014'='CUS_yr',
'agraria_2015&Drake2s_2015'='X',
'mira_2015&Drake2s_2015'='X',
'vilarinho_2015&Drake2s_2015'='X',
'Drake2_2004&Drake2s_2015'='cc_yr',
'Drake2_2014&Drake2s_2015'='cc_yr',
'Drake2_2015&Drake2s_2015'='cc',
'Drake2s_2004&Drake2s_2015'='CUS_yr',
'Drake2s_2014&Drake2s_2015'='CUS_yr',
'agraria_2015&Drake3_2004'='X_yr',
'mira_2015&Drake3_2004'='X_yr',
'vilarinho_2015&Drake3_2004'='X_yr',
'Drake2_2004&Drake3_2004'='CUSX',
'Drake2_2014&Drake3_2004'='CUSX_yr',
'Drake2_2015&Drake3_2004'='CUSX_yr',
'Drake2s_2004&Drake3_2004'='CUSX',
'Drake2s_2014&Drake3_2004'='CUSX_yr',
'Drake2s_2015&Drake3_2004'='CUSX_yr',
'agraria_2015&Drake3_2014'='X_yr',
'mira_2015&Drake3_2014'='X_yr',
'vilarinho_2015&Drake3_2014'='X_yr',
'Drake2_2004&Drake3_2014'='CUSX_yr',
'Drake2_2014&Drake3_2014'='CUSX',
'Drake2_2015&Drake3_2014'='CUSX_yr',
'Drake2s_2004&Drake3_2014'='CUSX_yr',
'Drake2s_2014&Drake3_2014'='CUSX',
'Drake2s_2015&Drake3_2014'='CUSX_yr',
'Drake3_2004&Drake3_2014'='CUS_yr',
'agraria_2015&Drake3_2015'='X',
'mira_2015&Drake3_2015'='X',
'vilarinho_2015&Drake3_2015'='X',
'Drake2_2004&Drake3_2015'='CUSX_yr',
'Drake2_2014&Drake3_2015'='CUSX_yr',
'Drake2_2015&Drake3_2015'='CUSX',
'Drake2s_2004&Drake3_2015'='CUSX_yr',
'Drake2s_2014&Drake3_2015'='CUSX_yr',
'Drake2s_2015&Drake3_2015'='CUSX',
'Drake3_2004&Drake3_2015'='CUS_yr',
'Drake3_2014&Drake3_2015'='CUS_yr'))

fst.mx.5.subpops.mlt[,c("value", "pair")] %>% 
  dplyr::group_by(pair) %>% 
  dplyr::summarise(
    n = n(),
    mean=mean(value),
    median=median(value),
    sd=sd(value),
    se=sd/sqrt(n),
    max=max(value),
    min=min(value)) 

fst.mx.5.subpops.mlt[,c("value", "pair", "subpair")] %>% 
  dplyr::group_by(subpair) %>% 
  dplyr::summarise(
    n = n(),
    mean=mean(value),
    median=median(value),
    sd=sd(value),
    se=sd/sqrt(n),
    max=max(value),
    min=min(value)) 


gplots::heatmap.2(fst.mx.5.subpops,           # cell labeling
                  cellnote=round(fst.mx.5.subpops, digits =3),
                  notecex=1, notecol="black", trace="none", key=F, cexRow = .5, cexCol = .5, dendrogram="none", symm=F ,Rowv=FALSE, Colv=FALSE, srtCol=45)


