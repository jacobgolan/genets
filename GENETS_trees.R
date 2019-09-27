
################################ ################################ ################################ 
################################ ################################ ################################ 
################################ ################################ ################################ 
################################ ################################ ################################ 

# Trees

setwd("~/Dropbox/Genet Size/Jacob Analyses")

library(phangorn)
library(hierfstat)
library(poppr)
library(adegenet)
library("RColorBrewer")
library(gplots)
library(vcfR)
library(magrittr)
library(ggtree)


#AFLP
hier <- read.genalex("AFLP_genalex_AMANITABASED.csv")

hier.gi<-as.genind(as.matrix(hier))
#saveRDS(hier.gi,"hier.gi.rds")
#hier.gl<-as.genlight(as.matrix(hier))
aflp.tree<-upgma(dist(hier.gi))
aflp.tree.nj<-nj(dist(hier.gi))

#aflp.tree.boot<-aboot(hier.gl, dist = bitwise.dist, sample = 200, tree = "upgma", cutoff = 0, quiet = TRUE)


tips.aflp<-as.data.frame(aflp.tree$tip.label)
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

pops.aflp.df<-as.data.frame(pops.aflp)
meta.aflp<-cbind(tips.aflp, pops.aflp.df)


dd <- meta.aflp
row.names(dd) <- NULL
colnames(dd)<-c("tip", "pop")
print(dd)


cols <- rainbow(unique(dd$pop))
colman<-c("#ff0000", "#0039e6", "#00b300",  "#ff6600", "#006699",
          "#7300e6", "#cc00cc", "#00b3b3", "#e6e600", "#666666",
          "#66cc00","#800040", "#669900", "#b366ff", "#000080", "#669999", "#66ffb3", "#ff3399", "#00ff00")


plot.phylo(aflp.tree, cex = 0.5, font = 2, adj = 0, tip.color = colman[dd$pop],
           label.offset = 0.001, type="fan", align.tip.label = T)
nodelabels(aflp.tree$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.8,
           font = 3, xpd = TRUE)
axisPhylo(3)


cols[dd$pop]



p <- ggtree(aflp.tree)
c <- ggtree(aflp.tree, layout="circular")


p <- p %<+% dd + geom_tippoint(aes(color=pop))
p+theme(legend.position="right") + geom_tiplab(size=2)+ggtitle("AFLP Mushrooms: Euclidean Distance")


p <- p %<+% dd + geom_tippoint(aes(color=pop), size=.4)
p+theme(legend.position="right") + geom_tiplab(size=1.25, offset = .002)+ggtitle("AFLP Mushrooms: Euclidean Distance")

c <- c %<+% dd + geom_tippoint(aes(color=pop))
c+theme(legend.position="right") + geom_tiplab2(aes(angle=angle), offset(-.2), size=2)+ggtitle("SNP Mushrooms: Euclidean Distance")


c+theme(legend.position="right") + geom_tiplab2(aes(angle=angle), offset(-.2), size=2)+ggtitle("SNP Mushrooms: Euclidean Distance")


c <- ggtree(aflp.tree, layout="circular")

myBoots <- boot.phylo(aflp.tree, hier.gi$tab, 1000, FUN = function(xx) upgma(dist(xx)))
plot(aflp.tree, show.tip=FALSE, edge.width=2)
nodelabels(myBoots/1000, cex=.6)


aflp.cores=values=c("#6600cc", "#d966ff", "#9900cc", "#993366", "#ff6699", "#ff0000", "#ff9900", "#cc5200", "#ffff00", "#999900","#e6b800","#0000ff",
                    "#ffcccc", "#8080ff", "#0099ff", "#006bb3", "#b3e0ff", "#df9fbf", "#4d0066")

c <- ggtree(aflp.tree, layout="circular")
c <- c %<+% dd + geom_tippoint(aes(color=pop))
c+theme(legend.position="right") + geom_tiplab2(aes(angle=angle), offset(-.2), size=2.5)+
  scale_color_manual(values=c("#6600cc", "#d966ff", "#9900cc", "#993366", "#ff6699", "#ff0000", "#ff9900", "#cc5200", "#ffff00", "#999900","#e6b800","#0000ff",
                              "#ffcccc", "#8080ff", "#0099ff", "#006bb3", "#b3e0ff", "#df9fbf", "#4d0066"))+
  ggtitle("AFLP Mushrooms: Euclidean Distance")+geom_treescale()+geom_nodelab(label=((myBoots/1000)*100), size=2)


################################################################################################
####### Just D2 + D3, 2004
hier.D.2004<-popsub(hier.gi, sublist=c("Drake_2_2004","Drake_3_2004"))
D.2004.tree<-upgma(dist(hier.D.2004))

myBoots.D.2004 <- boot.phylo(D.2004.tree, hier.D.2004$tab, 1000, FUN = function(xx) upgma(dist(xx)))
plot(D.2004.tree, show.tip=FALSE, edge.width=2)
nodelabels((myBoots.D.2004/1000)*100, cex=.6)
tiplabels(D.2004.tree$tip.label)

####### Just Invasive
hier.g.invi<-hier.gi


#aflp.tree.boot<-aboot(hier.gl, dist = bitwise.dist, sample = 200, tree = "upgma", cutoff = 50, quiet = TRUE)


pops.aflp.inv<-c(
  replicate(13, "Rochester_NY_2007_1"),
  replicate(7, "Rochester_NY_2007_2"),
  replicate(2, "Rochester_NY_2007_3"),
  replicate(19, "Drake_1_2004"),
  replicate(13, "Drake_2_2004"),
  replicate(8, "Drake_3_2004"),
  replicate(6, "Drake_4_2004"),
  replicate(7, "Hearts_Desire_1_2004"),
  replicate(7, "Hearts_Desire_2_2004"),
  replicate(4, "Hearts_Desire_3_2004"),
  replicate(44, "Jakes_Landing_2006"),
  replicate(1, "Monterrey_Pisto_2006"),
  replicate(11, "New_Jersey_str_2006")
)

pop(hier.g.invi)<-pops.aflp.inv

aflp.invas.gi<-popsub(hier.g.invi, sublist = c("Rochester_NY_2007_1","Rochester_NY_2007_2","Rochester_NY_2007_3",
                                               "Drake_1_2004","Drake_2_2004","Drake_3_2004","Drake_4_2004",
                                               "Hearts_Desire_1_2004","Hearts_Desire_2_2004","Hearts_Desire_3_2004",
                                               "Jakes_Landing_2006","New_Jersey_str_2006",
                                               "Monterrey_Pisto_2006"))

aflp.tree.inv<-upgma(dist(aflp.invas.gi))
tips.inv.aflp<-as.data.frame(aflp.tree.inv$tip.label)


pops.aflp.inv.df<-as.data.frame(pops.aflp.inv)
meta.aflp.inv<-cbind(tips.inv.aflp, pops.aflp.inv.df)



p <- ggtree(aflp.tree.inv)

dd <- meta.aflp.inv
row.names(dd) <- NULL
colnames(dd)<-c("tip", "pop")
print(dd)

cols <- rainbow(unique(dd$pop))

plot.phylo(aflp.tree.boot, cex = 0.4, font = 2, adj = 0, tip.color = cols[dd$pop],
           label.offset = 0.000125)
nodelabels(aflp.tree.boot$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.8,
           font = 3, xpd = TRUE)
axisPhylo(3)


p <- p %<+% dd + geom_tippoint(aes(color=pop))
p+theme(legend.position="right") + geom_tiplab( size=2)+ggtitle("AFLP Mushrooms: Euclid Distance")

p <- p %<+% dd + geom_tippoint(aes(color=pop), size=.4)
p+theme(legend.position="right") + geom_tiplab(size=1.25, offset = .002)+ggtitle("AFLP Mushrooms: Nei's Distance")

######## Just CA


pops.aflp.CA<-c(
  replicate(19, "Drake_1_2004"),
  replicate(13, "Drake_2_2004"),
  replicate(8, "Drake_3_2004"),
  replicate(6, "Drake_4_2004"),
  replicate(7, "Hearts_Desire_1_2004"),
  replicate(7, "Hearts_Desire_2_2004"),
  replicate(4, "Hearts_Desire_3_2004")
  
)

#pop(hier.gi)<-pops.aflp

hier.g.CA<-hier.gi


aflp.CA.gi<-popsub(hier.gi, sublist = c("Drake_1_2004","Drake_2_2004","Drake_3_2004","Drake_4_2004",
                                        "Hearts_Desire_1_2004","Hearts_Desire_2_2004","Hearts_Desire_3_2004"))

aflp.tree.CA<-nj(dist(aflp.CA.gi))
tips.CA.aflp<-as.data.frame(aflp.tree.CA$tip.label)


pops.aflp.CA.df<-as.data.frame(pops.aflp.CA)
meta.aflp.CA<-cbind(tips.CA.aflp, pops.aflp.CA.df)



p <- ggtree(aflp.tree.CA)

dd <- meta.aflp.CA
row.names(dd) <- NULL
colnames(dd)<-c("tip", "pop")
print(dd)

cols <- rainbow(unique(dd$pop))

plot.phylo(aflp.tree.boot, cex = 0.4, font = 2, adj = 0, tip.color = cols[dd$pop],
           label.offset = 0.000125)
nodelabels(aflp.tree.boot$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.8,
           font = 3, xpd = TRUE)
axisPhylo(3)


p <- p %<+% dd + geom_tippoint(aes(color=pop))
p+theme(legend.position="right") + geom_tiplab(size=1.25, offset = .002)+ggtitle("AFLP Mushrooms: Euclid Distance")







#SNPS

setwd("~/Dropbox/finalizing_genets")

maf006yall<-read.PLINK("HESS_fix/mx0.5/vcf.MMDP60.maf0.006.minQ30.recode.raw", 
                            map.file = "HESS_fix/mx0.5/vcf.MMDP60.maf0.006.minQ30.recode.map", quiet = FALSE,
                            parallel = require("parallel"), n.cores = NULL)


snp.tree<-upgma(dist((maf006yall)))
#snp.tree.boot<-aboot(maf006, dist = bitwise.dist, sample = 200, tree = "upgma", cutoff = 50, quiet = TRUE)


tips.snp<-as.data.frame(snp.tree$tip.label)

myPops <- as.factor(c('old_EU',	'old_EU',	'old_EU',	'old_EU',	'old_EU',	'old_EU',	'CA_old',	'old_EU',	'D2.04',	'D2.04',	
                      'D2.04',	'D2.04',	'D2.04',	'D2.04',	'D2.04',	'D2.04',	'D2.04',	'D2.04',	'D2.04',	'D2.04',
                      'D2.04',	'D3.04',	'D3.04',	'D3.04',	'D3.04',	'D3.04',	'old_EU',	'D2.14',	'D2.14',	'D2.14',
                      'D2.14',	'D2.14',	'D2.14',	'D2.14',	'D2.14',	'D2.14',	'D2.14',	'D2.14',	'D2.14',	'D2.14',
                      'D2.14',	'D2.14',	'D2.14',	'D2.14',	'D2.14',	'D2.14',	'D2.14',	'D2.14',	'D2.14',	'D2.14',
                      'D2.14',	'D2.14',	'D3.14',	'D3.14',	'D3.14',	'D3.14',	'D3.14',	'D3.14',	'D3.14',	'D3.14',
                      'D3.14',	'Agraria',	'Agraria',	'Vilarinho',	'Vilarinho',	'Vilarinho',	'Vilarinho',	'Vilarinho',
                      'Mira',	'Mira',	'Mira',	'Mira',	'D2.15',	'D2.15',	'D2.15',	'D2.15',	'D2.15',	'D2.15',	'D2.15',
                      'D2.15',	'D2.15',	'D2.15',	'D2.15',	'D3.15',	'D3.15',	'D3.15'))

pop(maf006yall)<-myPops


#two closest indvs
MIRA<-popsub(maf006yall, sublist=c("Mira"))
MIRAcl<-MIRA[2:3,]
MIRAcl<-as.genind(as.matrix(MIRAcl))
MIRAnm<-missingno(MIRAcl, type="loci", quiet=T) #takes like 20 minutes!
dist(MIRAnm, method="manhattan") #24157
dist(MIRAcl, method="manhattan") #24250.29

dist(MIRAcl, method="manhattan")/dist(MIRAnm, method="manhattan")


head(as.genind(as.matrix(MIRAnm))$tab[,1:10])
dist(gl.filter.monomorphs(MIRAnm), method = "manhattan")
informloci(as.genind(as.matrix(MIRAnm)))

bitwise.dist(MIRAnm, euclidean = T)
dist(MIRAnm)

na.index <- is.na(MIRAcl$tab[1,]) | is.na(MIRAcl$tab[2,])
dist(rbind(MIRAcl$tab[1,][!na.index], MIRAcl$tab[1,][!na.index])) * sqrt(length(MIRAcl$tab[1,]) / length(MIRAcl$tab[1,][!na.index]))

dist(rbind(c(0,1,2.5),c(0,1,1), method="m"))

#View(as.matrix(dist(maf006yall, method="manhattan")))


D3.04<-popsub(maf006yall, sublist=c("D3.04"))
mat<-as.matrix(dist(D3.04, method="manhattan"))
mat[lower.tri(mat, diag=T)] <- NA
matm<-melt(mat)
ja<-(matm[complete.cases(matm), ])
mean(ja$value)

pops.snp.df<-as.data.frame(myPops)
meta.snp<-cbind(tips.snp, pops.snp.df)


library(ggtree)
c <- ggtree(snp.tree, layout = "circular")

dd <- meta.snp
row.names(dd) <- NULL
colnames(dd)<-c("tip", "pop")
print(dd)

cols <- rainbow(unique(dd$pop))

plot.phylo(snp.tree, cex = 0.2, font = 2, adj = 0, tip.color = cols[dd$pop],
           label.offset = 0.0125)
nodelabels(snp.tree$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.8,
           font = 3, xpd = TRUE)
axisPhylo(3)


c <- c %<+% dd + geom_tippoint(aes(color=pop))
c+theme(legend.position="right") + geom_tiplab2(aes(angle=angle), offset(-.2), size=2)+ggtitle("SNP Mushrooms: Euclidean Distance")

c <- ggtree(snp.tree, layout = "circular")

myBoots.snp <- boot.phylo(treeUPGMA.snp, maf006.mx.5, 500, FUN = function(xx) upgma(dist(xx)))
plot(treeUPGMA.snp, show.tip=FALSE, edge.width=.5)
nodelabels(myBoots.snp, cex=.6)

c <- ggtree(snp.tree, layout="circular")
c <- c %<+% dd + geom_tippoint(aes(color=pop))
c+theme(legend.position="right") + geom_tiplab2(aes(angle=angle), offset(-.2), size=3.5)+
  scale_color_manual(values=c("#cc33ff", "#663300", "#ff0000", "#b30000", "#ff6666", "#ff9900", "#cc7a00", "#cc9900", "#9966ff", "#660066", "#6600ff"))+
  ggtitle("SNP Mushrooms: Euclidean Distance")+geom_treescale()+geom_nodelab(label=(myBoots.snp/500)*100, size=2)

myBoots.snp <- boot.phylo(treeUPGMA.snp, maf006.mx.5, 100, FUN = function(xx) upgma(dist(xx)))
plot(treeUPGMA.snp, show.tip=FALSE, edge.width=.5)
nodelabels((myBoots.snp/500)*100, cex=.6, frame="none")


##################################
# Just Europe

EU.aflp<-rownames(hier$tab[c(23:85,132,195,197:210),])

hier.EU<-hier[i=EU.aflp]
aflp.tree.EU<-upgma(dist(hier.EU))


tips.aflp.EU<-as.data.frame(aflp.tree.EU$tip.label)
pops.aflp.EU<-c(
  replicate(1, "Austria"),
  replicate(37, "CESAC_2006"),
  replicate(25, "CESAC_2002"),
  replicate(1, "Spain_2006"),
  replicate(1, "Arfons_2007"),
  replicate(14, "Serbia_2007")
)

pops.aflp.EU.df<-as.data.frame(pops.aflp.EU)
meta.aflp.EU<-cbind(tips.aflp.EU, pops.aflp.EU.df)

p <- ggtree(aflp.tree.EU, layout = "circular")
dd <- meta.aflp.EU
row.names(dd) <- NULL
colnames(dd)<-c("tip", "pop")
print(dd)
p <- p %<+% dd + geom_tippoint(aes(color=pop))
p+theme(legend.position="right") + geom_tiplab2(aes(angle=angle), offset(-.2), size=2)+ggtitle("AFLP Europe Mushrooms: Euclidean Distance")

q <- ggtree(aflp.tree.EU)
dd <- meta.aflp.EU
row.names(dd) <- NULL
colnames(dd)<-c("tip", "pop")
print(dd)
q <- q %<+% dd + geom_tippoint(aes(color=pop))
q+theme(legend.position="right") + geom_tiplab(size=2)+ggtitle("AFLP Europe Mushrooms: Euclidean Distance")



maf.012
maf.012.EU.ind<-as.character(maf.012$ind.names[c(1:6,8,27,62:72)])

maf.012.EU<-maf.012[i=as.factor(maf.012.EU.ind)]

maf.012.EU.tree<-upgma(dist(maf.012.EU))
tips.snp.EU<-as.data.frame(maf.012.EU.tree$tip.label)


myPops.EU <- as.factor(c('old_EU',	'old_EU',	'old_EU',	'old_EU',	'old_EU',	'old_EU',	'old_EU',	'old_EU',	
                         'Agraria',	'Agraria',	'Vilarinho',	'Vilarinho',	'Vilarinho',	'Vilarinho',	'Vilarinho',
                         'Mira',	'Mira',	'Mira',	'Mira'))




pops.snp.df.EU<-as.data.frame(myPops.EU)
meta.snp.EU<-cbind(tips.snp.EU, pops.snp.df.EU)

p <- ggtree(maf.012.EU.tree, layout = "circular")
dd <- meta.snp.EU
row.names(dd) <- NULL
colnames(dd)<-c("tip", "pop")
print(dd)
p <- p %<+% dd + geom_tippoint(aes(color=pop))
p+theme(legend.position="right") + geom_tiplab2(aes(angle=angle), offset(-.2), size=2)+ggtitle("AFLP Europe Mushrooms: Euclidean Distance")

q <- ggtree(maf.012.EU.tree)
dd <- meta.snp.EU
row.names(dd) <- NULL
colnames(dd)<-c("tip", "pop")
print(dd)
q <- q %<+% dd + geom_tippoint(aes(color=pop))
q+theme(legend.position="right") + geom_tiplab(size=3)+ggtitle("SNP Europe Mushrooms: Euclidean Distance")


###############################################################################
####### Just D2 + D3, 2004
snp.D.2004<-popsub(maf006yall, sublist=c("D2.04","D3.04"))
snp.D.2004.tree<-upgma(dist(snp.D.2004))

myBoots.snp.D.2004 <- boot.phylo(snp.D.2004.tree, snp.D.2004$tab, 1000, FUN = function(xx) upgma(dist(xx)))
plot(snp.D.2004.tree, show.tip=FALSE, edge.width=2)
nodelabels((myBoots.snp.D.2004/1000)*100, cex=.6)
tiplabels(snp.D.2004.tree$tip.label)



################################ ################################ ################################ 
################################ ################################ ################################ 
################################ ################################ ################################ 
###### # Comapre number of SNPs among individuals

#df<-maxmiss.NOMISS.maf.049.recode.FIN.vcf.gz.rds$tab
#setdiff((as.character(df[14,])), (as.character(df[14,])))

#df<-df[c(1:5),]

#ones<-vector()
#twos<-vector()
#zeros<-vector()
#ind1<-vector()
#ind2<-vector()

#z<-1
#x<-1
#while(x<nrow(df)+1){
#  y<-x+1
#  while(y<nrow(df)+1){
#    ones[z]<-sum(colSums(df[c(x:y),]==1,na.rm = T)==2)
#    twos[z]<-sum(colSums(df[c(x:y),]==2,na.rm = T)==2)
#    zeros[z]<-sum(colSums(df[c(x:y),]==0,na.rm = T)==2)
#    ind1[z]<-row.names(df)[x]
#    ind2[z]<-row.names(df)[y]
#    y<-y+1
#    z<-z+1
#  }

#  x<-x+1
#}



#out.df<-as.data.frame(cbind(ind1,ind2,ones,twos,zeros),stringsAsFactors = F)

#count number of 0, 1, and 2 in each individual

setwd("~/Dropbox/Genet Size/Jacob Analyses")
maf006yall<-readRDS("maf006.gi.rds")
maf006yall.tab<-maf006yall$tab

maf0.006.df<-as.data.frame(as.matrix(maf006yall.tab))
zeros<-rowSums(maf0.006.df == 0, na.rm=T)
ones<-rowSums(maf0.006.df == 1, na.rm=T)
twos<-rowSums(maf0.006.df == 2, na.rm=T)
nadas<-rowSums(is.na(maf0.006.df))

maf0.006.alleles<-cbind(zeros, ones, twos, nadas)
#write.csv(maf0.006.alleles, "maf0.006.all.012.csv")

###
#maf0.012.no.df<-as.data.frame(as.matrix(maf0.012.no))
#zeros<-rowSums(maf0.012.no.df == 0, na.rm=T)
#ones<-rowSums(maf0.012.no.df == 1, na.rm=T)
#twos<-rowSums(maf0.012.no.df == 2, na.rm=T)
#nadas<-rowSums(is.na(maf0.012.no.df))

#maf0.012.no.alleles<-cbind(zeros, ones, twos, nadas)
#write.csv(maf0.012.no.alleles, "maf0.012.no.012.csv")

###
#maf0.49.nomiss.df<-as.data.frame(as.matrix(maf0.49.maxmiss0.9new.all))
#zeros<-rowSums(maf0.49.nomiss.df == 0, na.rm=T)
#ones<-rowSums(maf0.49.nomiss.df == 1, na.rm=T)
#twos<-rowSums(maf0.49.nomiss.df == 2, na.rm=T)
#nadas<-rowSums(is.na(maf0.49.nomiss.df))

#maf0.49.alleles<-cbind(zeros, ones, twos, nadas)
#write.csv(maf0.49.alleles, "maf0.049.all.012.csv")

###########################################################
###########################################################
###########################################################
####### Hwe
setwd("~/Dropbox/finalizing_genets")
d2.14.hwetest<-read.table("D2.14.hwe.txt", head=T)
d2.14.hwetest.cl<-d2.14.hwetest[d2.14.hwetest$E.HOM1.HET.HOM2 !=   "nan/nan/nan",]

ggplot(d2.14.hwetest.cl, aes(x=P_HWE))+geom_density()

vil.hwetest<-read.table("vilarinho.hwe.txt", head=T)
vil.hwetest.cl<-vil.hwetest[vil.hwetest$E.HOM1.HET.HOM2 !=   "nan/nan/nan",]

ggplot(vil.hwetest.cl, aes(x=P_HWE))+geom_density()

################################ ################################ ################################ ################################ 
################################ ################################ ################################ ################################ 
################################ ################################ ################################ ################################ 
################################ ################################ ################################ ################################ 
################################ ################################ ################################ ################################ 
#####################################
# Cophylo of Drakes
#note, poppr::popsub drops monomorphic alleles unless specified
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

hier.drake23<-popsub(hier.gi, sublist =c("Drake_2_2004", "Drake_3_2004"))

library(phangorn)
dr.aflp.tr<-(upgma(dist(hier.drake23)))









#setwd("~/Dropbox/finalizing_genets")
setwd("~/Dropbox/finalizing_genets")
maf006rds<-readRDS("Robjs.ma006.mx.5/maf006.mx.5.gi.rds")

myPops <- as.factor(c('old_EU',	'old_EU',	'old_EU',	'old_EU',	'old_EU',	'old_EU',	'CA_old',	'old_EU',	'D2.04',	'D2.04',	
                      'D2.04',	'D2.04',	'D2.04',	'D2.04',	'D2.04',	'D2.04',	'D2.04',	'D2.04',	'D2.04',	'D2.04',
                      'D2.04',	'D3.04',	'D3.04',	'D3.04',	'D3.04',	'D3.04',	'old_EU',	'D2.14',	'D2.14',	'D2.14',
                      'D2.14',	'D2.14',	'D2.14',	'D2.14',	'D2.14',	'D2.14',	'D2.14',	'D2.14',	'D2.14',	'D2.14',
                      'D2.14',	'D2.14',	'D2.14',	'D2.14',	'D2.14',	'D2.14',	'D2.14',	'D2.14',	'D2.14',	'D2.14',
                      'D2.14',	'D2.14',	'D3.14',	'D3.14',	'D3.14',	'D3.14',	'D3.14',	'D3.14',	'D3.14',	'D3.14',
                      'D3.14',	'Agraria',	'Agraria',	'Vilarinho',	'Vilarinho',	'Vilarinho',	'Vilarinho',	'Vilarinho',
                      'Mira',	'Mira',	'Mira',	'Mira',	'D2.15',	'D2.15',	'D2.15',	'D2.15',	'D2.15',	'D2.15',	'D2.15',
                      'D2.15',	'D2.15',	'D2.15',	'D2.15',	'D3.15',	'D3.15',	'D3.15'))


pop(maf006rds)<-myPops
maf006.dr<-popsub(maf006rds, sublist =c("D2.04", "D3.04"))


snp.dr.tr<-(upgma(dist(maf006.dr)))

library(phytools)
obj<-cophylo(dr.aflp.tr,snp.dr.tr,rotate=F)

plot(obj)
add.scale.bar()



