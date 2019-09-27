#############################################################################################################
#############################################################################################################
#############################################################################################################
################################# Genetic Distance Density Distr    #########################################
#############################################################################################################
#############################################################################################################
#############################################################################################################
library(dartR)
library(adegenet)

setwd("~/Dropbox/finalizing_genets")

# SNP data
############ 
############ 
maf006yall<-read.PLINK("HESS_fix/mx0.5/vcf.MMDP60.maf0.006.minQ30.recode.raw", 
                map.file = "HESS_fix/mx0.5/vcf.MMDP60.maf0.006.minQ30.recode.map", quiet = FALSE,
                parallel = require("parallel"), n.cores = NULL)



# maf006yall.dist<-dist(maf006yall, method="manhattan")
maf006yall.dist<-bitwise.dist(maf006yall, percent = F, missing_match = T, scale_missing = T)

#NOTE: decimals because dist takes missing into consideration by scaling by the number of rows used in the comparison

library(reshape2)
library(ggplot2)
maf006.dist<-as.matrix(maf006yall.dist)
maf006.dist[lower.tri(maf006.dist, diag=T)] <- "boob"

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

# info_table(maf006yall, type = c("missing"), percent = TRUE,
#            plot = FALSE, df = FALSE, returnplot = FALSE, low = "blue",
#            high = "red", plotlab = TRUE, scaled = TRUE)

labs<-cbind(  colnames(maf006.dist), as.data.frame(myPops))
labs$nl<-paste(labs$`colnames(maf006.dist)`, labs$myPops, sep="*")

colnames(maf006.dist)<-rownames(maf006.dist)<-labs$nl

maf006.dist.mlt<-melt(maf006.dist)
maf006.dist.mlt<-maf006.dist.mlt[maf006.dist.mlt$value != "boob",]
#write.csv(maf006.dist.mlt, "SNP-euclid.dist_melted.csv")

library(tidyr)
sample.ind1<-as.data.frame(separate(data = maf006.dist.mlt, col = Var1, into = c("sample1", "pop1"), sep = "\\*"))
sample.ind2<-as.data.frame(separate(data = maf006.dist.mlt, col = Var2, into = c("sample2", "pop2"), sep = "\\*"))

maf006.build<-cbind(maf006.dist.mlt, sample.ind1[, c("sample1", "pop1")]  , sample.ind2[, c("sample2", "pop2")])
maf006.build$pop.pair<-paste(maf006.build$pop1, maf006.build$pop2)


maf006.build$value <- as.numeric(as.character(maf006.build$value))


w.in.dr<-as.vector(c("D2.04 D2.04", "D2.14 D2.14", "D2.15 D2.15", "D3.04 D3.04" ,"D3.14 D3.14","D3.15 D3.15" ))
#betw.in.dr<-as.vector(c("D2.04 D3.04", "D2.04 D3.14","D2.04 D3.15","D2.14 D3.14","D2.14 D3.15", "D2.15 D3.15"))
#europe<-as.vector(c("Agraria Agraria","Agraria Mira",  "Agraria Vilarinho", "Vilarinho Mira","Vilarinho Vilarinho", "Mira Mira", "old_EU Vilarinho", "old_EU Agraria", "old_EU Mira", "old_EU old_EU"))
portugal<-as.vector(c("Agraria Agraria","Vilarinho Vilarinho", "Mira Mira"))

maf006.build.w.in.drakes<-maf006.build[maf006.build$pop.pair %in% w.in.dr,]
#maf006.build.betw.in.drakes<-maf006.build[maf006.build$pop.pair %in% betw.in.dr,]
#maf006.build.EU<-maf006.build[maf006.build$pop.pair %in% europe,]
maf006.build.port<-maf006.build[maf006.build$pop.pair %in% portugal,]


w.inpops.acrossyrs<-as.vector(c("D2.04 D2.04", "D2.14 D2.14", "D2.15 D2.15",
                                "D3.04 D3.04" ,"D3.14 D3.14", "D3.15 D3.15",
                                "D2.04 D2.14", "D2.04 D2.15", "D2.14 D2.15",
                                "D2.14 D2.04", "D2.15 D2.04", "D2.15 D2.14",
                                "D3.04 D3.14" ,"D3.04 D3.15", "D3.14 D3.15",
                                "D3.14 D3.04" ,"D3.15 D3.04", "D3.15 D3.14"))
maf006.Drakes<-maf006.build[maf006.build$pop.pair %in% w.inpops.acrossyrs,]


SNPtogo<-rbind(maf006.build.w.in.drakes, maf006.build.port)
SNPtogo$pop.pair <- factor(SNPtogo$pop.pair,  levels = c("D2.04 D2.04", "D2.14 D2.14", "D2.15 D2.15", "D3.04 D3.04", "D3.14 D3.14", "D3.15 D3.15", "Agraria Agraria", "Mira Mira", "Vilarinho Vilarinho"))

#plots

ggplot(SNPtogo,     aes(value, fill = pop.pair)) + 
  geom_histogram(center=-1, bins=90, position="stack")+theme_bw()+labs(y="Frequency", x="Genetic Distance", title="SNP Data", fill="")+
  scale_x_continuous(limits = c(0, 90000))+
  scale_fill_manual( values = c("#ff0000","#b30000", "#ff6666", "#ff9900", "#cc7a00", "#cc9900","#cc33ff","#9966ff","#6600ff"))

ggplot(SNPtogo,     aes(value, fill = pop.pair)) + 
  geom_histogram(center=-1, bins=90, position="stack")+theme_bw()+labs(y="Frequency", x="Genetic Distance", title="SNP Data", fill="")+
  scale_x_continuous(limits = c(0, 90000))+
  scale_fill_manual( values = c("#ff0000","#b30000", "#ff6666", "#ff9900", "#cc7a00", "#cc9900","#cc33ff","#9966ff","#6600ff"))+
  theme(legend.position = "none")

SNPmanh<-ggplot(SNPtogo,     aes(value, fill = pop.pair)) + 
  geom_histogram(center=-1, bins=90, position="stack")+theme_bw()+labs(y="Frequency", x="Genetic Distance", title="SNP Data", fill="")+
  scale_x_continuous(limits = c(0, 90000))+
  scale_fill_manual( values = c("#ff0000","#b30000", "#ff6666", "#ff9900", "#cc7a00", "#cc9900","#cc33ff","#9966ff","#6600ff"))

SNPmanh.nl<-ggplot(SNPtogo,     aes(value, fill = pop.pair)) + 
  geom_histogram(center=-1, bins=90, position="stack")+theme_bw()+labs(y="Frequency", x="Genetic Distance", title="SNP Data", fill="")+
  scale_x_continuous(limits = c(0, 90000))+
  scale_fill_manual( values = c("#ff0000","#b30000", "#ff6666", "#ff9900", "#cc7a00", "#cc9900","#cc33ff","#9966ff","#6600ff"))+
  theme(legend.position = "none")

##############################
############ 
############ 
############ AFLP
############ 
############ 
##############################
setwd("~/Dropbox/Genet Size/Jacob Analyses")
library(poppr)
hier <- read.genalex("AFLP_genalex_AMANITABASED.csv")

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

library(reshape2)
library(ggplot2)
hier.dist<-as.matrix(diss.dist(hier.gi, percent=F))

hier.dist[lower.tri(hier.dist, diag=T)] <- "boob"

labs<-cbind(  colnames(hier.dist), as.data.frame(pops.aflp))
labs$nl<-paste(labs$`colnames(hier.dist)`, labs$pops.aflp, sep="*")

colnames(hier.dist)<-rownames(hier.dist)<-labs$nl

hier.dist.mlt<-melt(hier.dist)
hier.dist.mlt<-hier.dist.mlt[hier.dist.mlt$value != "boob",]
#write.csv(maf006.dist.mlt, "SNP-euclid.dist_melted.csv")

library(tidyr)
sample.ind1<-as.data.frame(separate(data = hier.dist.mlt, col = Var1, into = c("sample1", "pop1"), sep = "\\*"))
sample.ind2<-as.data.frame(separate(data = hier.dist.mlt, col = Var2, into = c("sample2", "pop2"), sep = "\\*"))

hier.dist.build<-cbind(hier.dist.mlt, sample.ind1[, c("sample1", "pop1")]  , sample.ind2[, c("sample2", "pop2")])
hier.dist.build$pop.pair<-paste(hier.dist.build$pop1, hier.dist.build$pop2)


hier.dist.build$value <- as.numeric(as.character(hier.dist.build$value))

roch<-as.vector(c("Rochester_NY_2007_1 Rochester_NY_2007_1","Rochester_NY_2007_2 Rochester_NY_2007_2", "Rochester_NY_2007_3 Rochester_NY_2007_3"))
hier.dist.build.roch<-hier.dist.build[hier.dist.build$pop.pair %in% roch,]

drakes.win<-as.vector(c("Drake_1_2004 Drake_1_2004",  "Drake_2_2004 Drake_2_2004" , "Drake_3_2004 Drake_3_2004", "Drake_4_2004 Drake_4_2004" ))
hier.dist.build.drakes.win<-hier.dist.build[hier.dist.build$pop.pair %in% drakes.win,]

#drakes.wout<-as.vector(c("Drake_1_2004 Drake_2_2004","Drake_1_2004 Drake_3_2004","Drake_1_2004 Drake_4_2004","Drake_2_2004 Drake_3_2004","Drake_2_2004 Drake_4_2004","Drake_3_2004 Drake_4_2004"))
#hier.dist.build.build.drakes.wout<-hier.dist.build[hier.dist.build$pop.pair %in% drakes.wout,]

nj.winwout<-as.vector(c("Jakes_Landing_2006 Jakes_Landing_2006","New_Jersey_str_2006 New_Jersey_str_2006"))
hier.dist.build.build.nj.winwout<-hier.dist.build[hier.dist.build$pop.pair %in% nj.winwout,]

hd.winwout<-as.vector(c("Hearts_Desire_1_2004 Hearts_Desire_1_2004", 
                        "Hearts_Desire_2_2004 Hearts_Desire_2_2004","Hearts_Desire_3_2004 Hearts_Desire_3_2004"))
hier.dist.build.build.hd.winwout<-hier.dist.build[hier.dist.build$pop.pair %in% hd.winwout,]

eu.winwout<-as.vector(c("CESAC_2002 CESAC_2002","CESAC_2006 CESAC_2006","Serbia_2007 Serbia_2007"))
hier.dist.build.build.eu.winwout<-hier.dist.build[hier.dist.build$pop.pair %in% eu.winwout,]




AFLPtogo<-rbind(hier.dist.build.roch, hier.dist.build.drakes.win, hier.dist.build.build.nj.winwout, hier.dist.build.build.hd.winwout, hier.dist.build.build.eu.winwout)
SNPtogo$pop.pair <- factor(SNPtogo$pop.pair,  levels = c(
  "Rochester_NY_2007_1 Rochester_NY_2007_1","Rochester_NY_2007_2 Rochester_NY_2007_2","Rochester_NY_2007_3 Rochester_NY_2007_3",
  "Drake_1_2004 Drake_1_2004","Drake_2_2004 Drake_2_2004","Drake_3_2004 Drake_3_2004","Drake_4_2004 Drake_4_2004",
  "Jakes_Landing_2006 Jakes_Landing_2006","New_Jersey_str_2006 New_Jersey_str_2006",
  "Hearts_Desire_1_2004 Hearts_Desire_1_2004","Hearts_Desire_2_2004 Hearts_Desire_2_2004","Hearts_Desire_3_2004 Hearts_Desire_3_2004",
  "CESAC_2002 CESAC_2002","CESAC_2006 CESAC_2006","Serbia_2007 Serbia_2007"
))


# Plots
ggplot(AFLPtogo,     aes(value, fill = pop.pair)) + 
  geom_histogram(bins=33, position="stack")+theme_bw()+labs(y="Frequency", x="Genetic Distance", title="AFLP Data", fill="")+
  #scale_x_continuous(limits = c(-.5, 35))+
  scale_fill_manual( values = c(
    "#9900cc","#993366","#ff6699","#ff0000","#ff9900","#cc5200","#ffff00","#999900","#e6b800","#0000ff","#8080ff","#0099ff","#006bb3","#b3e0ff","#df9fbf"
  ))


ggplot(AFLPtogo,     aes(value, fill = pop.pair)) + 
  geom_histogram(bins=33, position="stack")+theme_bw()+labs(y="Frequency", x="Genetic Distance", title="AFLP Data", fill="")+
  #scale_x_continuous(limits = c(-.5, 40))+
  scale_fill_manual( values = c(
    "#9900cc","#993366","#ff6699","#ff0000","#ff9900","#cc5200","#ffff00","#999900","#e6b800","#0000ff","#8080ff","#0099ff","#006bb3","#b3e0ff","#df9fbf"
  ))+theme(legend.position = "none")

AFLPmanh<-ggplot(AFLPtogo,     aes(value, fill = pop.pair)) + 
  geom_histogram(bins=33, position="stack")+theme_bw()+labs(y="Frequency", x="Genetic Distance", title="AFLP Data", fill="")+
  #scale_x_continuous(limits = c(-.5, 35))+
  scale_fill_manual( values = c(
    "#9900cc","#993366","#ff6699","#ff0000","#ff9900","#cc5200","#ffff00","#999900","#e6b800","#0000ff","#8080ff","#0099ff","#006bb3","#b3e0ff","#df9fbf"
  ))


AFLPmanh.nl<-ggplot(AFLPtogo,     aes(value, fill = pop.pair)) + 
  geom_histogram(bins=33, position="stack")+theme_bw()+labs(y="Frequency", x="Genetic Distance", title="AFLP Data", fill="")+
  #scale_x_continuous(limits = c(-.5, 40))+
  scale_fill_manual( values = c(
    "#9900cc","#993366","#ff6699","#ff0000","#ff9900","#cc5200","#ffff00","#999900","#e6b800","#0000ff","#8080ff","#0099ff","#006bb3","#b3e0ff","#df9fbf"
  ))+theme(legend.position = "none")


# Create Panel

library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
r <- rectGrob(gp=gpar(fill="white"))

grid.arrange(AFLPmanh.nl, r, SNPmanh.nl,r, nrow=2, ncol=2)
grid.arrange(AFLPmanh, r, SNPmanh,r, nrow=2, ncol=2)
# grid.arrange(a.l,b.l,d.l,e.l,f.l, ncol=1, nrow=5)



###################################################################
# averages

#simulate
test<-glSim(n.ind=2, n.snp.nonstruc = 10)
View(as.matrix(test))

bitwise.dist(test, percent = F)

e1071::hamming.distance(rbind( c(0,1,1,1,0,1,0,0,0,0), 
                               c(1,1,0,1,0,2,0,0,1,1) ))

x <- new("genlight", list(indiv1=c(1,1,0,1,1,0), 
                          indiv2=c(1,1,1,2,0,0)), ploidy=2)

bitwise.dist(x, percent = F, missing_match=T, scale_missing = T)

y <- new("genlight", list(indiv1=c(1,1, 0), 
                          indiv2=c(2, 2, 2)), ploidy=2)
bitwise.dist(y, percent = F, missing_match=T, scale_missing = T)

#MIRA
MIRA<-popsub(maf006yall, sublist=c("Mira"))
MIRAcl<-MIRA[2:3,]

bitwise.dist(MIRAcl, percent = F, missing_match=T, scale_missing = T) #24250.29

#same because I am ignoring missing data, decimal because scaling for missing data

head(as.matrix(MIRAnm[,1:10]))
unique(as.vector(as.matrix(MIRAcl[1,])))
unique(as.vector(as.matrix(MIRAcl[2,])))




# average dist in pops and clsuters
setwd("~/Dropbox/Genet Size/Jacob Analyses")
coords.relab<-read.csv("subpops.csv", head=T, stringsAsFactors = F)

maf006yall.dr2s<-popsub(maf006yall, sublist =c("D2.04", "D2.14", "D2.15"))

#d2.14 cluster 1
maf006yall.dr2s.cl1 <- maf006yall.dr2s[(indNames(maf006yall.dr2s) %in% coords.relab[coords.relab$subpop %in% c("D2.14"), ]$Sample_Name)]
maf006yall.dr2s.cl1.mat<-as.matrix(bitwise.dist(maf006yall.dr2s.cl1, percent = F, missing_match=T, scale_missing = T))
maf006yall.dr2s.cl1.mat[lower.tri(maf006yall.dr2s.cl1.mat, diag=T)] <- NA
maf006yall.dr2s.cl1.mat.melt<-melt(maf006yall.dr2s.cl1.mat)
maf006yall.dr2s.cl1.mat.melt.cl<-na.omit(maf006yall.dr2s.cl1.mat.melt)
mean(maf006yall.dr2s.cl1.mat.melt.cl$value)
hist(maf006yall.dr2s.cl1.mat.melt.cl$value)
plotrix::std.error((maf006yall.dr2s.cl1.mat.melt.cl$value))
sd((maf006yall.dr2s.cl1.mat.melt.cl$value))

#d2.14 cluster 2
maf006yall.dr2s.cl2 <- maf006yall.dr2s[(indNames(maf006yall.dr2s) %in% coords.relab[coords.relab$subpop %in% c("D2s.14"), ]$Sample_Name)]
maf006yall.dr2s.cl2.mat<-as.matrix(bitwise.dist(maf006yall.dr2s.cl2, percent = F, missing_match=T, scale_missing = T))
maf006yall.dr2s.cl2.mat[lower.tri(maf006yall.dr2s.cl2.mat, diag=T)] <- NA
maf006yall.dr2s.cl2.mat.melt<-melt(maf006yall.dr2s.cl2.mat)
maf006yall.dr2s.cl2.mat.melt.cl<-na.omit(maf006yall.dr2s.cl2.mat.melt)
mean(maf006yall.dr2s.cl2.mat.melt.cl$value)
hist(maf006yall.dr2s.cl2.mat.melt.cl$value)
plotrix::std.error((maf006yall.dr2s.cl2.mat.melt.cl$value))
sd((maf006yall.dr2s.cl2.mat.melt.cl$value))

#d2.15 cluster 1
maf006yall.dr2s15.cl1 <- maf006yall.dr2s[(indNames(maf006yall.dr2s) %in% coords.relab[coords.relab$subpop %in% c("D2.15"), ]$Sample_Name)]
maf006yall.dr2s15.cl1.mat<-as.matrix(bitwise.dist(maf006yall.dr2s15.cl1, percent = F, missing_match=T, scale_missing = T))
maf006yall.dr2s15.cl1.mat[lower.tri(maf006yall.dr2s15.cl1.mat, diag=T)] <- NA
maf006yall.dr2s15.cl1.mat.melt<-melt(maf006yall.dr2s15.cl1.mat)
maf006yall.dr2s15.cl1.mat.melt.cl<-na.omit(maf006yall.dr2s15.cl1.mat.melt)
mean(maf006yall.dr2s15.cl1.mat.melt.cl$value)
hist(maf006yall.dr2s15.cl1.mat.melt.cl$value)
plotrix::std.error((maf006yall.dr2s15.cl1.mat.melt.cl$value))
sd((maf006yall.dr2s15.cl1.mat.melt.cl$value))

#d2.15 cluster 2
maf006yall.dr2s15.cl2 <- maf006yall.dr2s[(indNames(maf006yall.dr2s) %in% coords.relab[coords.relab$subpop %in% c("D2s.15"), ]$Sample_Name)]
maf006yall.dr2s15.cl2.mat<-as.matrix(bitwise.dist(maf006yall.dr2s15.cl2, percent = F, missing_match=T, scale_missing = T))
maf006yall.dr2s15.cl2.mat[lower.tri(maf006yall.dr2s15.cl2.mat, diag=T)] <- NA
maf006yall.dr2s15.cl2.mat.melt<-melt(maf006yall.dr2s15.cl2.mat)
maf006yall.dr2s15.cl2.mat.melt.cl<-na.omit(maf006yall.dr2s15.cl2.mat.melt)
mean(maf006yall.dr2s15.cl2.mat.melt.cl$value)
hist(maf006yall.dr2s15.cl2.mat.melt.cl$value)
plotrix::std.error((maf006yall.dr2s15.cl2.mat.melt.cl$value))
sd((maf006yall.dr2s15.cl2.mat.melt.cl$value))






















