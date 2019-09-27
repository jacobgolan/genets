##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
########                                          MANTEL TESTS
######## Realized I can only do mantel tests within populations because we only map within pops, not between


setwd("~/Dropbox/finalizing_genets")

MAF006.mira.gi.rds<-read.PLINK("HESS_fix/pops/vcf.MMDP60.maf0.006.minQ30.recode.mira.HWE.raw", 
                                                    map.file = "HESS_fix/pops/vcf.MMDP60.maf0.006.minQ30.recode.mira.HWE.map", quiet = FALSE,
                                                    parallel = require("parallel"), n.cores = NULL)
MAF006.agraria.gi.rds<-read.PLINK("HESS_fix/pops/vcf.MMDP60.maf0.006.minQ30.recode.agraria.HWE.raw", 
                                                    map.file = "HESS_fix/pops/vcf.MMDP60.maf0.006.minQ30.recode.agraria.HWE.map", quiet = FALSE,
                                                    parallel = require("parallel"), n.cores = NULL)
MAF006.vilarinho.gi.rds<-read.PLINK("HESS_fix/pops/vcf.MMDP60.maf0.006.minQ30.recode.vilarinho.HWE.raw", 
                                                    map.file = "HESS_fix/pops/vcf.MMDP60.maf0.006.minQ30.recode.vilarinho.HWE.map", quiet = FALSE,
                                                    parallel = require("parallel"), n.cores = NULL)

MAF006.D2.04.gi.rds<-read.PLINK("HESS_fix/pops/vcf.MMDP60.maf0.006.minQ30.recode.D2.04.HWE.raw", 
                                                    map.file = "HESS_fix/pops/vcf.MMDP60.maf0.006.minQ30.recode.D2.04.HWE.map", quiet = FALSE,
                                                    parallel = require("parallel"), n.cores = NULL)
MAF006.D2.14.gi.rds<-read.PLINK("HESS_fix/pops/vcf.MMDP60.maf0.006.minQ30.recode.D2.14.HWE.raw", 
                                                    map.file = "HESS_fix/pops/vcf.MMDP60.maf0.006.minQ30.recode.D2.14.HWE.map", quiet = FALSE,
                                                    parallel = require("parallel"), n.cores = NULL)
MAF006.D2.15.gi.rds<-read.PLINK("HESS_fix/pops/vcf.MMDP60.maf0.006.minQ30.recode.D2.15.HWE.raw", 
                                                    map.file = "HESS_fix/pops/vcf.MMDP60.maf0.006.minQ30.recode.D2.15.HWE.map", quiet = FALSE,
                                                    parallel = require("parallel"), n.cores = NULL)

MAF006.D3.04.gi.rds<-read.PLINK("HESS_fix/pops/vcf.MMDP60.maf0.006.minQ30.recode.D3.04.HWE.raw", 
                                                    map.file = "HESS_fix/pops/vcf.MMDP60.maf0.006.minQ30.recode.D3.04.HWE.map", quiet = FALSE,
                                                    parallel = require("parallel"), n.cores = NULL)
MAF006.D3.14.gi.rds<-read.PLINK("HESS_fix/pops/vcf.MMDP60.maf0.006.minQ30.recode.D3.14.HWE.raw", 
                                                    map.file = "HESS_fix/pops/vcf.MMDP60.maf0.006.minQ30.recode.D3.14.HWE.map", quiet = FALSE,
                                                    parallel = require("parallel"), n.cores = NULL)
MAF006.D3.15.gi.rds<-read.PLINK("HESS_fix/pops/vcf.MMDP60.maf0.006.minQ30.recode.D3.15.HWE.raw", 
                                                    map.file = "HESS_fix/pops/vcf.MMDP60.maf0.006.minQ30.recode.D3.15.HWE.map", quiet = FALSE,
                                                    parallel = require("parallel"), n.cores = NULL)



#########################################
# Genetic distance

#SNP
MAF006.mira.dist<-bitwise.dist(MAF006.mira.gi.rds, percent = F, missing_match = T, scale_missing = T)
MAF006.agraria.dist<-bitwise.dist(MAF006.agraria.gi.rds, percent = F, missing_match = T, scale_missing = T)

MAF006.D2.04.dist<-bitwise.dist(MAF006.D2.04.gi.rds, percent = F, missing_match = T, scale_missing = T)
MAF006.D2.14.dist<-bitwise.dist(MAF006.D2.14.gi.rds, percent = F, missing_match = T, scale_missing = T)
MAF006.D2.15.dist<-bitwise.dist(MAF006.D2.15.gi.rds, percent = F, missing_match = T, scale_missing = T)
MAF006.D3.04.dist<-bitwise.dist(MAF006.D3.04.gi.rds, percent = F, missing_match = T, scale_missing = T)
MAF006.D3.14.dist<-bitwise.dist(MAF006.D3.14.gi.rds, percent = F, missing_match = T, scale_missing = T)
MAF006.D3.15.dist<-bitwise.dist(MAF006.D3.15.gi.rds, percent = F, missing_match = T, scale_missing = T)

#####################
#AFLP
setwd("~/Dropbox/Genet Size/Jacob Analyses")
library(poppr)
hier<-read.genalex("AFLP_genalex_AMANITABASED.csv")
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

hier.gi<-as.genind(as.matrix(hier))

pop(hier.gi)<-pops.aflp
hier.gc<-as.genclone(hier.gi)

hier.CESAC02<-popsub(hier.gi, sublist = "CESAC_2002")
hier.CESAC06<-popsub(hier.gi, sublist = "CESAC_2006")

hier.Roch1<-popsub(hier.gi, sublist = "Rochester_NY_2007_1")
hier.Roch2<-popsub(hier.gi, sublist = "Rochester_NY_2007_2")
hier.Roch2.sub<-hier.Roch2[c("10124" , "10103"  ,"10107",  "10105" , "10122", "10123"),] 
hier.Roch3<-popsub(hier.gi, sublist = "Rochester_NY_2007_3")

hier.d1<-popsub(hier.gi, sublist = "Drake_1_2004")
hier.d2<-popsub(hier.gi, sublist = "Drake_2_2004")
hier.d3<-popsub(hier.gi, sublist = "Drake_3_2004")
hier.d4<-popsub(hier.gi, sublist = "Drake_4_2004")

hier.hd1<-popsub(hier.gi, sublist = "Hearts_Desire_1_2004")
hier.hd2<-popsub(hier.gi, sublist = "Hearts_Desire_2_2004")
hier.hd3<-popsub(hier.gi, sublist = "Hearts_Desire_3_2004")

hier.JL<-popsub(hier.gi, sublist = "Jakes_Landing_2006")
hier.str<-popsub(hier.gi, sublist = "New_Jersey_str_2006")

hier.S<-popsub(hier.gi, sublist = "Serbia_2007")
hier.S.sub <- hier.S[c("Aph265", "Aph274", "Aph266", "Aph267", "Aph268", "Aph269", "Aph270", "Aph271", "Aph272", "Aph273"),]

#________________________

hier.CESAC02.dist<-diss.dist(hier.CESAC02, percent = F)
hier.CESAC06.dist<-diss.dist(hier.CESAC06, percent = F)

hier.Roch1.dist<-diss.dist(hier.Roch1, percent = F)
hier.Roch2.dist<-diss.dist(hier.Roch2, percent = F)
hier.Roch2.sub.dist<-diss.dist(hier.Roch2.sub, percent = F)
hier.Roch3.dist<-diss.dist(hier.Roch3, percent = F)

hier.d1.dist<-diss.dist(hier.d1, percent = F)
hier.d2.dist<-diss.dist(hier.d2, percent = F)
hier.d3.dist<-diss.dist(hier.d3, percent = F,)
hier.d4.dist<-diss.dist(hier.d4, percent = F)

hier.hd1.dist<-diss.dist(hier.hd1, percent = F)
hier.hd2.dist<-diss.dist(hier.hd2, percent = F)
hier.hd3.dist<-diss.dist(hier.hd3, percent = F)

hier.JL.dist<-diss.dist(hier.JL, percent = F)
hier.str.dist<-diss.dist(hier.str, percent = F)

hier.S.dist<-diss.dist(hier.S, percent = F)
hier.S.sub.dist<-diss.dist(hier.S.sub, percent = F)

# SNP
#########################################
#geo distance
setwd("~/Dropbox/GENETS_finishing")
coords<-read.csv("ALL_COORDS_geo.csv", head=T, stringsAsFactors = F)

coords.d1<-coords[coords$pop =="Drake1",]
coords.d1.aflp<-coords.d1

coords.d2<-coords[coords$pop =="Drake2",]
coords.d2.04<-coords.d2[coords.d2$year =="Sample_2004",]
coords.d2.04.aflp<-coords.d2.04[-6,]
coords.d2.04.snp<-coords.d2.04[-1,]


coords.d2.14<-coords.d2[coords.d2$year =="Sample_2014",]
coords.d2.14.snp<-coords.d2.14


coords.d2.15<-coords.d2[coords.d2$year =="Sample_2015",]
coords.d2.15.snp<-coords.d2.15[-9,]

coords.d3<-coords[coords$pop =="Drake3",]
coords.d3.04<-coords.d3[coords.d3$year =="Sample_2004",]
coords.d3.04.aflp<-coords.d3.04
coords.d3.04.snp<-coords.d3.04[3:7,]

coords.d3.14<-coords.d3[coords.d3$year =="Sample_2014",]
coords.d3.14.snp<-coords.d3.14

coords.d3.15<-coords.d3[coords.d3$year =="Sample_2015",]
coords.d3.15.snp<-coords.d3.15


coords.d4<-coords[coords$pop =="Drake4",]
coords.d4.aflp<-coords.d4[-2,]

coords.hd1<-coords[coords$pop =="HD1",]
coords.hd1.aflp<-coords.hd1
coords.hd2<-coords[coords$pop =="HD2",]
coords.hd2.aflp<-coords.hd2[-4,]
coords.hd3<-coords[coords$pop =="HD3",]
coords.hd3.aflp<-coords.hd3

coords.mira<-coords[coords$pop =="mira",]
coords.mira.snp<-coords.mira
coords.agraria<-coords[coords$pop =="agraria",]
coords.agraria.snp<-coords.agraria

coords.JL<-coords[coords$pop =="JL",]
coords.JL.aflp<-coords.JL

coords.CESAC02<-coords[coords$pop =="CESAC02",]
coords.CESAC02.aflp<-coords.CESAC02
coords.CESAC06<-coords[coords$pop =="CESAC06",]
coords.CESAC06.aflp<-coords.CESAC06

coords.S<-coords[coords$pop =="S",]
coords.S.aflp<-coords.S # note S 11-14 do not have coords, need to remove from genind

coords.str<-coords[coords$pop =="str",]
coords.str.aflp<-coords.str[-c(2,6,8,12,15,16),]

coords.roch<-coords[coords$pop =="roch",]
coords.roch1<-coords.roch[1:13,]
coords.roch1.aflp<-coords.roch1
coords.roch2<-coords.roch[14:20,]
coords.roch2.aflp<-coords.roch2[-4,] #need to remove 2_6b from genind file
coords.roch3<-coords.roch[21:30,]
coords.roch3.aflp<-coords.roch3[3:4,]

####
#AFLP
coords.d1.aflp
coords.d2.04.aflp
coords.d3.04.aflp
coords.d4.aflp

coords.hd1.aflp
coords.hd2.aflp
coords.hd3.aflp

coords.JL.aflp

coords.CESAC02.aflp
coords.CESAC06.aflp

coords.S.aflp

coords.str.aflp

coords.roch1.aflp

coords.roch2.aflp

coords.roch3.aflp


#SNP

coords.d2.04.snp
coords.d2.14.snp
coords.d2.15.snp

coords.d3.04.snp
coords.d3.14.snp
coords.d3.15.snp

coords.mira.snp
coords.agraria.snp


#rownames

rownames(coords.d1.aflp)<-coords.d1.aflp$FINALlab
rownames(coords.d2.04.aflp)<-coords.d2.04.aflp$FINALlab
rownames(coords.d3.04.aflp)<-coords.d3.04.aflp$FINALlab
rownames(coords.d4.aflp)<-coords.d4.aflp$FINALlab

rownames(coords.hd1.aflp)<-coords.hd1.aflp$FINALlab
rownames(coords.hd2.aflp)<-coords.hd2.aflp$FINALlab
rownames(coords.hd3.aflp)<-coords.hd3.aflp$FINALlab

rownames(coords.JL.aflp)<-coords.JL.aflp$FINALlab

rownames(coords.CESAC02.aflp)<-coords.CESAC02.aflp$FINALlab
rownames(coords.CESAC06.aflp)<-coords.CESAC06.aflp$FINALlab

rownames(coords.S.aflp)<-coords.S.aflp$FINALlab

rownames(coords.str.aflp)<-coords.str.aflp$FINALlab

rownames(coords.roch1.aflp)<-coords.roch1.aflp$FINALlab

rownames(coords.roch2.aflp)<-coords.roch2.aflp$FINALlab

rownames(coords.roch3.aflp)<-coords.roch3.aflp$FINALlab



rownames(coords.d2.04.snp)<-coords.d2.04.snp$FINALlab
rownames(coords.d2.14.snp)<-coords.d2.14.snp$FINALlab
rownames(coords.d2.15.snp)<-coords.d2.15.snp$FINALlab

rownames(coords.d3.04.snp)<-coords.d3.04.snp$FINALlab
rownames(coords.d3.14.snp)<-coords.d3.14.snp$FINALlab
rownames(coords.d3.15.snp)<-coords.d3.15.snp$FINALlab

rownames(coords.mira.snp)<-coords.mira.snp$FINALlab
rownames(coords.agraria.snp)<-coords.agraria.snp$FINALlab

# Euclidean distance of x and y

coords.d1.aflp.dist<-as.matrix(dist(coords.d1.aflp[,5:6]))
colnames(coords.d1.aflp.dist)<-rownames(coords.d1.aflp.dist)<-coords.d1.aflp$FINAlab

coords.d2.04.aflp.dist<-as.matrix(dist(coords.d2.04.aflp[,5:6]))
colnames(coords.d2.04.aflp.dist)<-rownames(coords.d2.04.aflp.dist)<-coords.d2.04.aflp$FINAlab

coords.d3.04.aflp.dist<-as.matrix(dist(coords.d3.04.aflp[,5:6]))
colnames(coords.d3.04.aflp.dist)<-rownames(coords.d3.04.aflp.dist)<-coords.d3.04.aflp$FINAlab


coords.d4.aflp.dist<-as.matrix(dist(coords.d4.aflp[,5:6]))
colnames(coords.d4.aflp.dist)<-rownames(coords.d4.aflp.dist)<-coords.d4.aflp$FINAlab


coords.hd1.aflp.dist<-as.matrix(dist(coords.hd1.aflp[,5:6]))
colnames(coords.hd1.aflp.dist)<-rownames(coords.hd1.aflp.dist)<-coords.hd1.aflp$FINAlab

coords.hd2.aflp.dist<-as.matrix(dist(coords.hd2.aflp[,5:6]))
colnames(coords.hd2.aflp.dist)<-rownames(coords.hd2.aflp.dist)<-coords.hd2.aflp$FINAlab

coords.hd3.aflp.dist<-as.matrix(dist(coords.hd3.aflp[,5:6]))
colnames(coords.hd3.aflp.dist)<-rownames(coords.hd3.aflp.dist)<-coords.hd3.aflp$FINAlab


coords.JL.aflp.dist<-as.matrix(dist(coords.JL.aflp[,5:6]))
colnames(coords.JL.aflp.dist)<-rownames(coords.JL.aflp.dist)<-coords.JL.aflp$FINAlab


coords.CESAC02.aflp.dist<-as.matrix(dist(coords.CESAC02.aflp[,5:6]))
colnames(coords.CESAC02.aflp.dist)<-rownames(coords.CESAC02.aflp.dist)<-coords.CESAC02.aflp$FINAlab

coords.CESAC06.aflp.dist<-as.matrix(dist(coords.CESAC06.aflp[,5:6]))
colnames(coords.CESAC06.aflp.dist)<-rownames(coords.CESAC06.aflp.dist)<-coords.CESAC06.aflp$FINAlab


coords.S.aflp.dist<-as.matrix(dist(coords.S.aflp[,5:6]))
colnames(coords.S.aflp.dist)<-rownames(coords.S.aflp.dist)<-coords.S.aflp$FINAlab


coords.str.aflp.dist<-as.matrix(dist(coords.str.aflp[,5:6]))
colnames(coords.str.aflp.dist)<-rownames(coords.str.aflp.dist)<-coords.str.aflp$FINAlab


coords.roch1.aflp.dist<-as.matrix(dist(coords.roch1.aflp[,5:6]))
colnames(coords.roch1.aflp.dist)<-rownames(coords.roch1.aflp.dist)<-coords.roch1.aflp$FINAlab


coords.roch2.aflp.dist<-as.matrix(dist(coords.roch2.aflp[,5:6]))
colnames(coords.roch2.aflp.dist)<-rownames(coords.roch2.aflp.dist)<-coords.roch2.aflp$FINAlab


coords.roch3.aflp.dist<-as.matrix(dist(coords.roch3.aflp[,5:6]))
colnames(coords.roch3.aflp.dist)<-rownames(coords.roch3.aflp.dist)<-coords.roch3.aflp$FINAlab




coords.d2.04.snp.dist<-as.matrix(dist(coords.d2.04.snp[,5:6]))
colnames(coords.d2.04.snp.dist)<-rownames(coords.d2.04.snp.dist)<-coords.d2.04.snp$FINAlab

coords.d2.14.snp.dist<-as.matrix(dist(coords.d2.14.snp[,5:6]))
colnames(coords.d2.14.snp.dist)<-rownames(coords.d2.14.snp.dist)<-coords.d2.14.snp$FINAlab

coords.d2.15.snp.dist<-as.matrix(dist(coords.d2.15.snp[,5:6]))
colnames(coords.d2.15.snp.dist)<-rownames(coords.d2.15.snp.dist)<-coords.d2.15.snp$FINAlab


coords.d3.04.snp.dist<-as.matrix(dist(coords.d3.04.snp[,5:6]))
colnames(coords.d3.04.snp.dist)<-rownames(coords.d3.04.snp.dist)<-coords.d3.04.snp$FINAlab

coords.d3.14.snp.dist<-as.matrix(dist(coords.d3.14.snp[,5:6]))
colnames(coords.d3.14.snp.dist)<-rownames(coords.d3.14.snp.dist)<-coords.d3.14.snp$FINAlab


coords.d3.15.snp.dist<-as.matrix(dist(coords.d3.15.snp[,5:6]))
colnames(coords.d3.15.snp.dist)<-rownames(coords.d3.15.snp.dist)<-coords.d3.15.snp$FINAlab


coords.mira.snp.dist<-as.matrix(dist(coords.mira.snp[,5:6]))
colnames(coords.mira.snp.dist)<-rownames(coords.mira.snp.dist)<-coords.mira.snp$FINAlab

coords.agraria.snp.dist<-as.matrix(dist(coords.agraria.snp[,5:6]))
colnames(coords.agraria.snp.dist)<-rownames(coords.agraria.snp.dist)<-coords.agraria.snp$FINAlab



##############################################################################################
library(ape)
library(dartR)
library(vegan)

mantel.test(as.matrix(coords.d1.aflp.dist),as.matrix(hier.d1.dist), alternative = "greater")
gl.ibd(gl=NULL, Dgen=(hier.d1.dist), Dgeo=as.dist(coords.d1.aflp.dist))
mantel(coords.d1.aflp.dist, hier.d1.dist, method="p")

mantel.test(as.matrix(coords.d2.04.aflp.dist),as.matrix(hier.d2.dist), alternative = "greater")
gl.ibd(gl=NULL, Dgen=(hier.d2.dist), Dgeo=as.dist(coords.d2.04.aflp.dist))
mantel(coords.d2.04.aflp.dist, hier.d2.dist, method="p")

mantel.test(as.matrix(coords.d3.04.aflp.dist),as.matrix(hier.d3.dist), alternative = "greater")
gl.ibd(gl=NULL, Dgen=(hier.d3.dist), Dgeo=as.dist(coords.d3.04.aflp.dist))
mantel(coords.d3.04.aflp.dist, hier.d3.dist, method="p")

mantel.test(as.matrix(coords.d4.aflp.dist),as.matrix(hier.d4.dist), alternative = "greater")
gl.ibd(gl=NULL, Dgen=(hier.d4.dist), Dgeo=as.dist(coords.d4.aflp.dist)) 
mantel(coords.d4.aflp.dist, hier.d4.dist, method="p")

mantel.test(as.matrix(coords.hd1.aflp.dist),as.matrix(hier.hd1.dist), alternative = "greater")
gl.ibd(gl=NULL, Dgen=(hier.hd1.dist), Dgeo=as.dist(coords.hd1.aflp.dist))
mantel(coords.hd1.aflp.dist, hier.hd1.dist, method="p")



#need to reorder matrices so that they conform
hd2.geo<-as.matrix(coords.hd2.aflp.dist)
hd2.gen<-as.matrix(hier.hd2.dist)

hd2.gen <- hd2.gen[, order(as.integer(colnames(hd2.gen)))]
hd2.gen <- hd2.gen[order(as.integer(colnames(hd2.gen))),]

hd2.geo <- hd2.geo[, order(as.integer(colnames(hd2.geo)))]
hd2.geo <- hd2.geo[order(as.integer(colnames(hd2.geo))),]


mantel.test(hd2.gen, hd2.geo, alternative = "greater", nperm=10000)
gl.ibd(gl=NULL, Dgen=as.dist(hd2.gen), Dgeo=as.dist(hd2.geo), permutations = 10000)
mantel(hd2.geo, hd2.gen, method="p", permutations = 10000)


mantel.test(as.matrix(coords.hd3.aflp.dist),as.matrix(hier.hd3.dist), alternative = "greater")
gl.ibd(gl=NULL, Dgen=(hier.hd3.dist), Dgeo=as.dist(coords.hd3.aflp.dist))
mantel(coords.hd3.aflp.dist, hier.hd3.dist, method="p")


mantel.test(as.matrix(coords.JL.aflp.dist),as.matrix(hier.JL.dist), alternative = "greater")
gl.ibd(gl=NULL, Dgen=(hier.JL.dist), Dgeo=as.dist(coords.JL.aflp.dist))
mantel(coords.JL.aflp.dist, hier.JL.dist, method="p")

mantel.test(as.matrix(coords.CESAC02.aflp.dist),as.matrix(hier.CESAC02.dist), alternative = "greater", nperm=10000)
gl.ibd(gl=NULL, Dgen=(hier.CESAC02.dist), Dgeo=as.dist(coords.CESAC02.aflp.dist), permutations = 10000)
mantel(coords.CESAC02.aflp.dist, hier.CESAC02.dist, method="p", permutations = 10000)

mantel.test(as.matrix(coords.CESAC06.aflp.dist),as.matrix(hier.CESAC06.dist), alternative = "greater")
gl.ibd(gl=NULL, Dgen=(hier.CESAC06.dist), Dgeo=as.dist(coords.CESAC06.aflp.dist))
mantel(coords.CESAC06.aflp.dist, hier.CESAC06.dist, method="p")

mantel.test(as.matrix(coords.S.aflp.dist),as.matrix(hier.S.sub.dist), alternative = "greater")
gl.ibd(gl=NULL, Dgen=(hier.S.sub.dist), Dgeo=as.dist(coords.S.aflp.dist))
mantel(coords.S.aflp.dist, hier.S.sub.dist, method="p")

mantel.test(as.matrix(coords.str.aflp.dist),as.matrix(hier.str.dist), alternative = "greater")
gl.ibd(gl=NULL, Dgen=(hier.str.dist), Dgeo=as.dist(coords.str.aflp.dist))
mantel(coords.str.aflp.dist, hier.str.dist, method="p")

mantel.test(as.matrix(coords.roch1.aflp.dist),as.matrix(hier.Roch1.dist), alternative = "greater")
gl.ibd(gl=NULL, Dgen=(hier.Roch1.dist), Dgeo=as.dist(coords.roch1.aflp.dist))
mantel(coords.roch1.aflp.dist, hier.Roch1.dist, method="p")

mantel.test(as.matrix(coords.roch2.aflp.dist),as.matrix(hier.Roch2.sub.dist), alternative = "greater", nperm = 10000)
gl.ibd(gl=NULL, Dgen=(hier.Roch2.sub.dist), Dgeo=as.dist(coords.roch2.aflp.dist), permutations=10000)
mantel(coords.roch2.aflp.dist, hier.Roch2.sub.dist, method="p", permutations=10000)

mantel.test(as.matrix(coords.roch3.aflp.dist),as.matrix(hier.Roch3.dist), alternative = "greater")
gl.ibd(gl=NULL, Dgen=(hier.Roch3.dist), Dgeo=as.dist(coords.roch3.aflp.dist))
mantel(coords.roch3.aflp.dist, hier.Roch3.dist, method="p")

# SNP
mantel.test(as.matrix(coords.d2.04.snp.dist),as.matrix(MAF006.D2.04.dist), alternative = "greater")
gl.ibd(gl=NULL, Dgen=(MAF006.D2.04.dist), Dgeo=as.dist(coords.d2.04.snp.dist))
mantel(coords.d2.04.snp.dist, MAF006.D2.04.dist, method="p")


#need to reorder
d2.14.gen<-as.matrix(MAF006.D2.14.dist)
d2.14.geo<-as.matrix(coords.d2.14.snp.dist)

d2.14.gen <- d2.14.gen[, order(as.integer(colnames(d2.14.gen)))]
d2.14.gen <- d2.14.gen[order(as.integer(colnames(d2.14.gen))),]

d2.14.geo <- d2.14.geo[, order(as.integer(colnames(d2.14.geo)))]
d2.14.geo <- d2.14.geo[order(as.integer(colnames(d2.14.geo))),]

mantel.test((d2.14.geo),(d2.14.gen), alternative = "greater")
gl.ibd(gl=NULL, Dgen=as.dist(d2.14.gen), Dgeo=as.dist(d2.14.geo))
mantel(d2.14.geo, d2.14.gen, method="p")


mantel.test(as.matrix(coords.d2.15.snp.dist),as.matrix(MAF006.D2.15.dist), alternative = "greater")
gl.ibd(gl=NULL, Dgen=(MAF006.D2.15.dist), Dgeo=as.dist(coords.d2.15.snp.dist))
mantel(coords.d2.15.snp.dist, MAF006.D2.15.dist, method="p")


mantel.test(as.matrix(coords.d3.04.snp.dist),as.matrix(MAF006.D3.04.dist), alternative = "greater")
gl.ibd(gl=NULL, Dgen=(MAF006.D3.04.dist), Dgeo=as.dist(coords.d3.04.snp.dist))
mantel(coords.d3.04.snp.dist, MAF006.D3.04.dist, method="p")


#need to reorder
d3.14.gen<-as.matrix(MAF006.D3.14.dist)
d3.14.geo<-as.matrix(coords.d3.14.snp.dist)

d3.14.gen <- d3.14.gen[, order(as.integer(colnames(d3.14.gen)))]
d3.14.gen <- d3.14.gen[order(as.integer(colnames(d3.14.gen))),]

d3.14.geo <- d3.14.geo[, order(as.integer(colnames(d3.14.geo)))]
d3.14.geo <- d3.14.geo[order(as.integer(colnames(d3.14.geo))),]

mantel.test((d3.14.geo),(d3.14.gen), alternative = "greater", nperm=100000)
gl.ibd(gl=NULL, Dgen=as.dist(d3.14.gen), Dgeo=as.dist(d3.14.geo), permutations=100000)
mantel(d3.14.geo, d3.14.gen, method="p", permutations = 100000)


mantel.test(as.matrix(coords.d3.15.snp.dist),as.matrix(MAF006.D3.15.dist), alternative = "greater")
gl.ibd(gl=NULL, Dgen=(MAF006.D3.15.dist), Dgeo=as.dist(coords.d3.15.snp.dist))
mantel(coords.d3.15.snp.dist, MAF006.D3.15.dist, method="p")


mantel.test(as.matrix(coords.mira.snp.dist),as.matrix(MAF006.mira.dist), alternative = "greater", nperm = 10000)
gl.ibd(gl=NULL, Dgen=(MAF006.mira.dist), Dgeo=as.dist(coords.mira.snp.dist), permutations = 10000)
mantel(coords.mira.snp.dist, MAF006.mira.dist, method="p", permutations = 10000)


mantel.test(as.matrix(coords.agraria.snp.dist),as.matrix(MAF006.agraria.dist), alternative = "greater")
gl.ibd(gl=NULL, Dgen=(MAF006.agraria.dist), Dgeo=as.dist(coords.agraria.snp.dist))
mantel(coords.agraria.snp.dist, MAF006.agraria.dist, method="p")


#Make figure of d2.2.2014, D3.2004, JL, HD2, CESAC2006
library(gridExtra)
library(dartR)
library(reshape2)
library(ggplot2)
#pdf("Rplot.pdf",width = 4, height = 20)
#png("Rplot.png", width = 200, height=800, units = "px", pointsize = 20)

#par(mfrow=c(5,1))

S<-gl.ibd(gl=NULL, Dgen=(hier.S.sub.dist), Dgeo=as.dist(coords.S.aflp.dist))
SDgeo<-as.matrix(S$Dgeo)
SDgen<-as.matrix(S$Dgen)

SDgeo[upper.tri(SDgeo, diag=T)]<-NA
SDgen[upper.tri(SDgen, diag=T)]<-NA


Sdf<-cbind(melt(SDgeo/100),melt(SDgen))
Sdf<-na.omit(Sdf)

colnames(Sdf)[c(3,6)]<-c("geo","gen")
summary(lm(Sdf$gen~Sdf$geo))
cor.test(Sdf$gen, Sdf$geo, method="p", alternative = "g")

plS<-ggplot(Sdf[,c(3,6)], aes(x=Sdf$geo, y=Sdf$gen))+
  #stat_density2d(aes(alpha=..level.., fill=..level..), size=2, 
  #bins=8, geom="polygon") + 
  #scale_fill_gradient(low = "yellow", high = "red") +
  #scale_alpha(range = c(0, 1), guide = FALSE) +
  #geom_density2d(colour="white", bins=20) +
  geom_point(data = Sdf[,c(3,6)]) +
  stat_smooth(method="lm", level=.95, size=1, color="black")+
  #lims(x=c(-0.5,5), y=c(-1,4))+
  theme(axis.text.y   = element_text(size=10),
        axis.text.x   = element_text(size=10),
        axis.title.y  = element_text(size=12),
        axis.title.x  = element_text(size=12),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = "none",
  ) + labs(x="Distance (m)", y="Genetic Distance")



A<-gl.ibd(gl=NULL, Dgen=(hier.CESAC06.dist), Dgeo=as.dist(coords.CESAC06.aflp.dist))
ADgeo<-as.matrix(A$Dgeo)
ADgen<-as.matrix(A$Dgen)

ADgeo[upper.tri(ADgeo, diag=T)]<-NA
ADgen[upper.tri(ADgen, diag=T)]<-NA


Adf<-cbind(melt(ADgeo),melt(ADgen))
Adf<-na.omit(Adf)
colnames(Adf)[c(3,6)]<-c("geo","gen")
summary(lm(Adf$gen~Adf$geo))
cor.test(Adf$gen,Adf$geo,method="p")

plCE<-ggplot(Adf[,c(3,6)], aes(x=Adf$geo, y=Adf$gen))+
  #stat_density2d(aes(alpha=..level.., fill=..level..), size=2, 
  #               bins=100, geom="polygon") + 
  #scale_fill_gradient(low = "yellow", high = "red") +
  #scale_alpha(range = c(0, 1), guide = FALSE) +
  geom_density2d(colour="white", bins=20) +
  geom_point(data = Adf[,c(3,6)]) +
  stat_smooth(method="lm", level=.95, size=1, color="black")+
  #lims(x=c(min(Adf$geo),max(Adf$geo)), y=c(min(Adf$gen),max(Adf$gen)))+
  theme(axis.text.y   = element_text(size=10),
        axis.text.x   = element_text(size=10),
        axis.title.y  = element_text(size=12),
        axis.title.x  = element_text(size=12),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = "none",
  ) + labs(x="Distance (m)", y="Genetic Distance")



B<-gl.ibd(gl=NULL, Dgen=(hier.JL.dist), Dgeo=as.dist(coords.JL.aflp.dist))
BDgeo<-as.matrix(B$Dgeo)
BDgen<-as.matrix(B$Dgen)

BDgeo[upper.tri(BDgeo, diag=T)]<-NA
BDgen[upper.tri(BDgen, diag=T)]<-NA


Bdf<-cbind(melt(BDgeo),melt(BDgen))
Bdf<-na.omit(Bdf)
colnames(Bdf)[c(3,6)]<-c("geo","gen")
summary(lm(Bdf$gen~Bdf$geo))
cor.test(Bdf$gen, Bdf$geo, method="p")

plJL<-ggplot(Bdf[,c(3,6)], aes(x=Bdf$geo, y=Bdf$gen))+
  # stat_density2d(aes(alpha=..level.., fill=..level..), size=2, 
  #                bins=100, geom="polygon") + 
  # scale_fill_gradient(low = "yellow", high = "red") +
  # scale_alpha(range = c(0, 1), guide = FALSE) +
  # geom_density2d(colour="white", bins=20) +
  geom_point(data = Bdf[,c(3,6)]) +
  stat_smooth(method="lm", level=.95, size=1, color="black")+
  #lims(x=c(min(Bdf$geo),max(Bdf$geo)), y=c(min(Bdf$gen),max(Bdf$gen)))+
  theme(axis.text.y   = element_text(size=10),
        axis.text.x   = element_text(size=10),
        axis.title.y  = element_text(size=12),
        axis.title.x  = element_text(size=12),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = "none",
  ) + labs(x="Distance (m)", y="Genetic Distance")


E<-gl.ibd(gl=NULL, Dgen=as.dist(d2.14.gen), Dgeo=as.dist(d2.14.geo))
EDgeo<-as.matrix(E$Dgeo)
EDgen<-as.matrix(E$Dgen)

EDgeo[upper.tri(EDgeo, diag=T)]<-NA
EDgen[upper.tri(EDgen, diag=T)]<-NA


Edf<-cbind(melt(EDgeo),melt(EDgen))
Edf<-na.omit(Edf)
colnames(Edf)[c(3,6)]<-c("geo","gen")
summary(lm(Edf$gen~Edf$geo))
cor.test( Edf$geo,Edf$gen, method="p")

pldr214<-ggplot(Edf[,c(3,6)], aes(x=Edf$geo, y=Edf$gen))+
  # stat_density2d(aes(alpha=..level.., fill=..level..), size=2, 
  #                bins=100, geom="polygon") + 
  # scale_fill_gradient(low = "yellow", high = "red") +
  # scale_alpha(range = c(0, 1), guide = FALSE) +
  # geom_density2d(colour="white", bins=20) +
  geom_point(data = Edf[,c(3,6)]) +
  stat_smooth(method="lm", level=.95, size=1, color="black")+
  #lims(x=c(min(Edf$geo),max(Edf$geo)), y=c(min(Edf$gen),max(Edf$gen)))+
  theme(axis.text.y   = element_text(size=10),
        axis.text.x   = element_text(size=10),
        axis.title.y  = element_text(size=12),
        axis.title.x  = element_text(size=12),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = "none",
  ) + labs(x="Distance (m)", y="Genetic Distance")
#+scale_y_continuous(labels = function(x) round(x/100))

library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
r <- rectGrob(gp=gpar(fill="white"))
grid.arrange(plS, r, plCE, r, plJL,r, pldr214, r, ncol=2, nrow=4)
grid.arrange(plS, plCE, plJL, pldr214, ncol=2, nrow=2)
#################################################################
#################################################################
#All for suppl
library(reshape2)
library(vegan)
scaleFUN <- function(x) sprintf("%.2f", x)

#serbia
S<-gl.ibd(gl=NULL, Dgen=(hier.S.sub.dist), Dgeo=as.dist(coords.S.aflp.dist))
SDgeo<-as.matrix(S$Dgeo)
SDgen<-as.matrix(S$Dgen)

SDgeo[upper.tri(SDgeo, diag=T)]<-NA
SDgen[upper.tri(SDgen, diag=T)]<-NA


Sdf<-cbind(melt(SDgeo/100),melt(SDgen))
Sdf<-na.omit(Sdf)

colnames(Sdf)[c(3,6)]<-c("geo","gen")
summary(lm(Sdf$gen~Sdf$geo))
cor.test(y=Sdf$gen, x=Sdf$geo, method="pearson", alternative = "g")


plS<-ggplot(Sdf[,c(3,6)], aes(x=Sdf$geo, y=Sdf$gen))+
  #stat_density2d(aes(alpha=..level.., fill=..level..), size=2, 
  #              bins=8, geom="polygon") + 
  scale_fill_gradient(low = "yellow", high = "red") +
  scale_alpha(range = c(0, 1), guide = FALSE) +
  #geom_density2d(colour="white", bins=20) +
  geom_point(data = Sdf[,c(3,6)]) +
  stat_smooth(method="lm", level=.95, size=1, color="black")+
  #lims(x=c(-0.5,5), y=c(-1,4))+
  theme(axis.text.y   = element_text(size=10, angle = 90, hjust = .5),
        axis.text.x   = element_text(size=10),
        axis.title.y  = element_text(size=12),
        axis.title.x  = element_text(size=12),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = "none",
  ) + labs(x="Distance (m)", y="Genetic Distance")+coord_cartesian(xlim=c(min(Sdf$geo),max(Sdf$geo)), ylim=c(min(Sdf$gen),max(Sdf$gen)))+
  scale_y_continuous(labels=scaleFUN)+scale_x_continuous(labels=scaleFUN)+ggtitle("Serbia 2007")



# cesac02
A02<-gl.ibd(gl=NULL, Dgen=(hier.CESAC02.dist), Dgeo=as.dist(coords.CESAC02.aflp.dist))
A02Dgeo<-as.matrix(A02$Dgeo)
A02Dgen<-as.matrix(A02$Dgen)

A02Dgeo[upper.tri(A02Dgeo, diag=T)]<-NA
A02Dgen[upper.tri(A02Dgen, diag=T)]<-NA


A02df<-cbind(melt(A02Dgeo),melt(A02Dgen))
A02df<-na.omit(A02df)
colnames(A02df)[c(3,6)]<-c("geo","gen")
summary(lm(A02df$gen~A02df$geo))
cor.test(y=A02df$gen, x=A02df$geo, method="pearson", alternative = "g")

plCE02<-ggplot(A02df[,c(3,6)], aes(x=A02df$geo, y=A02df$gen))+
  #stat_density2d(aes(alpha=..level.., fill=..level..), size=2, 
  #            bins=100, geom="polygon") + 
  scale_fill_gradient(low = "yellow", high = "red") +
  scale_alpha(range = c(0, 1), guide = FALSE) +
  #geom_density2d(colour="white", bins=20) +
  geom_point(data = A02df[,c(3,6)]) +
  stat_smooth(method="lm", level=.95, size=1, color="black")+
  #lims(x=c(min(Adf$geo),max(A02df$geo)), y=c(min(A02df$gen),max(A02df$gen)))+
  theme(axis.text.y   = element_text(size=10, angle = 90, hjust = .5),
        axis.text.x   = element_text(size=10),
        axis.title.y  = element_text(size=12),
        axis.title.x  = element_text(size=12),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = "none",
  ) + labs(x="Distance (m)", y="Genetic Distance")+coord_cartesian(xlim=c(min(A02df$geo),max(A02df$geo)), ylim=c(min(A02df$gen),max(A02df$gen)))+
  scale_y_continuous(labels=scaleFUN)+scale_x_continuous(labels=scaleFUN)+ggtitle("CESAC 2002")



#cesac 06
A<-gl.ibd(gl=NULL, Dgen=(hier.CESAC06.dist), Dgeo=as.dist(coords.CESAC06.aflp.dist))
ADgeo<-as.matrix(A$Dgeo)
ADgen<-as.matrix(A$Dgen)

ADgeo[upper.tri(ADgeo, diag=T)]<-NA
ADgen[upper.tri(ADgen, diag=T)]<-NA


Adf<-cbind(melt(ADgeo),melt(ADgen))
Adf<-na.omit(Adf)
colnames(Adf)[c(3,6)]<-c("geo","gen")
summary(lm(Adf$gen~Adf$geo))
cor.test(y=Adf$gen, x=Adf$geo, method="p", alternative = "g")

plCE<-ggplot(Adf[,c(3,6)], aes(x=Adf$geo, y=Adf$gen))+
  #stat_density2d(aes(alpha=..level.., fill=..level..), size=2, 
  #             bins=100, geom="polygon") + 
  scale_fill_gradient(low = "yellow", high = "red") +
  scale_alpha(range = c(0, 1), guide = FALSE) +
  #geom_density2d(colour="white", bins=20) +
  geom_point(data = Adf[,c(3,6)]) +
  stat_smooth(method="lm", level=.95, size=1, color="black")+
  #lims(x=c(min(Adf$geo),max(Adf$geo)), y=c(min(Adf$gen),max(Adf$gen)))+
  theme(axis.text.y   = element_text(size=10, angle = 90, hjust = .5),
        axis.text.x   = element_text(size=10),
        axis.title.y  = element_text(size=12),
        axis.title.x  = element_text(size=12),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = "none",
  ) + labs(x="Distance (m)", y="Genetic Distance")+coord_cartesian(xlim=c(min(Adf$geo),max(Adf$geo)), ylim=c(min(Adf$gen),max(Adf$gen)))+
  scale_y_continuous(labels=scaleFUN)+scale_x_continuous(labels=scaleFUN)+ggtitle("CESAC 2006")



#jl
B<-gl.ibd(gl=NULL, Dgen=(hier.JL.dist), Dgeo=as.dist(coords.JL.aflp.dist))
BDgeo<-as.matrix(B$Dgeo)
BDgen<-as.matrix(B$Dgen)

BDgeo[upper.tri(BDgeo, diag=T)]<-NA
BDgen[upper.tri(BDgen, diag=T)]<-NA


Bdf<-cbind(melt(BDgeo),melt(BDgen))
Bdf<-na.omit(Bdf)
colnames(Bdf)[c(3,6)]<-c("geo","gen")
summary(lm(Bdf$gen~Bdf$geo))
cor.test(y=Bdf$gen, x=Bdf$geo,method="pearson", alternative = "g")

plJL<-ggplot(Bdf[,c(3,6)], aes(x=Bdf$geo, y=Bdf$gen))+
  #stat_density2d(aes(alpha=..level.., fill=..level..), size=2, 
  #              bins=100, geom="polygon") + 
  scale_fill_gradient(low = "yellow", high = "red") +
  scale_alpha(range = c(0, 1), guide = FALSE) +
  #geom_density2d(colour="white", bins=20) +
  geom_point(data = Bdf[,c(3,6)]) +
  stat_smooth(method="lm", level=.95, size=1, color="black")+
  #lims(x=c(min(Bdf$geo),max(Bdf$geo)), y=c(min(Bdf$gen),max(Bdf$gen)))+
  theme(axis.text.y   = element_text(size=10, angle = 90, hjust = .5),
        axis.text.x   = element_text(size=10),
        axis.title.y  = element_text(size=12),
        axis.title.x  = element_text(size=12),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = "none",
  ) + labs(x="Distance (m)", y="Genetic Distance")+coord_cartesian(xlim=c(min(Bdf$geo),max(Bdf$geo)), ylim=c(min(Bdf$gen),max(Bdf$gen)))+
  scale_y_continuous(labels=scaleFUN)+scale_x_continuous(labels=scaleFUN)+ggtitle("Jake's Landing 2006")


#str
str<-gl.ibd(gl=NULL, Dgen=(hier.str.dist), Dgeo=as.dist(coords.str.aflp.dist))
strDgeo<-as.matrix(str$Dgeo)
strDgen<-as.matrix(str$Dgen)

strDgeo[upper.tri(strDgeo, diag=T)]<-NA
strDgen[upper.tri(strDgen, diag=T)]<-NA


strdf<-cbind(melt(strDgeo),melt(strDgen))
strdf<-na.omit(strdf)
colnames(strdf)[c(3,6)]<-c("geo","gen")
summary(lm(strdf$gen~strdf$geo))
cor.test(y=strdf$gen, x=strdf$geo, method="p", alternative = "g")

plstr<-ggplot(strdf[,c(3,6)], aes(x=strdf$geo, y=strdf$gen))+
  #stat_density2d(aes(alpha=..level.., fill=..level..), size=2, 
  #          bins=100, geom="polygon") + 
  scale_fill_gradient(low = "yellow", high = "red") +
  scale_alpha(range = c(0, 1), guide = FALSE) +
  #geom_density2d(colour="white", bins=20) +
  geom_point(data = strdf[,c(3,6)]) +
  stat_smooth(method="lm", level=.95, size=1, color="black")+
  #lims(x=c(min(strdf$geo),max(strdf$geo)), y=c(min(strdf$gen),max(strdf$gen)))+
  theme(axis.text.y   = element_text(size=10, angle = 90, hjust = .5),
        axis.text.x   = element_text(size=10),
        axis.title.y  = element_text(size=12),
        axis.title.x  = element_text(size=12),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = "none",
  ) + labs(x="Distance (m)", y="Genetic Distance")+coord_cartesian(xlim=c(min(strdf$geo),max(strdf$geo)), ylim=c(min(strdf$gen),max(strdf$gen)))+
  scale_y_continuous(labels=scaleFUN)+scale_x_continuous(labels=scaleFUN)+ggtitle("Round Valley 2006")


#roch1
roch1<-gl.ibd(gl=NULL, Dgen=(hier.Roch1.dist), Dgeo=as.dist(coords.roch1.aflp.dist))
roch1Dgeo<-as.matrix(roch1$Dgeo)
roch1Dgen<-as.matrix(roch1$Dgen)

roch1Dgeo[upper.tri(roch1Dgeo, diag=T)]<-NA
roch1Dgen[upper.tri(roch1Dgen, diag=T)]<-NA


roch1df<-cbind(melt(roch1Dgeo),melt(roch1Dgen))
roch1df<-na.omit(roch1df)
colnames(roch1df)[c(3,6)]<-c("geo","gen")
summary(lm(roch1df$gen~roch1df$geo))
cor.test(y=roch1df$gen,x=roch1df$geo, method="p", alternative = "g")

plroch1<-ggplot(roch1df[,c(3,6)], aes(x=roch1df$geo, y=roch1df$gen))+
  #stat_density2d(aes(alpha=..level.., fill=..level..), size=2, 
  #         bins=100, geom="polygon") + 
  scale_fill_gradient(low = "yellow", high = "red") +
  scale_alpha(range = c(0, 1), guide = FALSE) +
  #geom_density2d(colour="white", bins=20) +
  geom_point(data = roch1df[,c(3,6)]) +
  stat_smooth(method="lm", level=.95, size=1, color="black")+
  # lims(x=c(min(roch1df$geo),max(roch1df$geo)), y=c(min(roch1df$gen),max(roch1df$gen)))+
  theme(axis.text.y   = element_text(size=10, angle = 90, hjust = .5),
        axis.text.x   = element_text(size=10),
        axis.title.y  = element_text(size=12),
        axis.title.x  = element_text(size=12),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = "none",
  ) + labs(x="Distance (m)", y="Genetic Distance")+coord_cartesian(xlim=c(min(roch1df$geo),max(roch1df$geo)), ylim=c(min(roch1df$gen),max(roch1df$gen)))+
  scale_y_continuous(labels=scaleFUN)+scale_x_continuous(labels=scaleFUN)+ggtitle("Rochester 1 2007")


#roch 2
roch2<-gl.ibd(gl=NULL, Dgen=(hier.Roch2.sub.dist), Dgeo=as.dist(coords.roch2.aflp.dist))
roch2Dgeo<-as.matrix(roch2$Dgeo)
roch2Dgen<-as.matrix(roch2$Dgen)

roch2Dgeo[upper.tri(roch2Dgeo, diag=T)]<-NA
roch2Dgen[upper.tri(roch2Dgen, diag=T)]<-NA


roch2df<-cbind(melt(roch2Dgeo),melt(roch2Dgen))
roch2df<-na.omit(roch2df)
colnames(roch2df)[c(3,6)]<-c("geo","gen")
summary(lm(roch2df$gen~roch2df$geo))
cor.test(y=roch2df$gen,x=roch2df$geo, method="p", alternative = "g")

plroch2<-ggplot(roch2df[,c(3,6)], aes(x=roch2df$geo, y=roch2df$gen))+
  #stat_density2d(aes(alpha=..level.., fill=..level..), size=2, 
  #               bins=100, geom="polygon") + 
  scale_fill_gradient(low = "yellow", high = "red") +
  scale_alpha(range = c(0, 1), guide = FALSE) +
  #geom_density2d(colour="white", bins=20) +
  geom_point(data = roch2df[,c(3,6)]) +
  stat_smooth(method="lm", level=.95, size=1, color="black", fullrange = T)+
  #lims(x=c(min(roch2df$geo),max(roch2df$geo)), y=c(min(roch2df$gen),max(roch2df$gen)))+
  theme(axis.text.y   = element_text(size=10, angle = 90, hjust = .5),
        axis.text.x   = element_text(size=10),
        axis.title.y  = element_text(size=12),
        axis.title.x  = element_text(size=12),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = "none",
  ) + labs(x="Distance (m)", y="Genetic Distance")+coord_cartesian(xlim=c(min(roch2df$geo),max(roch2df$geo)), ylim=c(min(roch2df$gen),max(roch2df$gen)))+
  scale_y_continuous(labels=scaleFUN)+scale_x_continuous(labels=scaleFUN)+ggtitle("Rcohester 2 2007")


#hd1
hd1<-gl.ibd(gl=NULL, Dgen=(hier.hd1.dist), Dgeo=as.dist(coords.hd1.aflp.dist))
hd1Dgeo<-as.matrix(hd1$Dgeo)
hd1Dgen<-as.matrix(hd1$Dgen)

hd1Dgeo[upper.tri(hd1Dgeo, diag=T)]<-NA
hd1Dgen[upper.tri(hd1Dgen, diag=T)]<-NA


hd1df<-cbind(melt(hd1Dgeo),melt(hd1Dgen))
hd1df<-na.omit(hd1df)
colnames(hd1df)[c(3,6)]<-c("geo","gen")
summary(lm(hd1df$gen~hd1df$geo))
cor.test(y=hd1df$gen,x=hd1df$geo,method="p", alternative = "g")

plhd1<-ggplot(hd1df[,c(3,6)], aes(x=hd1df$geo, y=hd1df$gen))+
  #stat_density2d(aes(alpha=..level.., fill=..level..), size=2, 
  #           bins=100, geom="polygon") + 
  scale_fill_gradient(low = "yellow", high = "red") +
  scale_alpha(range = c(0, 1), guide = FALSE) +
  #geom_density2d(colour="white", bins=20) +
  geom_point(data = hd1df[,c(3,6)]) +
  stat_smooth(method="lm", level=.95, size=1, color="black")+
  #lims(x=c(min(hd1df$geo),max(hd1df$geo)), y=c(min(hd1df$gen),max(hd1df$gen)))+
  theme(axis.text.y   = element_text(size=10, angle = 90, hjust = .5),
        axis.text.x   = element_text(size=10),
        axis.title.y  = element_text(size=12),
        axis.title.x  = element_text(size=12),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = "none",
  ) + labs(x="Distance (m)", y="Genetic Distance")+coord_cartesian(xlim=c(min(hd1df$geo),max(hd1df$geo)), ylim=c(min(hd1df$gen),max(hd1df$gen)))+
  scale_y_continuous(labels=scaleFUN)+scale_x_continuous(labels=scaleFUN)+ggtitle("Heart's Desire 1 2004")



#hd2
hd2<-gl.ibd(gl=NULL, Dgen=as.dist(hd2.gen), Dgeo=as.dist(hd2.geo))
hd2Dgeo<-as.matrix(hd2$Dgeo)
hd2Dgen<-as.matrix(hd2$Dgen)

hd2Dgeo[upper.tri(hd2Dgeo, diag=T)]<-NA
hd2Dgen[upper.tri(hd2Dgen, diag=T)]<-NA


hd2df<-cbind(melt(hd2Dgeo),melt(hd2Dgen))
hd2df<-na.omit(hd2df)
colnames(hd2df)[c(3,6)]<-c("geo","gen")
summary(lm(hd2df$gen~hd2df$geo))
cor.test(y=hd2df$gen, x=hd2df$geo, method="p", alternative = 'g')

plhd2<-ggplot(hd2df[,c(3,6)], aes(x=hd2df$geo, y=hd2df$gen))+
  #stat_density2d(aes(alpha=..level.., fill=..level..), size=2, 
  #               bins=100, geom="polygon") + 
  scale_fill_gradient(low = "yellow", high = "red") +
  scale_alpha(range = c(0, 1), guide = FALSE) +
  #geom_density2d(colour="white", bins=20) +
  geom_point(data = hd2df[,c(3,6)]) +
  stat_smooth(method="lm", level=.95, size=1, color="black")+
  #lims(x=c(min(hd2df$geo),max(hd2df$geo)), y=c(min(hd2df$gen),max(hd2df$gen)))+
  theme(axis.text.y   = element_text(size=10, angle = 90, hjust = .5),
        axis.text.x   = element_text(size=10),
        axis.title.y  = element_text(size=12),
        axis.title.x  = element_text(size=12),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = "none",
  ) + labs(x="Distance (m)", y="Genetic Distance")+coord_cartesian(xlim=c(min(hd2df$geo),max(hd2df$geo)), ylim=c(min(hd2df$gen),max(hd2df$gen)))+
  scale_y_continuous(labels=scaleFUN)+scale_x_continuous(labels=scaleFUN)+ggtitle("Heart's Desire 2 2004")


#hd3
hd3<-gl.ibd(gl=NULL, Dgen=(hier.hd3.dist), Dgeo=as.dist(coords.hd3.aflp.dist))
hd3Dgeo<-as.matrix(hd3$Dgeo)
hd3Dgen<-as.matrix(hd3$Dgen)

hd3Dgeo[upper.tri(hd3Dgeo, diag=T)]<-NA
hd3Dgen[upper.tri(hd3Dgen, diag=T)]<-NA


hd3df<-cbind(melt(hd3Dgeo),melt(hd3Dgen))
hd3df<-na.omit(hd3df)
colnames(hd3df)[c(3,6)]<-c("geo","gen")
summary(lm(hd3df$gen~hd3df$geo))
cor.test(y=hd3df$gen, x=hd3df$geo, method="p", alternative = "g")

plhd3<-ggplot(hd3df[,c(3,6)], aes(x=hd3df$geo, y=hd3df$gen))+
  #stat_density2d(aes(alpha=..level.., fill=..level..), size=2, 
  #               bins=100, geom="polygon") + 
  scale_fill_gradient(low = "yellow", high = "red") +
  scale_alpha(range = c(0, 1), guide = FALSE) +
  #geom_density2d(colour="white", bins=20) +
  geom_point(data = hd3df[,c(3,6)]) +
  stat_smooth(method="lm", level=.95, size=1, color="black")+
  #coord_cartesian(xlim=c(min(hd3df$geo),max(hd3df$geo)), ylim=c(min(hd3df$gen),max(hd3df$gen)))+
  theme(axis.text.y   = element_text(size=10, angle = 90, hjust = .5),
        axis.text.x   = element_text(size=10),
        axis.title.y  = element_text(size=12),
        axis.title.x  = element_text(size=12),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = "none",
  ) + labs(x="Distance (m)", y="Genetic Distance")+coord_cartesian(xlim=c(min(hd3df$geo),max(hd3df$geo)), ylim=c(min(hd3df$gen),max(hd3df$gen)))+
  scale_y_continuous(labels=scaleFUN)+scale_x_continuous(labels=scaleFUN)+ggtitle("Heart's Desire 3 2004")




#d1
d1<-gl.ibd(gl=NULL, Dgen=(hier.d1.dist), Dgeo=as.dist(coords.d1.aflp.dist))
d1Dgeo<-as.matrix(d1$Dgeo)
d1Dgen<-as.matrix(d1$Dgen)

d1Dgeo[upper.tri(d1Dgeo, diag=T)]<-NA
d1Dgen[upper.tri(d1Dgen, diag=T)]<-NA


d1df<-cbind(melt(d1Dgeo),melt(d1Dgen))
d1df<-na.omit(d1df)
colnames(d1df)[c(3,6)]<-c("geo","gen")
summary(lm(d1df$gen~d1df$geo))
cor.test(y=d1df$gen, x=d1df$geo, method="p", alternative = 'g')

pld1<-ggplot(d1df[,c(3,6)], aes(x=d1df$geo, y=d1df$gen))+
  #stat_density2d(aes(alpha=..level.., fill=..level..), size=2, 
  #              bins=100, geom="polygon") + 
  scale_fill_gradient(low = "yellow", high = "red") +
  scale_alpha(range = c(0, 1), guide = FALSE) +
  #geom_density2d(colour="white", bins=20) +
  geom_point(data = d1df[,c(3,6)]) +
  stat_smooth(method="lm", level=.95, size=1, color="black")+
  #lims(x=c(min(d1df$geo),max(d1df$geo)), y=c(min(d1df$gen),max(d1df$gen)))+
  theme(axis.text.y   = element_text(size=10, angle = 90, hjust = .5),
        axis.text.x   = element_text(size=10),
        axis.title.y  = element_text(size=12),
        axis.title.x  = element_text(size=12),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = "none",
  ) + labs(x="Distance (m)", y="Genetic Distance")+coord_cartesian(xlim=c(min(d1df$geo),max(d1df$geo)), ylim=c(min(d1df$gen),max(d1df$gen)))+
  scale_y_continuous(labels=scaleFUN)+scale_x_continuous(labels=scaleFUN)+ggtitle("Drake 1 2004")



#d2
d2<-gl.ibd(gl=NULL, Dgen=(hier.d2.dist), Dgeo=as.dist(coords.d2.04.aflp.dist))
d2Dgeo<-as.matrix(d2$Dgeo)
d2Dgen<-as.matrix(d2$Dgen)

d2Dgeo[upper.tri(d2Dgeo, diag=T)]<-NA
d2Dgen[upper.tri(d2Dgen, diag=T)]<-NA


d2df<-cbind(melt(d2Dgeo),melt(d2Dgen))
d2df<-na.omit(d2df)
colnames(d2df)[c(3,6)]<-c("geo","gen")
summary(lm(d2df$gen~d2df$geo))
cor.test(y=d2df$gen, x=d2df$geo, method="p", alternative = 'g')

pld2<-ggplot(d2df[,c(3,6)], aes(x=d2df$geo, y=d2df$gen))+
  #stat_density2d(aes(alpha=..level.., fill=..level..), size=2, 
  #              bins=100, geom="polygon") + 
  scale_fill_gradient(low = "yellow", high = "red") +
  scale_alpha(range = c(0, 1), guide = FALSE) +
  #geom_density2d(colour="white", bins=20) +
  geom_point(data = d2df[,c(3,6)]) +
  stat_smooth(method="lm", level=.95, size=1, color="black")+
  #lims(x=c(min(d2df$geo),max(d2df$geo)), y=c(min(d2df$gen),max(d2df$gen)))+
  theme(axis.text.y   = element_text(size=10, angle = 90, hjust = .5),
        axis.text.x   = element_text(size=10),
        axis.title.y  = element_text(size=12),
        axis.title.x  = element_text(size=12),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = "none",
  ) + labs(x="Distance (m)", y="Genetic Distance")+coord_cartesian(xlim=c(min(d2df$geo),max(d2df$geo)), ylim=c(min(d2df$gen),max(d2df$gen)))+
  scale_y_continuous(labels=scaleFUN)+scale_x_continuous(labels=scaleFUN)+ggtitle("Drake 2 2004")


#d3
d3<-gl.ibd(gl=NULL, Dgen=(hier.d3.dist), Dgeo=as.dist(coords.d3.04.aflp.dist))
d3Dgeo<-as.matrix(d3$Dgeo)
d3Dgen<-as.matrix(d3$Dgen)

d3Dgeo[upper.tri(d3Dgeo, diag=T)]<-NA
d3Dgen[upper.tri(d3Dgen, diag=T)]<-NA


d3df<-cbind(melt(d3Dgeo),melt(d3Dgen))
d3df<-na.omit(d3df)
colnames(d3df)[c(3,6)]<-c("geo","gen")
summary(lm(d3df$gen~d3df$geo))
cor.test(y=d3df$gen,x=d3df$geo, method="p", alternative = 'g')

pld3<-ggplot(d3df[,c(3,6)], aes(x=d3df$geo, y=d3df$gen))+
  #stat_density2d(aes(alpha=..level.., fill=..level..), size=2, 
  #          bins=100, geom="polygon") + 
  scale_fill_gradient(low = "yellow", high = "red") +
  scale_alpha(range = c(0, 1), guide = FALSE) +
  #geom_density2d(colour="white", bins=20) +
  geom_point(data = d3df[,c(3,6)]) +
  stat_smooth(method="lm", level=.95, size=1, color="black")+
  #lims(x=c(min(d3df$geo),max(d3df$geo)), y=c(min(d3df$gen),max(d3df$gen)))+
  theme(axis.text.y   = element_text(size=10, angle = 90, hjust = .5),
        axis.text.x   = element_text(size=10),
        axis.title.y  = element_text(size=12),
        axis.title.x  = element_text(size=12),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = "none",
  ) + labs(x="Distance (m)", y="Genetic Distance")+coord_cartesian(xlim=c(min(d3df$geo),max(d3df$geo)), ylim=c(min(d3df$gen),max(d3df$gen)))+
  scale_y_continuous(labels=scaleFUN)+scale_x_continuous(labels=scaleFUN)+ggtitle("Drake 3 2004")




#d4
d4<-gl.ibd(gl=NULL, Dgen=(hier.d4.dist), Dgeo=as.dist(coords.d4.aflp.dist))
d4Dgeo<-as.matrix(d4$Dgeo)
d4Dgen<-as.matrix(d4$Dgen)

d4Dgeo[upper.tri(d4Dgeo, diag=T)]<-NA
d4Dgen[upper.tri(d4Dgen, diag=T)]<-NA


d4df<-cbind(melt(d4Dgeo),melt(d4Dgen))
d4df<-na.omit(d4df)
colnames(d4df)[c(3,6)]<-c("geo","gen")
summary(lm(d4df$gen~d4df$geo))
cor.test(y=d4df$gen,x=d4df$geo, method="p", alternative = "g")

pld4<-ggplot(d4df[,c(3,6)], aes(x=d4df$geo, y=d4df$gen))+
  #stat_density2d(aes(alpha=..level.., fill=..level..), size=2, 
  #             bins=100, geom="polygon") + 
  scale_fill_gradient(low = "yellow", high = "red") +
  scale_alpha(range = c(0, 1), guide = FALSE) +
  #geom_density2d(colour="white", bins=20) +
  geom_point(data = d4df[,c(3,6)]) +
  stat_smooth(method="lm", level=.95, size=1, color="black")+
  #lims(x=c(min(d4df$geo),max(d4df$geo)), y=c(min(d4df$gen),max(d4df$gen)))+
  theme(axis.text.y   = element_text(size=10, angle = 90, hjust = .5),
        axis.text.x   = element_text(size=10),
        axis.title.y  = element_text(size=12),
        axis.title.x  = element_text(size=12),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = "none",
  ) + labs(x="Distance (m)", y="Genetic Distance")+coord_cartesian(xlim=c(min(d4df$geo),max(d4df$geo)), ylim=c(min(d4df$gen),max(d4df$gen)))+
  scale_y_continuous(labels=scaleFUN)+scale_x_continuous(labels=scaleFUN)+ggtitle("Drake 4 2004")

library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
r <- rectGrob(gp=gpar(fill="white"))
#grid.arrange(plS, plCE02, plCE, plJL, plstr, r, plroch1, plroch2, r, plhd1, plhd2, plhd3, pld1, pld2, r, pld3, pld4, r,  ncol=3, nrow=6)

grid.arrange(plS, plCE02, plCE, plJL, plstr, r, plroch1, plroch2, r,  ncol=3, nrow=3)
             
grid.arrange(plhd1, plhd2, plhd3, pld1, pld2, r, pld3, pld4, r,  ncol=3, nrow=3)


#############
#############
#SNP

#mira
mira<-gl.ibd(gl=NULL, Dgen=(MAF006.mira.dist), Dgeo=as.dist(coords.mira.snp.dist))
miraDgeo<-as.matrix(mira$Dgeo)
miraDgen<-as.matrix(mira$Dgen)

miraDgeo[upper.tri(miraDgeo, diag=T)]<-NA
miraDgen[upper.tri(miraDgen, diag=T)]<-NA


miradf<-cbind(melt(miraDgeo),melt(miraDgen))
miradf<-na.omit(miradf)
colnames(miradf)[c(3,6)]<-c("geo","gen")
summary(lm(miradf$gen~miradf$geo))
cor.test(y=miradf$gen,x=miradf$geo, method="p", alternative = "g")

plmira<-ggplot(miradf[,c(3,6)], aes(x=miradf$geo, y=miradf$gen))+
  #stat_density2d(aes(alpha=..level.., fill=..level..), size=2, 
  #            bins=100, geom="polygon") + 
  scale_fill_gradient(low = "yellow", high = "red") +
  scale_alpha(range = c(0, 1), guide = FALSE) +
  #geom_density2d(colour="white", bins=20) +
  geom_point(data = miradf[,c(3,6)]) +
  stat_smooth(method="lm", level=.95, size=1, color="black")+
  #lims(x=c(min(miradf$geo),max(miradf$geo)), y=c(min(miradf$gen),max(miradf$gen)))+
  theme(axis.text.y   = element_text(size=10, angle = 90, hjust = .5),
        axis.text.x   = element_text(size=10),
        axis.title.y  = element_text(size=12),
        axis.title.x  = element_text(size=12),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = "none",
  ) + labs(x="Distance (m)", y="Genetic Distance")+coord_cartesian(xlim=c(min(miradf$geo),max(miradf$geo)), ylim=c(min(miradf$gen),max(miradf$gen)))+
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE))+scale_x_continuous(labels=scaleFUN)+ggtitle("Mira 2015")



#agraria
#agr<-gl.ibd(gl=NULL, Dgen=(MAF006.agraria.dist), Dgeo=as.dist(coords.agraria.snp.dist))
#agrDgeo<-as.matrix(agr$Dgeo)
#agrDgen<-as.matrix(agr$Dgen)

#agrDgeo[upper.tri(agrDgeo, diag=T)]<-NA
#agrDgen[upper.tri(agrDgen, diag=T)]<-NA


#agrdf<-cbind(melt(agrDgeo),melt(agrDgen))
#agrdf<-na.omit(agrdf)
#colnames(agrdf)[c(3,6)]<-c("geo","gen")
#summary(lm(agrdf$gen~agrdf$geo))

#plagr<-ggplot(agrdf[,c(3,6)], aes(x=agrdf$geo, y=agrdf$gen))+
#  stat_density2d(aes(alpha=..level.., fill=..level..), size=2, 
#                 bins=100, geom="polygon") + 
#  scale_fill_gradient(low = "yellow", high = "red") +
#  scale_alpha(range = c(0, 1), guide = FALSE) +
#  geom_density2d(colour="white", bins=20) +
#  geom_point(data = agrdf[,c(3,6)]) +
#  stat_smooth(method="lm", level=.95, size=1, color="black")+
#  lims(x=c(min(agrdf$geo),max(agrdf$geo)), y=c(min(agrdf$gen),max(agrdf$gen)))+
#  theme(axis.text.y   = element_text(size=10),
#        axis.text.x   = element_text(size=10),
#        axis.title.y  = element_text(size=12),
#        axis.title.x  = element_text(size=12),
#        panel.background = element_blank(),
#        panel.grid.major = element_blank(), 
#        panel.grid.minor = element_blank(),
#        axis.line = element_line(colour = "black"),
#        panel.border = element_rect(colour = "black", fill=NA, size=1),
#        legend.position = "none",
#  ) + labs(x="Distance (m)", y="Genetic Distance")


#dr2.04
dr204<-gl.ibd(gl=NULL, Dgen=(MAF006.D2.04.dist), Dgeo=as.dist(coords.d2.04.snp.dist))
dr204Dgeo<-as.matrix(dr204$Dgeo)
dr204Dgen<-as.matrix(dr204$Dgen)

dr204Dgeo[upper.tri(dr204Dgeo, diag=T)]<-NA
dr204Dgen[upper.tri(dr204Dgen, diag=T)]<-NA


dr204df<-cbind(melt(dr204Dgeo),melt(dr204Dgen))
dr204df<-na.omit(dr204df)
colnames(dr204df)[c(3,6)]<-c("geo","gen")
summary(lm(dr204df$gen~dr204df$geo))
cor.test(y=dr204df$gen, x=dr204df$geo,method="pearson", alternative = "g")

pldr204<-ggplot(dr204df[,c(3,6)], aes(x=dr204df$geo, y=dr204df$gen))+
  #stat_density2d(aes(alpha=..level.., fill=..level..), size=2, 
  #               bins=100, geom="polygon") + 
  scale_fill_gradient(low = "yellow", high = "red") +
  scale_alpha(range = c(0, 1), guide = FALSE) +
  #geom_density2d(colour="white", bins=20) +
  geom_point(data = dr204df[,c(3,6)]) +
  stat_smooth(method="lm", level=.95, size=1, color="black")+
  #lims(x=c(min(dr204df$geo),max(dr204df$geo)), y=c(min(dr204df$gen),max(dr204df$gen)))+
  theme(axis.text.y   = element_text(size=10, angle=90, hjust=0.5),
        axis.text.x   = element_text(size=10),
        axis.title.y  = element_text(size=12),
        axis.title.x  = element_text(size=12),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = "none",
  ) + labs(x="Distance (m)", y="Genetic Distance")+coord_cartesian(xlim=c(min(dr204df$geo),max(dr204df$geo)), ylim=c(min(dr204df$gen),max(dr204df$gen)))+
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE))+scale_x_continuous(labels=scaleFUN)+ggtitle("Drake 2 2004")


#d3.04
dr304<-gl.ibd(gl=NULL, Dgen=(MAF006.D3.04.dist), Dgeo=as.dist(coords.d3.04.snp.dist))
dr304Dgeo<-as.matrix(dr304$Dgeo)
dr304Dgen<-as.matrix(dr304$Dgen)

dr304Dgeo[upper.tri(dr304Dgeo, diag=T)]<-NA
dr304Dgen[upper.tri(dr304Dgen, diag=T)]<-NA


dr304df<-cbind(melt(dr304Dgeo),melt(dr304Dgen))
dr304df<-na.omit(dr304df)
colnames(dr304df)[c(3,6)]<-c("geo","gen")
summary(lm(dr304df$gen~dr304df$geo))
cor.test(y=dr304df$gen, x=dr304df$geo, method="p", alternative = "g")

pldr304<-ggplot(dr304df[,c(3,6)], aes(x=dr304df$geo, y=dr304df$gen))+
  #stat_density2d(aes(alpha=..level.., fill=..level..), size=2, 
  #               bins=100, geom="polygon") + 
  scale_fill_gradient(low = "yellow", high = "red") +
  scale_alpha(range = c(0, 1), guide = FALSE) +
  #geom_density2d(colour="white", bins=20) +
  geom_point(data = dr304df[,c(3,6)]) +
  stat_smooth(method="lm", level=.95, size=1, color="black")+
  #lims(x=c(min(dr304df$geo),max(dr304df$geo)), y=c(min(dr304df$gen),max(dr304df$gen)))+
  theme(axis.text.y   = element_text(size=10, angle=90, hjust=.5),
        axis.text.x   = element_text(size=10),
        axis.title.y  = element_text(size=12),
        axis.title.x  = element_text(size=12),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = "none",
  ) + labs(x="Distance (m)", y="Genetic Distance")+coord_cartesian(xlim=c(min(dr304df$geo),max(dr304df$geo)), ylim=c(min(dr304df$gen),max(dr304df$gen)))+
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE))+scale_x_continuous(labels=scaleFUN)+ggtitle("Drake 3 2004")



#d2.14
E<-gl.ibd(gl=NULL, Dgen=as.dist(d2.14.gen), Dgeo=as.dist(d2.14.geo))
EDgeo<-as.matrix(E$Dgeo)
EDgen<-as.matrix(E$Dgen)

EDgeo[upper.tri(EDgeo, diag=T)]<-NA
EDgen[upper.tri(EDgen, diag=T)]<-NA


Edf<-cbind(melt(EDgeo),melt(EDgen))
Edf<-na.omit(Edf)
colnames(Edf)[c(3,6)]<-c("geo","gen")
summary(lm(Edf$gen~Edf$geo))
cor.test(y= Edf$gen, x=Edf$geo, method="pearson", alternative = "g")

pldr214<-ggplot(Edf[,c(3,6)], aes(x=Edf$geo, y=Edf$gen))+
  #stat_density2d(aes(alpha=..level.., fill=..level..), size=2, 
  #              bins=100, geom="polygon") + 
  scale_fill_gradient(low = "yellow", high = "red") +
  scale_alpha(range = c(0, 1), guide = FALSE) +
  #geom_density2d(colour="white", bins=20) +
  geom_point(data = Edf[,c(3,6)]) +
  stat_smooth(method="lm", level=.95, size=1, color="black")+
  #lims(x=c(min(Edf$geo),max(Edf$geo)), y=c(min(Edf$gen),max(Edf$gen)))+
  theme(axis.text.y   = element_text(size=10, angle=90, hjust=0.5),
        axis.text.x   = element_text(size=10),
        axis.title.y  = element_text(size=12),
        axis.title.x  = element_text(size=12),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = "none",
  ) + labs(x="Distance (m)", y="Genetic Distance")+coord_cartesian(xlim=c(min(Edf$geo),max(Edf$geo)), ylim=c(min(Edf$gen),max(Edf$gen)))+
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE))+scale_x_continuous(labels=scaleFUN)+ggtitle("Drake 2 2014")


#d3.14
dr314<-gl.ibd(gl=NULL, Dgen=as.dist(d3.14.gen), Dgeo=as.dist(d3.14.geo))
dr314Dgeo<-as.matrix(dr314$Dgeo)
dr314Dgen<-as.matrix(dr314$Dgen)

dr314Dgeo[upper.tri(dr314Dgeo, diag=T)]<-NA
dr314Dgen[upper.tri(dr314Dgen, diag=T)]<-NA


dr314df<-cbind(melt(dr314Dgeo),melt(dr314Dgen))
dr314df<-na.omit(dr314df)
colnames(dr314df)[c(3,6)]<-c("geo","gen")
summary(lm(dr314df$gen~dr314df$geo))
cor.test(y=dr314df$gen,x=dr314df$geo, method="p", alternative = "g")

pldr314<-ggplot(dr314df[,c(3,6)], aes(x=dr314df$geo, y=dr314df$gen))+
  #stat_density2d(aes(alpha=..level.., fill=..level..), size=2, 
  #              bins=100, geom="polygon") + 
  scale_fill_gradient(low = "yellow", high = "red") +
  scale_alpha(range = c(0, 1), guide = FALSE) +
  #geom_density2d(colour="white", bins=20) +
  geom_point(data = dr314df[,c(3,6)]) +
  stat_smooth(method="lm", level=.95, size=1, color="black")+
  #lims(x=c(min(dr314df$geo),max(dr314df$geo)), y=c(min(dr314df$gen),max(dr314df$gen)))+
  theme(axis.text.y   = element_text(size=10, angle=90, hjust=0.5),
        axis.text.x   = element_text(size=10),
        axis.title.y  = element_text(size=12),
        axis.title.x  = element_text(size=12),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = "none",
  ) + labs(x="Distance (m)", y="Genetic Distance")+coord_cartesian(xlim=c(min(dr314df$geo),max(dr314df$geo)), ylim=c(min(dr314df$gen),max(dr314df$gen)))+
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE))+scale_x_continuous(labels=scaleFUN)+ggtitle("Drake 3 2014")


#d2.15
dr215<-gl.ibd(gl=NULL, Dgen=(MAF006.D2.15.dist), Dgeo=as.dist(coords.d2.15.snp.dist))
dr215Dgeo<-as.matrix(dr215$Dgeo)
dr215Dgen<-as.matrix(dr215$Dgen)

dr215Dgeo[upper.tri(dr215Dgeo, diag=T)]<-NA
dr215Dgen[upper.tri(dr215Dgen, diag=T)]<-NA


dr215df<-cbind(melt(dr215Dgeo),melt(dr215Dgen))
dr215df<-na.omit(dr215df)
colnames(dr215df)[c(3,6)]<-c("geo","gen")
summary(lm(dr215df$gen~dr215df$geo))
cor.test(y = dr215df$gen, x= dr215df$geo, method = "pearson", alternative = "g")

pldr215<-ggplot(dr215df[,c(3,6)], aes(x=dr215df$geo, y=dr215df$gen))+
  #stat_density2d(aes(alpha=..level.., fill=..level..), size=2, 
  #               bins=100, geom="polygon") + 
  scale_fill_gradient(low = "yellow", high = "red") +
  scale_alpha(range = c(0, 1), guide = FALSE) +
  #geom_density2d(colour="white", bins=20) +
  geom_point(data = dr215df[,c(3,6)]) +
  stat_smooth(method="lm", level=.95, size=1, color="black")+
  #lims(x=c(min(dr215df$geo),max(dr215df$geo)), y=c(min(dr215df$gen),max(dr215df$gen)))+
  theme(axis.text.y   = element_text(size=10, angle = 90, hjust = 0.5),
        axis.text.x   = element_text(size=10),
        axis.title.y  = element_text(size=12),
        axis.title.x  = element_text(size=12),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = "none",
  ) + labs(x="Distance (m)", y="Genetic Distance")+coord_cartesian(xlim=c(min(dr215df$geo),max(dr215df$geo)), ylim=c(min(dr215df$gen),max(dr215df$gen)))+
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE))+scale_x_continuous(labels=scaleFUN)+ggtitle("Drake 2 2015")

#d3.15
dr315<-gl.ibd(gl=NULL, Dgen=(MAF006.D3.15.dist), Dgeo=as.dist(coords.d3.15.snp.dist))
dr315Dgeo<-as.matrix(dr315$Dgeo)
dr315Dgen<-as.matrix(dr315$Dgen)

dr315Dgeo[upper.tri(dr315Dgeo, diag=T)]<-NA
dr315Dgen[upper.tri(dr315Dgen, diag=T)]<-NA


dr315df<-cbind(melt(dr315Dgeo),melt(dr315Dgen))
dr315df<-na.omit(dr315df)
colnames(dr315df)[c(3,6)]<-c("geo","gen")
summary(lm(dr315df$gen~dr315df$geo))
cor.test(y=dr315df$gen, x=dr315df$geo, method="p", alternative = "g")

pldr315<-ggplot(dr315df[,c(3,6)], aes(x=dr315df$geo, y=dr315df$gen))+
  #stat_density2d(aes(alpha=..level.., fill=..level..), size=2, 
  #               bins=100, geom="polygon") + 
  scale_fill_gradient(low = "yellow", high = "red") +
  scale_alpha(range = c(0, 1), guide = FALSE) +
  #geom_density2d(colour="white", bins=20) +
  geom_point(data = dr315df[,c(3,6)]) +
  stat_smooth(method="lm", level=.95, size=1, color="black")+
  #lims(x=c(min(dr315df$geo),max(dr315df$geo)), y=c(min(dr315df$gen),max(dr315df$gen)))+
  theme(axis.text.y   = element_text(size=10, angle=90, hjust=0.5),
        axis.text.x   = element_text(size=10),
        axis.title.y  = element_text(size=12),
        axis.title.x  = element_text(size=12),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = "none",
  ) + labs(x="Distance (m)", y="Genetic Distance")+coord_cartesian(xlim=c(min(dr315df$geo),max(dr315df$geo)), ylim=c(min(dr315df$gen),max(dr315df$gen)))+
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE))+scale_x_continuous(labels=scaleFUN)+ggtitle("Drake 3 2015")


grid.arrange(plmira, r, pldr204, pldr304, pldr214, pldr314, pldr215, pldr315, ncol=2, nrow=4)



############## Main Figure Autocorr

grid.arrange(plS, pldr204, plCE02, pldr214, plJL, pldr215, ncol=2, nrow=3)


#################################################################
#################################################################

dens <- MASS::kde2d(S$Dgeo, S$Dgen, n = 300)
myPal <- colorRampPalette(c("white", "blue", "gold", 
                            "orange", "red"))
plot(S$Dgeo, S$Dgen, pch = 20, cex = 0.8)
image(dens, col = transp(myPal(300), 0.7), add = TRUE)
points(S$Dgeo, S$Dgen, pch = 20, cex = 0.8)
abline(lm(S$Dgen ~ S$Dgeo))
title("Isolation by distance")

A<-gl.ibd(gl=NULL, Dgen=(hier.CESAC06.dist), Dgeo=as.dist(coords.CESAC06.aflp.dist))
B<-gl.ibd(gl=NULL, Dgen=(hier.JL.dist), Dgeo=as.dist(coords.JL.aflp.dist))
E<-gl.ibd(gl=NULL, Dgen=as.dist(d2.14.gen), Dgeo=as.dist(d2.14.geo))

###########################################
# Combine analyses

lat.dr1<-38.0545
long.dr1<--122.83343

lat.dr2<-38.054785
long.dr2<--122.833232

lat.dr3<-38.055212
long.dr3<--122.834134

lat.dr4<-38.054675
long.dr4<--122.836545


library(geosphere)
distm(c(long.dr1, lat.dr1), c(long.dr3, lat.dr3), fun = distHaversine)
distm(c(long.dr2, lat.dr2), c(long.dr3, lat.dr3), fun = distHaversine)
distm(c(long.dr4, lat.dr4), c(long.dr3, lat.dr3), fun = distHaversine)


lat.roch1<-43.233019
long.roch1<--77.554686

lat.roch2<-43.232747
long.roch2<--77.554917

lat.roch3<-43.232896
long.roch3<--77.554851

library(geosphere)

distm(c(long.roch1, lat.roch1), c(long.roch3, lat.roch3), fun = distHaversine)
distm(c(long.roch2, lat.roch2), c(long.roch3, lat.roch3), fun = distHaversine)

hier.gc.roch<-popsub(hier.gc, sublist=c("Rochester_NY_2007_1", "Rochester_NY_2007_2", "Rochester_NY_2007_3"))
plot(nj(dist(hier.gc.roch)))

####################################
####################################
####################################
#         Mantel combo


setwd("~/Dropbox/Genet Size/Jacob Analyses")
coords.trans<-read.csv("COMBINED_coords_Dr.HD.Roch.csv", head=T, stringsAsFactors = F)

coords.trans.dr<-coords.trans[1:97,]
ggplot(coords.trans.dr, aes(x=x.trans, y=y.trans, color=pop))+geom_point()

coords.trans.roch<-coords.trans[117:146,]
ggplot(coords.trans.roch, aes(x=x.trans, y=y.trans, color=pop))+geom_point()


coords.trans.HD<-coords.trans[98:116,]
ggplot(coords.trans.HD, aes(x=x.trans, y=y.trans, color=pop))+geom_point()
##############################################
# HD Mantel Combo
hier.gi.HD12<-popsub(hier.gi, sublist=c("Hearts_Desire_1_2004", "Hearts_Desire_2_2004"))
hier.gi.HD23<-popsub(hier.gi, sublist=c("Hearts_Desire_2_2004", "Hearts_Desire_3_2004"))
hier.gi.HD13<-popsub(hier.gi, sublist=c("Hearts_Desire_1_2004", "Hearts_Desire_3_2004"))

hier.gi.HD<-popsub(hier.gi, sublist=c("Hearts_Desire_1_2004","Hearts_Desire_2_2004", "Hearts_Desire_3_2004"))

hier.gi.HD12.dist<-as.matrix(diss.dist(hier.gi.HD12, percent = F))
hier.gi.HD23.dist<-as.matrix(diss.dist(hier.gi.HD23, percent = F))
hier.gi.HD13.dist<-as.matrix(diss.dist(hier.gi.HD13, percent = F))

hier.gi.HD23.dist <- hier.gi.HD23.dist[, order(as.double(colnames(hier.gi.HD23.dist)))]
hier.gi.HD23.dist <- hier.gi.HD23.dist[order(as.double(rownames(hier.gi.HD23.dist))),]

hier.gi.HD.dist<-as.matrix(diss.dist(hier.gi.HD, percent = F))
hier.gi.HD.dist <- hier.gi.HD.dist[, order(as.double(colnames(hier.gi.HD.dist)))]
hier.gi.HD.dist <- hier.gi.HD.dist[order(as.double(rownames(hier.gi.HD.dist))),]


coords.HD1<-coords.trans[98:104,]
coords.HD2<-coords.trans[105:112,]
coords.HD3<-coords.trans[113:116,]
coords.HD<-coords.trans[98:116,]

cords.HD12<-rbind(coords.HD1, coords.HD2)
cords.HD23<-rbind(coords.HD2, coords.HD3)
cords.HD13<-rbind(coords.HD1, coords.HD3)

coords.HD

cords.HD12.dist<-as.matrix(dist(cords.HD12[,c("x.trans", "y.trans")]))
colnames(cords.HD12.dist)<-rownames(cords.HD12.dist)<-cords.HD12$FINAlab
cords.HD12.dist<-cords.HD12.dist[-11,-11]

cords.HD23.dist<-as.matrix(dist(cords.HD23[,c("x.trans", "y.trans")]))
colnames(cords.HD23.dist)<-rownames(cords.HD23.dist)<-cords.HD23$FINAlab
cords.HD23.dist<-cords.HD23.dist[-4,-4]
cords.HD23.dist <- cords.HD23.dist[, order(as.integer(colnames(cords.HD23.dist)))]
cords.HD23.dist <- cords.HD23.dist[order(as.integer(colnames(cords.HD23.dist))),]

cords.HD13.dist<-as.matrix(dist(cords.HD13[,c("x.trans", "y.trans")]))
colnames(cords.HD13.dist)<-rownames(cords.HD13.dist)<-cords.HD13$FINAlab
cords.HD13.dist<-cords.HD13.dist


cords.HD.dist<-as.matrix(dist(coords.HD[,c("x.trans", "y.trans")]))
colnames(cords.HD.dist)<-rownames(cords.HD.dist)<-coords.HD$FINAlab
cords.HD.dist <- cords.HD.dist[, order(as.integer(colnames(cords.HD.dist)))]
cords.HD.dist <- cords.HD.dist[order(as.integer(colnames(cords.HD.dist))),]
cords.HD.dist<-cords.HD.dist[-11,-11]



mantel.test(cords.HD12.dist, hier.gi.HD12.dist, alternative = "greater")
gl.ibd(gl=NULL, Dgen=as.dist(hier.gi.HD12.dist), Dgeo=as.dist(cords.HD12.dist))

mantel.test(cords.HD23.dist, hier.gi.HD23.dist, alternative = "greater")
gl.ibd(gl=NULL, Dgen=as.dist(hier.gi.HD23.dist), Dgeo=as.dist(cords.HD23.dist))

mantel.test(cords.HD13.dist, hier.gi.HD13.dist, alternative = "greater")
gl.ibd(gl=NULL, Dgen=as.dist(hier.gi.HD13.dist), Dgeo=as.dist(cords.HD13.dist))


mantel.test(cords.HD.dist, hier.gi.HD.dist, alternative = "greater")
gl.ibd(gl=NULL, Dgen=as.dist(hier.gi.HD.dist), Dgeo=as.dist(cords.HD.dist))
mantel(cords.HD.dist, hier.gi.HD.dist, method = "p")




##############################################
# Rochester Mantel Combo
hier.gc.roch13<-popsub(hier.gc, sublist=c("Rochester_NY_2007_1", "Rochester_NY_2007_3"))
hier.gc.roch12<-popsub(hier.gc, sublist=c("Rochester_NY_2007_1", "Rochester_NY_2007_2"))
hier.gc.roch23<-popsub(hier.gc, sublist=c("Rochester_NY_2007_2", "Rochester_NY_2007_3"))

hier.gc.roch<-popsub(hier.gc, sublist=c("Rochester_NY_2007_1", "Rochester_NY_2007_2", "Rochester_NY_2007_3"))


hier.gc.roch13.dist<-as.matrix(diss.dist(hier.gc.roch13, percent = F))
hier.gc.roch12.dist<-as.matrix(diss.dist(hier.gc.roch12, percent = F))
hier.gc.roch12.dist<-hier.gc.roch12.dist[-19,-19]
hier.gc.roch23.dist<-as.matrix(diss.dist(hier.gc.roch23, percent = F))
hier.gc.roch23.dist<-hier.gc.roch23.dist[-6,-6]

hier.gc.roch.dist<-as.matrix(diss.dist(hier.gc.roch, percent = F))
hier.gc.roch.dist<-hier.gc.roch.dist[-19,-19]

coords.trans.roch
coords.roch1<-coords.trans.roch[1:13,]
coords.roch1.aflp<-coords.roch1

coords.roch2<-coords.trans.roch[14:20,]
coords.roch2.aflp<-coords.roch2[-4,] #need to remove 2_6b from genind file

coords.roch3<-coords.trans.roch[21:30,]
coords.roch3.aflp<-coords.roch3[3:4,]

coords.roch12<-rbind(coords.roch1.aflp, coords.roch2.aflp)
coords.roch23<-rbind(coords.roch2.aflp, coords.roch3.aflp)
coords.roch13<-rbind(coords.roch1.aflp, coords.roch3.aflp)

coords.roch<-rbind(coords.roch1.aflp, coords.roch2.aflp, coords.roch3.aflp)

coords.roch12.aflp.dist<-as.matrix(dist(coords.roch12[,c("x.trans", "y.trans")]))
colnames(coords.roch12.aflp.dist)<-rownames(coords.roch12.aflp.dist)<-coords.roch12$FINAlab

coords.roch23.aflp.dist<-as.matrix(dist(coords.roch23[,c("x.trans", "y.trans")]))
colnames(coords.roch23.aflp.dist)<-rownames(coords.roch23.aflp.dist)<-coords.roch23$FINAlab

coords.roch13.aflp.dist<-as.matrix(dist(coords.roch13[,c("x.trans", "y.trans")]))
colnames(coords.roch13.aflp.dist)<-rownames(coords.roch13.aflp.dist)<-coords.roch13$FINAlab

coords.roch.aflp.dist<-as.matrix(dist(coords.roch[,c("x.trans", "y.trans")]))
colnames(coords.roch.aflp.dist)<-rownames(coords.roch.aflp.dist)<-coords.roch$FINAlab

library(phangorn)
plot(upgma(dist(hier.gc.roch)))

mantel.test(coords.roch12.aflp.dist, hier.gc.roch12.dist, alternative = "greater")
gl.ibd(gl=NULL, Dgen=as.dist(hier.gc.roch12.dist), Dgeo=as.dist(coords.roch12.aflp.dist))

mantel.test(coords.roch23.aflp.dist, hier.gc.roch23.dist, alternative = "greater")
gl.ibd(gl=NULL, Dgen=as.dist(hier.gc.roch23.dist), Dgeo=as.dist(coords.roch23.aflp.dist))

mantel.test(coords.roch13.aflp.dist, hier.gc.roch13.dist, alternative = "greater")
gl.ibd(gl=NULL, Dgen=as.dist(hier.gc.roch13.dist), Dgeo=as.dist(coords.roch13.aflp.dist))

mantel.test(coords.roch.aflp.dist, hier.gc.roch.dist, alternative = "greater")
gl.ibd(gl=NULL, Dgen=as.dist(hier.gc.roch.dist), Dgeo=as.dist(coords.roch.aflp.dist))
mantel(coords.roch.aflp.dist, hier.gc.roch.dist, method="p")

##############################################
#Drake Mantel Combo
#AFLP

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

hier.gc<-as.genclone(hier.gi)


hier.gi.dr12<-popsub(hier.gi, sublist=c("Drake_1_2004", "Drake_2_2004"))
hier.gi.dr13<-popsub(hier.gi, sublist=c("Drake_1_2004", "Drake_3_2004"))
hier.gi.dr14<-popsub(hier.gi, sublist=c("Drake_1_2004", "Drake_4_2004"))
hier.gi.dr23<-popsub(hier.gi, sublist=c("Drake_2_2004", "Drake_3_2004"))
hier.gi.dr24<-popsub(hier.gi, sublist=c("Drake_2_2004", "Drake_4_2004"))
hier.gi.dr34<-popsub(hier.gi, sublist=c("Drake_3_2004", "Drake_4_2004"))

hier.gi.dr<-popsub(hier.gi, sublist=c("Drake_1_2004", "Drake_2_2004", "Drake_3_2004", "Drake_4_2004"))

hier.gi.dr12.dist<-as.matrix(diss.dist(hier.gi.dr12, percent = F))
hier.gi.dr13.dist<-as.matrix(diss.dist(hier.gi.dr13, percent = F))
hier.gi.dr14.dist<-as.matrix(diss.dist(hier.gi.dr14, percent = F))
hier.gi.dr23.dist<-as.matrix(diss.dist(hier.gi.dr23, percent = F))
hier.gi.dr24.dist<-as.matrix(diss.dist(hier.gi.dr24, percent = F))
hier.gi.dr34.dist<-as.matrix(diss.dist(hier.gi.dr34, percent = F))

hier.gi.dr.dist<-as.matrix(bitwise.dist(hier.gi.dr, percent = F, missing_match = T, scale_missing = T))

coords.trans.dr1.aflp<-coords.trans[1:19,]
coords.trans.dr2.aflp<-coords.trans[20:33,]
coords.trans.dr2.aflp<-coords.trans.dr2.aflp[-9,]
coords.trans.dr3.aflp<-coords.trans[74:81,]
coords.trans.dr4.aflp<-coords.trans[91:97,]
coords.trans.dr4.aflp<-coords.trans.dr4.aflp[-2,]

coords.trans.dr12<-rbind(coords.trans.dr1.aflp, coords.trans.dr2.aflp)
coords.trans.dr13<-rbind(coords.trans.dr1.aflp, coords.trans.dr3.aflp)
coords.trans.dr14<-rbind(coords.trans.dr1.aflp, coords.trans.dr4.aflp)
coords.trans.dr23<-rbind(coords.trans.dr2.aflp, coords.trans.dr3.aflp)
coords.trans.dr24<-rbind(coords.trans.dr2.aflp, coords.trans.dr4.aflp)
coords.trans.dr34<-rbind(coords.trans.dr3.aflp, coords.trans.dr4.aflp)

coords.trans.dr<-rbind(coords.trans.dr1.aflp, coords.trans.dr2.aflp, coords.trans.dr3.aflp, coords.trans.dr4.aflp)


coords.trans.dr12.dist<-as.matrix(dist(coords.trans.dr12[,c("x.trans", "y.trans")]))
colnames(coords.trans.dr12.dist)<-rownames(coords.trans.dr12.dist)<-coords.trans.dr12$FINAlab

coords.trans.dr13.dist<-as.matrix(dist(coords.trans.dr13[,c("x.trans", "y.trans")]))
colnames(coords.trans.dr13.dist)<-rownames(coords.trans.dr13.dist)<-coords.trans.dr13$FINAlab

coords.trans.dr14.dist<-as.matrix(dist(coords.trans.dr14[,c("x.trans", "y.trans")]))
colnames(coords.trans.dr14.dist)<-rownames(coords.trans.dr14.dist)<-coords.trans.dr14$FINAlab

coords.trans.dr23.dist<-as.matrix(dist(coords.trans.dr23[,c("x.trans", "y.trans")]))
colnames(coords.trans.dr23.dist)<-rownames(coords.trans.dr23.dist)<-coords.trans.dr23$FINAlab

coords.trans.dr24.dist<-as.matrix(dist(coords.trans.dr24[,c("x.trans", "y.trans")]))
colnames(coords.trans.dr24.dist)<-rownames(coords.trans.dr24.dist)<-coords.trans.dr24$FINAlab

coords.trans.dr34.dist<-as.matrix(dist(coords.trans.dr34[,c("x.trans", "y.trans")]))
colnames(coords.trans.dr34.dist)<-rownames(coords.trans.dr34.dist)<-coords.trans.dr34$FINAlab

coords.trans.dr.dist<-as.matrix(dist(coords.trans.dr[,c("x.trans", "y.trans")]))
colnames(coords.trans.dr.dist)<-rownames(coords.trans.dr.dist)<-coords.trans.dr$FINAlab




mantel.test(coords.trans.dr12.dist, hier.gi.dr12.dist, alternative="greater")
gl.ibd(gl=NULL, Dgen=as.dist(hier.gi.dr12.dist), Dgeo=as.dist(coords.trans.dr12.dist))

mantel.test(coords.trans.dr13.dist, hier.gi.dr13.dist, alternative="greater")
gl.ibd(gl=NULL, Dgen=as.dist(hier.gi.dr13.dist), Dgeo=as.dist(coords.trans.dr13.dist))

mantel.test(coords.trans.dr14.dist, hier.gi.dr14.dist, alternative="greater")
gl.ibd(gl=NULL, Dgen=as.dist(hier.gi.dr14.dist), Dgeo=as.dist(coords.trans.dr14.dist))

mantel.test(coords.trans.dr23.dist, hier.gi.dr23.dist, alternative="greater")
gl.ibd(gl=NULL, Dgen=as.dist(hier.gi.dr23.dist), Dgeo=as.dist(coords.trans.dr23.dist))

mantel.test(coords.trans.dr24.dist, hier.gi.dr24.dist, alternative="greater")
gl.ibd(gl=NULL, Dgen=as.dist(hier.gi.dr24.dist), Dgeo=as.dist(coords.trans.dr24.dist))

mantel.test(coords.trans.dr34.dist, hier.gi.dr34.dist, alternative="greater")
gl.ibd(gl=NULL, Dgen=as.dist(hier.gi.dr34.dist), Dgeo=as.dist(coords.trans.dr34.dist))

mantel.test(coords.trans.dr.dist, hier.gi.dr.dist, alternative="greater")
gl.ibd(gl=NULL, Dgen=as.dist(hier.gi.dr.dist), Dgeo=as.dist(coords.trans.dr.dist))
mantel(coords.trans.dr.dist, hier.gi.dr.dist, method="p")

#SNP
coords.trans

coords.trans.dr2.04<-coords.trans[20:34,]
coords.trans.dr2.04<-coords.trans.dr2.04[-c(1,15),]
coords.trans.dr3.04<-coords.trans[74:81,]
coords.trans.dr3.04<-coords.trans.dr3.04[-c(1,2,8),]

coords.trans.dr2.14<-coords.trans[46:70,]
coords.trans.dr2.14<-coords.trans.dr2.14[-8,]
coords.trans.dr2.14 <- coords.trans.dr2.14[order(coords.trans.dr2.14$FINAlab),] 
coords.trans.dr3.14<-coords.trans[82:90,]
coords.trans.dr3.14<-coords.trans.dr3.14[order(coords.trans.dr3.14$FINAlab),]

coords.trans.dr2.15<-coords.trans[34:45,]
coords.trans.dr2.15<-coords.trans.dr2.15[-8,]
coords.trans.dr3.15<-coords.trans[71:73,]

coords.trans.dr23.04.snp<-rbind(coords.trans.dr2.04, coords.trans.dr3.04)
coords.trans.dr23.14.snp<-rbind(coords.trans.dr2.14, coords.trans.dr3.14)
coords.trans.dr23.15.snp<-rbind(coords.trans.dr2.15, coords.trans.dr3.15)

coords.trans.dr23.04.snp.dist<-as.matrix(dist(coords.trans.dr23.04.snp[,c("x.trans", "y.trans")]))
colnames(coords.trans.dr23.04.snp.dist)<-rownames(coords.trans.dr23.04.snp.dist)<-coords.trans.dr23.04.snp$FINAlab

coords.trans.dr23.14.snp.dist<-as.matrix(dist(coords.trans.dr23.14.snp[,c("x.trans", "y.trans")]))
colnames(coords.trans.dr23.14.snp.dist)<-rownames(coords.trans.dr23.14.snp.dist)<-coords.trans.dr23.14.snp$FINAlab

coords.trans.dr23.15.snp.dist<-as.matrix(dist(coords.trans.dr23.15.snp[,c("x.trans", "y.trans")]))
colnames(coords.trans.dr23.15.snp.dist)<-rownames(coords.trans.dr23.15.snp.dist)<-coords.trans.dr23.15.snp$FINAlab


setwd("~/Dropbox/finalizing_genets")
maff006<-read.PLINK("HESS_fix/mx0.5/vcf.MMDP60.maf0.006.minQ30.recode.raw", 
                      map.file = "HESS_fix/mx0.5/vcf.MMDP60.maf0.006.minQ30.recode.map", quiet = FALSE,
                      parallel = require("parallel"), n.cores = NULL)

myPops <- as.factor(c('old_EU',	'old_EU',	'old_EU',	'old_EU',	'old_EU',	'old_EU',	'CA_old',	'old_EU',	'D2.04',	'D2.04',	
                      'D2.04',	'D2.04',	'D2.04',	'D2.04',	'D2.04',	'D2.04',	'D2.04',	'D2.04',	'D2.04',	'D2.04',
                      'D2.04',	'D3.04',	'D3.04',	'D3.04',	'D3.04',	'D3.04',	'old_EU',	'D2.14',	'D2.14',	'D2.14',
                      'D2.14',	'D2.14',	'D2.14',	'D2.14',	'D2.14',	'D2.14',	'D2.14',	'D2.14',	'D2.14',	'D2.14',
                      'D2.14',	'D2.14',	'D2.14',	'D2.14',	'D2.14',	'D2.14',	'D2.14',	'D2.14',	'D2.14',	'D2.14',
                      'D2.14',	'D2.14',	'D3.14',	'D3.14',	'D3.14',	'D3.14',	'D3.14',	'D3.14',	'D3.14',	'D3.14',
                      'D3.14',	'Agraria',	'Agraria',	'Vilarinho',	'Vilarinho',	'Vilarinho',	'Vilarinho',	'Vilarinho',
                      'Mira',	'Mira',	'Mira',	'Mira',	'D2.15',	'D2.15',	'D2.15',	'D2.15',	'D2.15',	'D2.15',	'D2.15',
                      'D2.15',	'D2.15',	'D2.15',	'D2.15',	'D3.15',	'D3.15',	'D3.15'))
pop(maff006)<-myPops

maff006.dr04<-popsub(maff006, sublist = c("D2.04", "D3.04"))
maff006.dr14<-popsub(maff006, sublist = c("D2.14", "D3.14"))
maff006.dr15<-popsub(maff006, sublist = c("D2.15", "D3.15"))


maff006.dr04.dist<-as.matrix(bitwise.dist(maff006.dr04, percent = F, missing_match = T, scale_missing = T))
maff006.dr14.dist<-as.matrix(bitwise.dist(maff006.dr14, percent = F, missing_match = T, scale_missing = T))
maff006.dr14.dist<-maff006.dr14.dist[-10,-10]
maff006.dr15.dist<-as.matrix(bitwise.dist(maff006.dr15, percent = F, missing_match = T, scale_missing = T))


mantel.test(coords.trans.dr23.04.snp.dist, maff006.dr04.dist, alternative = "greater")
gl.ibd(gl=NULL, Dgen=as.dist(maff006.dr04.dist), Dgeo=as.dist(coords.trans.dr23.04.snp.dist))
mantel(coords.trans.dr23.04.snp.dist, maff006.dr04.dist, method="p")

mantel.test(coords.trans.dr23.14.snp.dist, maff006.dr14.dist, alternative = "greater")
gl.ibd(gl=NULL, Dgen=as.dist(maff006.dr14.dist), Dgeo=as.dist(coords.trans.dr23.14.snp.dist))
mantel(coords.trans.dr23.14.snp.dist, maff006.dr14.dist, method="p")
hist(maff006.dr14.dist)
hist(coords.trans.dr23.14.snp.dist)

mantel.test(coords.trans.dr23.15.snp.dist, maff006.dr15.dist, alternative = "greater")
gl.ibd(gl=NULL, Dgen=as.dist(maff006.dr15.dist), Dgeo=as.dist(coords.trans.dr23.15.snp.dist))
mantel(coords.trans.dr23.15.snp.dist, maff006.dr15.dist, method="p")

####################################################
#Mantel between PCoA clusters
setwd("~/Dropbox/Genet Size/Jacob Analyses")


coords.relab<-read.csv("subpops.csv", head=T, stringsAsFactors = F)
# coords.relab$subpop<-as.character(coords.relab$Sample_Name)
# coords.relab$subpop<-revalue(coords.relab$subpop, c(
#   '10004'="EU",	'10007'="EU",	'10016'="EU",	'10018'="EU",	'10019'="EU",	'10169'="EU",	'10170'="TRH",	'10171'="EU",	'10221'="D2.04",	'10222'="D2.04",	'10223'="D2.04",	
#   '10224'="D2.04",	'10225'="D2.04",	'10226'="D2s.04",'10227'="D2.04",	'10228'="D2.04",	'10229'="D2.04",	'10230'="D2.04",	'10231'="D2.04",	'10232'="D2.04",	'10233'="D2.04",	'10237'="D3",	
#   '10238'="D3",	'10239'="D3",	'10240'="D3",	'10241'="D3",	'10277'="EU",	'10280'="D2s.14",'10281'="D2s.14",'10282'="D2s.14",'10283'="D2s.14",'10287'="D2.14",	'10288'="D2.14",	
#   '10292'="D2.14",	'10293'="D2.14",	'10294'="D2.14",	'10295'="D2.14",	'10298'="D2.14",	'10299'="D2.14",	'10300'="D2.14",	'10301'="D2.14",	'10303'="D2.14",	'10304'="D2s.14",'10306'="D2s.14",	
#   '10309'="D2.14",	'10326'="D2s.14",'10327'="D2s.14",'10328'="D2s.14",'10329'="D2s.14",'10330'="D2s.14",'10331'="D2s.14",'10334'="D2",'10347'="D3",'10348'="D3",'10349'="D3",	
#   '10350'="D3",	'10354'="D3",	'10355'="D3",	'10356'="D3",	'10380'="D3",	'10384'="D3",	'10502'="EU",	'10503'="EU",	'10504'="EU",	'10505'="EU",	'10506'="EU",	
#   '10508'="EU",	'10509'="EU",	'10510'="EU",	'10511'="EU",	'10512'="EU",	'10513'="EU",	'10707'="D2.15",	'10708'="D2.15",	'10709'="D2.15",	'10710'="D2s.15",'10711'="D2s.15",	
#   '10712'="D2s.15",'10713'="D2.15",	'10715'="D2s.15",'10716'="D2.15",	'10717'="D2.15",	'10718'="D2.15",	'10719'="D3",	'10720'="D3",	'10721'="D3"
# ))

maff006.dr04<-popsub(maff006, sublist = c("D2.04", "D3.04"))
maff006.dr14<-popsub(maff006, sublist = c("D2.14", "D3.14"))
maff006.dr15<-popsub(maff006, sublist = c("D2.15", "D3.15"))

setwd("~/Dropbox/Genet Size/Jacob Analyses")
coords.trans<-read.csv("COMBINED_coords_Dr.HD.Roch.csv", head=T, stringsAsFactors = F)



#d2 2014 cluster 1
maff006.dr14.cl1vd3 <- maff006.dr14[!(indNames(maff006.dr14) %in% coords.relab[coords.relab$subpop %in% c("D2s.14"), ]$Sample_Name)]
maff006.dr14.cl1vd3.dist<-as.matrix(bitwise.dist(maff006.dr14.cl1vd3, percent = F, missing_match = T, scale_missing = T))

coords.trans.cl1vd3<-coords.trans[coords.trans$FINAlab %in% indNames(maff006.dr14.cl1vd3),]
rownames(coords.trans.cl1vd3)<-coords.trans.cl1vd3$FINAlab
coords.trans.cl1vd3.dist<-as.matrix(dist(coords.trans.cl1vd3[,c("x.trans", "y.trans")]))

mantel.test(coords.trans.cl1vd3.dist, maff006.dr14.cl1vd3.dist, alternative = "greater")
gl.ibd(gl=NULL, Dgen=as.dist(maff006.dr14.cl1vd3.dist), Dgeo=as.dist(coords.trans.cl1vd3.dist))
mantel(coords.trans.cl1vd3.dist, maff006.dr14.cl1vd3.dist, method="p")

#d2 2014 cluster 2
maff006.dr14.cl2vd3 <- maff006.dr14[!(indNames(maff006.dr14) %in% coords.relab[coords.relab$subpop %in% c("D2.14"), ]$Sample_Name)]
maff006.dr14.cl2vd3.dist<-as.matrix(bitwise.dist(maff006.dr14.cl2vd3, percent = F, missing_match = T, scale_missing = T))

coords.trans.cl2vd3<-coords.trans[coords.trans$FINAlab %in% indNames(maff006.dr14.cl2vd3),]
rownames(coords.trans.cl2vd3)<-coords.trans.cl2vd3$FINAlab
coords.trans.cl2vd3.dist<-as.matrix(dist(coords.trans.cl2vd3[,c("x.trans", "y.trans")]))

mantel.test(coords.trans.cl2vd3.dist, maff006.dr14.cl2vd3.dist, alternative = "greater")
gl.ibd(gl=NULL, Dgen=as.dist(maff006.dr14.cl2vd3.dist), Dgeo=as.dist(coords.trans.cl2vd3.dist))
mantel(coords.trans.cl2vd3.dist, maff006.dr14.cl2vd3.dist, method="p")



#d1 2015 cluster 1
maff006.dr15.cl1vd3 <- maff006.dr15[!(indNames(maff006.dr15) %in% coords.relab[coords.relab$subpop %in% c("D2s.15"), ]$Sample_Name)]
maff006.dr15.cl1vd3.dist<-as.matrix(bitwise.dist(maff006.dr15.cl1vd3, percent = F, missing_match = T, scale_missing = T))

coords.trans.cl1vd3.15<-coords.trans[coords.trans$FINAlab %in% indNames(maff006.dr15.cl1vd3),]
rownames(coords.trans.cl1vd3.15)<-coords.trans.cl1vd3.15$FINAlab
coords.trans.cl1vd3.15.dist<-as.matrix(dist(coords.trans.cl1vd3.15[,c("x.trans", "y.trans")]))

mantel.test(coords.trans.cl1vd3.dist, maff006.dr14.cl1vd3.dist, alternative = "greater")
gl.ibd(gl=NULL, Dgen=as.dist(maff006.dr14.cl1vd3.dist), Dgeo=as.dist(coords.trans.cl1vd3.dist))
mantel(coords.trans.cl1vd3.dist, maff006.dr14.cl1vd3.dist, method="p")


#d1 2015 cluster 2
maff006.dr15.cl2vd3 <- maff006.dr15[!(indNames(maff006.dr15) %in% coords.relab[coords.relab$subpop %in% c("D2.15"), ]$Sample_Name)]
maff006.dr15.cl2vd3.dist<-as.matrix(bitwise.dist(maff006.dr15.cl2vd3, percent = F, missing_match = T, scale_missing = T))

coords.trans.cl2vd3.15<-coords.trans[coords.trans$FINAlab %in% indNames(maff006.dr15.cl2vd3),]
rownames(coords.trans.cl2vd3.15)<-coords.trans.cl2vd3.15$FINAlab
coords.trans.cl2vd3.15.dist<-as.matrix(dist(coords.trans.cl2vd3.15[,c("x.trans", "y.trans")]))

mantel.test(coords.trans.cl2vd3.dist, maff006.dr14.cl2vd3.dist, alternative = "greater")
gl.ibd(gl=NULL, Dgen=as.dist(maff006.dr14.cl2vd3.dist), Dgeo=as.dist(coords.trans.cl2vd3.dist))
mantel(coords.trans.cl2vd3.dist, maff006.dr14.cl2vd3.dist, method="p")

# coords.relab.drake<-coords.relab[coords.relab$pop %in% c("D2.04", "D2.14", "D2.15"), ]
# write.csv(coords.relab.drake[,c("Sample_Name", "subpop")], "subpopsD2.4denny.csv")



#within each cluster
#d2.14 cluster 1
maff006.dr14.cl1 <- maff006.dr14[(indNames(maff006.dr14) %in% coords.relab[coords.relab$subpop %in% c("D2.14"), ]$Sample_Name)]
maff006.dr14.cl1.dist<-as.matrix(bitwise.dist(maff006.dr14.cl1, percent = F, missing_match = T, scale_missing = T))

coords.trans.cl1<-coords.trans[coords.trans$FINAlab %in% indNames(maff006.dr14.cl1),]
rownames(coords.trans.cl1)<-coords.trans.cl1$FINAlab
coords.trans.cl1.dist<-as.matrix(dist(coords.trans.cl1[,c("x.trans", "y.trans")]))

mantel.test(coords.trans.cl1.dist, maff006.dr14.cl1.dist, alternative = "greater")
gl.ibd(gl=NULL, Dgen=as.dist(maff006.dr14.cl1.dist), Dgeo=as.dist(coords.trans.cl1.dist))
mantel(coords.trans.cl1.dist, maff006.dr14.cl1.dist, method="p")

#d2.14 cluster 2
maff006.dr14.cl2 <- maff006.dr14[(indNames(maff006.dr14) %in% coords.relab[coords.relab$subpop %in% c("D2s.14"), ]$Sample_Name)]
maff006.dr14.cl2.dist<-as.matrix(bitwise.dist(maff006.dr14.cl2, percent = F, missing_match = T, scale_missing = T))

coords.trans.cl2<-coords.trans[coords.trans$FINAlab %in% indNames(maff006.dr14.cl2),]
rownames(coords.trans.cl2)<-coords.trans.cl2$FINAlab
coords.trans.cl2.dist<-as.matrix(dist(coords.trans.cl2[,c("x.trans", "y.trans")]))

mantel.test(coords.trans.cl2.dist, maff006.dr14.cl2.dist, alternative = "greater")
gl.ibd(gl=NULL, Dgen=as.dist(maff006.dr14.cl2.dist), Dgeo=as.dist(coords.trans.cl2.dist))
mantel(coords.trans.cl2.dist, maff006.dr14.cl2.dist, method="p")


#d2.15 cluster 1
maff006.dr15.cl1 <- maff006.dr15[(indNames(maff006.dr15) %in% coords.relab[coords.relab$subpop %in% c("D2.15"), ]$Sample_Name)]
maff006.dr15.cl1.dist<-as.matrix(bitwise.dist(maff006.dr15.cl1, percent = F, missing_match = T, scale_missing = T))

coords.trans.cl1.15<-coords.trans[coords.trans$FINAlab %in% indNames(maff006.dr15.cl1),]
rownames(coords.trans.cl1.15)<-coords.trans.cl1.15$FINAlab
coords.trans.cl1.15.dist<-as.matrix(dist(coords.trans.cl1.15[,c("x.trans", "y.trans")]))

mantel.test(coords.trans.cl1.15.dist, maff006.dr15.cl1.dist, alternative = "greater")
gl.ibd(gl=NULL, Dgen=as.dist(maff006.dr15.cl1.dist), Dgeo=as.dist(coords.trans.cl1.15.dist))
mantel(coords.trans.cl1.15.dist, maff006.dr15.cl1.dist, method="p")

#d2.15 cluster 2
maff006.dr15.cl2 <- maff006.dr15[(indNames(maff006.dr15) %in% coords.relab[coords.relab$subpop %in% c("D2s.15"), ]$Sample_Name)]
maff006.dr15.cl2.dist<-as.matrix(bitwise.dist(maff006.dr15.cl2, percent = F, missing_match = T, scale_missing = T))

coords.trans.cl2.15<-coords.trans[coords.trans$FINAlab %in% indNames(maff006.dr15.cl2),]
rownames(coords.trans.cl2.15)<-coords.trans.cl2.15$FINAlab
coords.trans.cl2.15.dist<-as.matrix(dist(coords.trans.cl2.15[,c("x.trans", "y.trans")]))

mantel.test(coords.trans.cl2.15.dist, maff006.dr15.cl2.dist, alternative = "greater")
gl.ibd(gl=NULL, Dgen=as.dist(maff006.dr15.cl2.dist), Dgeo=as.dist(coords.trans.cl2.15.dist))
mantel(coords.trans.cl2.15.dist, maff006.dr15.cl2.dist, method="p")









#######################################################################################
#### time
maff006.dr2.0414<-popsub(maff006, sublist = c("D2.04", "D2.14"))
maff006.dr2.0415<-popsub(maff006, sublist = c("D2.04", "D2.15"))
maff006.dr2.1415<-popsub(maff006, sublist = c("D2.14", "D2.15"))

maff006.dr3.0414<-popsub(maff006, sublist = c("D3.04", "D3.14"))
maff006.dr3.0415<-popsub(maff006, sublist = c("D3.04", "D3.15"))
maff006.dr3.1415<-popsub(maff006, sublist = c("D3.14", "D3.15"))


maff006.dr2.0414.dist<-as.matrix(dist(maff006.dr2.0414))
maff006.dr2.0415.dist<-as.matrix(dist(maff006.dr2.0415))
maff006.dr2.1415.dist<-as.matrix(dist(maff006.dr2.1415))

maff006.dr3.0414.dist<-as.matrix(dist(maff006.dr3.0414))
maff006.dr3.0415.dist<-as.matrix(dist(maff006.dr3.0415))
maff006.dr3.1415.dist<-as.matrix(dist(maff006.dr3.1415))

coords.trans.dr2.0414.snp<-rbind(coords.trans.dr2.04, coords.trans.dr2.14)
coords.trans.dr2.0415.snp<-rbind(coords.trans.dr2.04, coords.trans.dr2.15)
coords.trans.dr2.1415.snp<-rbind(coords.trans.dr2.14, coords.trans.dr2.15)

coords.trans.dr3.0414.snp<-rbind(coords.trans.dr3.04, coords.trans.dr3.14)
coords.trans.dr3.0415.snp<-rbind(coords.trans.dr3.04, coords.trans.dr3.15)
coords.trans.dr3.1415.snp<-rbind(coords.trans.dr3.14, coords.trans.dr3.15)

coords.trans.dr2.0414.snp.dist<-as.matrix(dist(coords.trans.dr2.0414.snp[,5:6]))
colnames(coords.trans.dr2.0414.snp.dist)<-rownames(coords.trans.dr2.0414.snp.dist)<-coords.trans.dr2.0414.snp$FINAlab
coords.trans.dr2.0415.snp.dist<-as.matrix(dist(coords.trans.dr2.0415.snp[,5:6]))
colnames(coords.trans.dr2.0415.snp.dist)<-rownames(coords.trans.dr2.0415.snp.dist)<-coords.trans.dr2.0415.snp$FINAlab
coords.trans.dr2.1415.snp.dist<-as.matrix(dist(coords.trans.dr2.1415.snp[,5:6]))
colnames(coords.trans.dr2.1415.snp.dist)<-rownames(coords.trans.dr2.1415.snp.dist)<-coords.trans.dr2.1415.snp$FINAlab

coords.trans.dr3.0414.snp.dist<-as.matrix(dist(coords.trans.dr3.0414.snp[,5:6]))
colnames(coords.trans.dr3.0414.snp.dist)<-rownames(coords.trans.dr3.0414.snp.dist)<-coords.trans.dr3.0414.snp$FINAlab
coords.trans.dr3.0415.snp.dist<-as.matrix(dist(coords.trans.dr3.0415.snp[,5:6]))
colnames(coords.trans.dr3.0415.snp.dist)<-rownames(coords.trans.dr3.0415.snp.dist)<-coords.trans.dr3.0415.snp$FINAlab
coords.trans.dr3.1415.snp.dist<-as.matrix(dist(coords.trans.dr3.1415.snp[,5:6]))
colnames(coords.trans.dr3.1415.snp.dist)<-rownames(coords.trans.dr3.1415.snp.dist)<-coords.trans.dr3.1415.snp$FINAlab


maff006.dr2.0414.dist<-maff006.dr2.0414.dist[-23,-23]
mantel.test(coords.trans.dr2.0414.snp.dist, maff006.dr2.0414.dist, alternative = "greater")
mantel.test(coords.trans.dr2.0415.snp.dist, maff006.dr2.0415.dist, alternative = "greater")
maff006.dr2.1415.dist<-maff006.dr2.1415.dist[-10,-10]
mantel.test(coords.trans.dr2.1415.snp.dist, maff006.dr2.1415.dist, alternative = "greater")

mantel.test(coords.trans.dr3.0414.snp.dist, maff006.dr3.0414.dist, alternative = "greater")
mantel.test(coords.trans.dr3.0415.snp.dist, maff006.dr3.0415.dist, alternative = "greater")
mantel.test(coords.trans.dr3.1415.snp.dist, maff006.dr3.1415.dist, alternative = "greater")

################################################################################################
################################################################################################
################################################################################################
#######################     Pop binding box, follows from MANTEL loadings      #################
################################################################################################
################################################################################################
################################################################################################
################################################################################################



library(alphahull)  
library(pracma)
# Exposes ashape()
MBR <- function(points) {
  # Analyze the convex hull edges                       
  a <- ashape(points, alpha=1000)                 # One way to get a convex hull...
  e <- a$edges[, 5:6] - a$edges[, 3:4]            # Edge directions
  norms <- apply(e, 1, function(x) sqrt(x %*% x)) # Edge lengths
  v <- diag(1/norms) %*% e                        # Unit edge directions
  w <- cbind(-v[,2], v[,1])                       # Normal directions to the edges
  
  # Find the MBR
  vertices <- (points) [a$alpha.extremes, 1:2]    # Convex hull vertices
  minmax <- function(x) c(min(x), max(x))         # Computes min and max
  x <- apply(vertices %*% t(v), 2, minmax)        # Extremes along edges
  y <- apply(vertices %*% t(w), 2, minmax)        # Extremes normal to edges
  areas <- (y[1,]-y[2,])*(x[1,]-x[2,])            # Areas
  k <- which.min(areas)                           # Index of the best edge (smallest area)
  
  # Form a rectangle from the extremes of the best edge
  cbind(x[c(1,2,2,1,1),k], y[c(1,1,2,2,1),k]) %*% rbind(v[k,], w[k,])
}


# Create sample data
set.seed(23)
points <- matrix(rnorm(20*2), ncol=2)                 # Random (normally distributed) points
mbr <- MBR(points)

# Plot the hull, the MBR, and the points
limits <- apply(mbr, 2, function(x) c(min(x),max(x))) # Plotting limits
plot(ashape(points, alpha=1000), col="Gray", pch=20, 
     xlim=limits[,1], ylim=limits[,2])                # The hull
lines(mbr, col="Blue", lwd=3)                         # The MBR
points(points, pch=19)                                # The points


polyarea(mbr[,1], mbr[,2])

#########################################################



d1.hull<-as.matrix(coords.d1.aflp[,5:6])
d1.mbr<-MBR(d1.hull)
limits <- apply(d1.mbr, 2, function(x) c(min(x),max(x))) # Plotting limits
plot(ashape(d1.hull, alpha=1000), col="Gray", pch=20, 
     xlim=limits[,1], ylim=limits[,2])                # The hull
lines(d1.mbr, col="Blue", lwd=3)                         # The MBR
points(d1.hull, pch=19)    
polyarea(d1.mbr[,1], d1.mbr[,2])


coords.d2.04.aflp.hull<-as.matrix(coords.d2.04.aflp[,5:6])
coords.d2.04.aflp.mbr<-MBR(coords.d2.04.aflp.hull)
limits <- apply(coords.d2.04.aflp.mbr, 2, function(x) c(min(x),max(x))) # Plotting limits
plot(ashape(coords.d2.04.aflp.hull, alpha=1000), col="Gray", pch=20, 
     xlim=limits[,1], ylim=limits[,2])                # The hull
lines(coords.d2.04.aflp.mbr, col="Blue", lwd=3)                         # The MBR
points(coords.d2.04.aflp.hull, pch=19)    
polyarea(coords.d2.04.aflp.mbr[,1], coords.d2.04.aflp.mbr[,2])


coords.d3.04.aflp.hull<-as.matrix(coords.d3.04.aflp[,5:6])
coords.d3.04.aflp.mbr<-MBR(coords.d3.04.aflp.hull)
limits <- apply(coords.d3.04.aflp.mbr, 2, function(x) c(min(x),max(x))) # Plotting limits
plot(ashape(coords.d3.04.aflp.hull, alpha=1000), col="Gray", pch=20, 
     xlim=limits[,1], ylim=limits[,2])                # The hull
lines(coords.d3.04.aflp.mbr, col="Blue", lwd=3)                         # The MBR
points(coords.d3.04.aflp.hull, pch=19)    
polyarea(coords.d3.04.aflp.mbr[,1], coords.d3.04.aflp.mbr[,2])


coords.d4.aflp.hull<-as.matrix(coords.d4.aflp[,5:6])
coords.d4.aflp.mbr<-MBR(coords.d4.aflp.hull)
limits <- apply(coords.d4.aflp.mbr, 2, function(x) c(min(x),max(x))) # Plotting limits
plot(ashape(coords.d4.aflp.hull, alpha=1000), col="Gray", pch=20, 
     xlim=limits[,1], ylim=limits[,2])                # The hull
lines(coords.d4.aflp.mbr, col="Blue", lwd=3)                         # The MBR
points(coords.d4.aflp.hull, pch=19)    
polyarea(coords.d4.aflp.mbr[,1], coords.d4.aflp.mbr[,2])


coords.hd1.aflp.hull<-as.matrix(coords.hd1.aflp[,5:6])
coords.hd1.aflp.mbr<-MBR(coords.hd1.aflp.hull)
limits <- apply(coords.hd1.aflp.mbr, 2, function(x) c(min(x),max(x))) # Plotting limits
plot(ashape(coords.hd1.aflp.hull, alpha=1000), col="Gray", pch=20, 
     xlim=limits[,1], ylim=limits[,2])                # The hull
lines(coords.hd1.aflp.mbr, col="Blue", lwd=3)                         # The MBR
points(coords.hd1.aflp.hull, pch=19)    
polyarea(coords.hd1.aflp.mbr[,1], coords.hd1.aflp.mbr[,2])


coords.hd2.aflp.hull<-as.matrix(coords.hd2.aflp[,5:6])
coords.hd2.aflp.mbr<-MBR(coords.hd2.aflp.hull)
limits <- apply(coords.hd2.aflp.mbr, 2, function(x) c(min(x),max(x))) # Plotting limits
plot(ashape(coords.hd2.aflp.hull, alpha=1000), col="Gray", pch=20, 
     xlim=limits[,1], ylim=limits[,2])                # The hull
lines(coords.hd2.aflp.mbr, col="Blue", lwd=3)                         # The MBR
points(coords.hd2.aflp.hull, pch=19)    
polyarea(coords.hd2.aflp.mbr[,1], coords.hd2.aflp.mbr[,2])


coords.hd3.aflp.hull<-as.matrix(coords.hd3.aflp[,5:6])
coords.hd3.aflp.mbr<-MBR(coords.hd3.aflp.hull)
limits <- apply(coords.hd3.aflp.mbr, 2, function(x) c(min(x),max(x))) # Plotting limits
plot(ashape(coords.hd3.aflp.hull, alpha=1000), col="Gray", pch=20, 
     xlim=limits[,1], ylim=limits[,2])                # The hull
lines(coords.hd3.aflp.mbr, col="Blue", lwd=3)                         # The MBR
points(coords.hd3.aflp.hull, pch=19)    
polyarea(coords.hd3.aflp.mbr[,1], coords.hd3.aflp.mbr[,2])


coords.JL.aflp.u<-coords.JL.aflp[-38,]

coords.JL.aflp.hull<-as.matrix(coords.JL.aflp.u[,c("x","y")])
coords.JL.aflp.mbr<-MBR(coords.JL.aflp.hull)
limits <- apply(coords.JL.aflp.mbr, 2, function(x) c(min(x),max(x))) # Plotting limits
plot(ashape(coords.JL.aflp.hull, alpha=1000), col="Gray", pch=20, 
     xlim=limits[,1], ylim=limits[,2])                # The hull
lines(coords.JL.aflp.mbr, col="Blue", lwd=3)                         # The MBR
points(coords.JL.aflp.hull, pch=19)    
polyarea(coords.JL.aflp.mbr[,1], coords.JL.aflp.mbr[,2])

coords.CESAC02.aflp.u<-coords.CESAC02.aflp[-c(3,21),]
coords.CESAC02.aflp.hull<-as.matrix(coords.CESAC02.aflp.u[,c("x","y")])
coords.CESAC02.aflp.mbr<-MBR(coords.CESAC02.aflp.hull)
limits <- apply(coords.CESAC02.aflp.mbr, 2, function(x) c(min(x),max(x))) # Plotting limits
plot(ashape(coords.CESAC02.aflp.hull, alpha=1000), col="Gray", pch=20, 
     xlim=limits[,1], ylim=limits[,2])                # The hull
lines(coords.CESAC02.aflp.mbr, col="Blue", lwd=3)                         # The MBR
points(coords.CESAC02.aflp.hull, pch=19)    
polyarea(coords.CESAC02.aflp.mbr[,1], coords.CESAC02.aflp.mbr[,2])


coords.CESAC06.aflp.u<-coords.CESAC06.aflp[-32,]
coords.CESAC06.aflp.hull<-as.matrix(coords.CESAC06.aflp.u[,c("x","y")])
coords.CESAC06.aflp.mbr<-MBR(coords.CESAC06.aflp.hull)
limits <- apply(coords.CESAC06.aflp.mbr, 2, function(x) c(min(x),max(x))) # Plotting limits
plot(ashape(coords.CESAC06.aflp.hull, alpha=1000), col="Gray", pch=20, 
     xlim=limits[,1], ylim=limits[,2])                # The hull
lines(coords.CESAC06.aflp.mbr, col="Blue", lwd=3)                         # The MBR
points(coords.CESAC06.aflp.hull, pch=19)    
polyarea(coords.CESAC06.aflp.mbr[,1], coords.CESAC02.aflp.mbr[,2])


coords.S.aflp.hull<-as.matrix(coords.S.aflp[,c("x","y")])
coords.S.aflp.mbr<-MBR(coords.S.aflp.hull)
limits <- apply(coords.S.aflp.mbr, 2, function(x) c(min(x),max(x))) # Plotting limits
plot(ashape(coords.S.aflp.hull, alpha=1000), col="Gray", pch=20, 
     xlim=limits[,1], ylim=limits[,2])                # The hull
lines(coords.S.aflp.mbr, col="Blue", lwd=3)                         # The MBR
points(coords.S.aflp.hull, pch=19)    
#polyarea(coords.S.aflp.mbr[,1], coords.S.aflp.mbr[,2])/100
polyarea(coords.S.aflp.mbr[,1]/100, coords.S.aflp.mbr[,2]/100)


coords.str.aflp.hull<-as.matrix(coords.str.aflp[,c("x","y")])
coords.str.aflp.mbr<-MBR(coords.str.aflp.hull)
limits <- apply(coords.str.aflp.mbr, 2, function(x) c(min(x),max(x))) # Plotting limits
plot(ashape(coords.str.aflp.hull, alpha=1000), col="Gray", pch=20, 
     xlim=limits[,1], ylim=limits[,2])                # The hull
lines(coords.str.aflp.mbr, col="Blue", lwd=3)                         # The MBR
points(coords.str.aflp.hull, pch=19)    
polyarea(coords.str.aflp.mbr[,1], coords.str.aflp.mbr[,2])


coords.roch1.aflp.hull<-as.matrix(coords.roch1.aflp[,c("x","y")])
coords.roch1.aflp.mbr<-MBR(coords.roch1.aflp.hull)
limits <- apply(coords.roch1.aflp.mbr, 2, function(x) c(min(x),max(x))) # Plotting limits
plot(ashape(coords.roch1.aflp.hull, alpha=1000), col="Gray", pch=20, 
     xlim=limits[,1], ylim=limits[,2])                # The hull
lines(coords.roch1.aflp.mbr, col="Blue", lwd=3)                         # The MBR
points(coords.roch1.aflp.hull, pch=19)    
polyarea(coords.roch1.aflp.mbr[,1], coords.roch1.aflp.mbr[,2])


coords.roch2.aflp.hull<-as.matrix(coords.roch2.aflp[,c("x","y")])
coords.roch2.aflp.mbr<-MBR(coords.roch2.aflp.hull)
limits <- apply(coords.roch2.aflp.mbr, 2, function(x) c(min(x),max(x))) # Plotting limits
plot(ashape(coords.roch2.aflp.hull, alpha=1000), col="Gray", pch=20, 
     xlim=limits[,1], ylim=limits[,2])                # The hull
lines(coords.roch2.aflp.mbr, col="Blue", lwd=3)                         # The MBR
points(coords.roch2.aflp.hull, pch=19)    
polyarea(coords.roch2.aflp.mbr[,1], coords.roch2.aflp.mbr[,2])


#coords.roch3.aflp.hull<-as.matrix(coords.roch3.aflp[,c("x","y")])
#coords.roch3.aflp.mbr<-MBR(coords.roch3.aflp.hull)
#limits <- apply(coords.roch3.aflp.mbr, 2, function(x) c(min(x),max(x))) # Plotting limits
#plot(ashape(coords.roch3.aflp.hull, alpha=1000), col="Gray", pch=20, 
#     xlim=limits[,1], ylim=limits[,2])                # The hull
#lines(coords.roch3.aflp.mbr, col="Blue", lwd=3)                         # The MBR
#points(coords.roch3.aflp.hull, pch=19)    
#polyarea(coords.roch3.aflp.mbr[,1], coords.roch3.aflp.mbr[,2])


coords.d2.04.snp.hull<-as.matrix(coords.d2.04.snp[,c("x","y")])
coords.d2.04.snp.mbr<-MBR(coords.d2.04.snp.hull)
limits <- apply(coords.d2.04.snp.mbr, 2, function(x) c(min(x),max(x))) # Plotting limits
plot(ashape(coords.d2.04.snp.hull, alpha=1000), col="Gray", pch=20, 
     xlim=limits[,1], ylim=limits[,2])                # The hull
lines(coords.d2.04.snp.mbr, col="Blue", lwd=3)                         # The MBR
points(coords.d2.04.snp.hull, pch=19)    
polyarea(coords.d2.04.snp.mbr[,1], coords.d2.04.snp.mbr[,2])


coords.d2.14.snp.hull<-as.matrix(coords.d2.14.snp[,c("x","y")])
coords.d2.14.snp.mbr<-MBR(coords.d2.14.snp.hull)
limits <- apply(coords.d2.14.snp.mbr, 2, function(x) c(min(x),max(x))) # Plotting limits
plot(ashape(coords.d2.14.snp.hull, alpha=1000), col="Gray", pch=20, 
     xlim=limits[,1], ylim=limits[,2])                # The hull
lines(coords.d2.14.snp.mbr, col="Blue", lwd=3)                         # The MBR
points(coords.d2.14.snp.hull, pch=19)    
polyarea(coords.d2.14.snp.mbr[,1], coords.d2.14.snp.mbr[,2])



coords.d2.15.snp.hull<-as.matrix(coords.d2.15.snp[,c("x","y")])
coords.d2.15.snp.mbr<-MBR(coords.d2.15.snp.hull)
limits <- apply(coords.d2.15.snp.mbr, 2, function(x) c(min(x),max(x))) # Plotting limits
plot(ashape(coords.d2.15.snp.hull, alpha=1000), col="Gray", pch=20, 
     xlim=limits[,1], ylim=limits[,2])                # The hull
lines(coords.d2.15.snp.mbr, col="Blue", lwd=3)                         # The MBR
points(coords.d2.15.snp.hull, pch=19)    
polyarea(coords.d2.15.snp.mbr[,1], coords.d2.15.snp.mbr[,2])


coords.d3.04.snp.hull<-as.matrix(coords.d3.04.snp[,c("x","y")])
coords.d3.04.snp.mbr<-MBR(coords.d3.04.snp.hull)
limits <- apply(coords.d3.04.snp.mbr, 2, function(x) c(min(x),max(x))) # Plotting limits
plot(ashape(coords.d3.04.snp.hull, alpha=1000), col="Gray", pch=20, 
     xlim=limits[,1], ylim=limits[,2])                # The hull
lines(coords.d3.04.snp.mbr, col="Blue", lwd=3)                         # The MBR
points(coords.d3.04.snp.hull, pch=19)    
polyarea(coords.d3.04.snp.mbr[,1], coords.d3.04.snp.mbr[,2])


coords.d3.14.snp.hull<-as.matrix(coords.d3.14.snp[,c("x","y")])
coords.d3.14.snp.mbr<-MBR(coords.d3.14.snp.hull)
limits <- apply(coords.d3.14.snp.mbr, 2, function(x) c(min(x),max(x))) # Plotting limits
plot(ashape(coords.d3.14.snp.hull, alpha=1000), col="Gray", pch=20, 
     xlim=limits[,1], ylim=limits[,2])                # The hull
lines(coords.d3.14.snp.mbr, col="Blue", lwd=3)                         # The MBR
points(coords.d3.14.snp.hull, pch=19)    
polyarea(coords.d3.14.snp.mbr[,1], coords.d3.14.snp.mbr[,2])



coords.d3.15.snp.hull<-as.matrix(coords.d3.15.snp[,c("x","y")])
coords.d3.15.snp.mbr<-MBR(coords.d3.15.snp.hull)
limits <- apply(coords.d3.15.snp.mbr, 2, function(x) c(min(x),max(x))) # Plotting limits
plot(ashape(coords.d3.15.snp.hull, alpha=1000), col="Gray", pch=20, 
     xlim=limits[,1], ylim=limits[,2])                # The hull
lines(coords.d3.15.snp.mbr, col="Blue", lwd=3)                         # The MBR
points(coords.d3.15.snp.hull, pch=19)    
polyarea(coords.d3.15.snp.mbr[,1], coords.d3.15.snp.mbr[,2])



coords.mira.hull<-as.matrix(coords.mira[,c("x","y")])
coords.mira.mbr<-MBR(coords.mira.hull)
limits <- apply(coords.mira.mbr, 2, function(x) c(min(x),max(x))) # Plotting limits
plot(ashape(coords.mira.hull, alpha=1000), col="Gray", pch=20, 
     xlim=limits[,1], ylim=limits[,2])                # The hull
lines(coords.mira.mbr, col="Blue", lwd=3)                         # The MBR
points(coords.mira.hull, pch=19)    
polyarea(coords.mira.mbr[,1], coords.mira.mbr[,2])

# # # # # # # # # # # Combo


cords.HD12.hull<-as.matrix(cords.HD12[,c("x.trans","y.trans")])
cords.HD12.mbr<-MBR(cords.HD12.hull)
limits <- apply(cords.HD12.mbr, 2, function(x) c(min(x),max(x))) # Plotting limits
plot(ashape(cords.HD12.hull, alpha=1000), col="Gray", pch=20, 
     xlim=limits[,1], ylim=limits[,2])                # The hull
lines(cords.HD12.mbr, col="Blue", lwd=3)                         # The MBR
points(cords.HD12.hull, pch=19)    
polyarea(cords.HD12.mbr[,1], cords.HD12.mbr[,2])


cords.HD23.hull<-as.matrix(cords.HD23[,c("x.trans","y.trans")])
cords.HD23.mbr<-MBR(cords.HD23.hull)
limits <- apply(cords.HD23.mbr, 2, function(x) c(min(x),max(x))) # Plotting limits
plot(ashape(cords.HD23.hull, alpha=1000), col="Gray", pch=20, 
     xlim=limits[,1], ylim=limits[,2])                # The hull
lines(cords.HD23.mbr, col="Blue", lwd=3)                         # The MBR
points(cords.HD23.hull, pch=19)    
polyarea(cords.HD23.mbr[,1], cords.HD23.mbr[,2])


cords.HD13.hull<-as.matrix(cords.HD13[,c("x.trans","y.trans")])
cords.HD13.mbr<-MBR(cords.HD13.hull)
limits <- apply(cords.HD13.mbr, 2, function(x) c(min(x),max(x))) # Plotting limits
plot(ashape(cords.HD13.hull, alpha=1000), col="Gray", pch=20, 
     xlim=limits[,1], ylim=limits[,2])                # The hull
lines(cords.HD13.mbr, col="Blue", lwd=3)                         # The MBR
points(cords.HD13.hull, pch=19)    
polyarea(cords.HD13.mbr[,1], cords.HD13.mbr[,2])


coords.HD.hull<-as.matrix(coords.HD[,c("x.trans","y.trans")])
coords.HD.mbr<-MBR(coords.HD.hull)
limits <- apply(coords.HD.mbr, 2, function(x) c(min(x),max(x))) # Plotting limits
plot(ashape(coords.HD.hull, alpha=1000), col="Gray", pch=20, 
     xlim=limits[,1], ylim=limits[,2])                # The hull
lines(coords.HD.mbr, col="Blue", lwd=3)                         # The MBR
points(coords.HD.hull, pch=19)    
polyarea(coords.HD.mbr[,1], coords.HD.mbr[,2])


coords.roch12.hull<-as.matrix(coords.roch12[,c("x.trans","y.trans")])
coords.roch12.mbr<-MBR(coords.roch12.hull)
limits <- apply(coords.roch12.mbr, 2, function(x) c(min(x),max(x))) # Plotting limits
plot(ashape(coords.roch12.hull, alpha=1000), col="Gray", pch=20, 
     xlim=limits[,1], ylim=limits[,2])                # The hull
lines(coords.roch12.mbr, col="Blue", lwd=3)                         # The MBR
points(coords.roch12.hull, pch=19)    
polyarea(coords.roch12.mbr[,1], coords.roch12.mbr[,2])


coords.roch23.hull<-as.matrix(coords.roch23[,c("x.trans","y.trans")])
coords.roch23.mbr<-MBR(coords.roch23.hull)
limits <- apply(coords.roch23.mbr, 2, function(x) c(min(x),max(x))) # Plotting limits
plot(ashape(coords.roch23.hull, alpha=1000), col="Gray", pch=20, 
     xlim=limits[,1], ylim=limits[,2])                # The hull
lines(coords.roch23.mbr, col="Blue", lwd=3)                         # The MBR
points(coords.roch23.hull, pch=19)    
polyarea(coords.roch23.mbr[,1], coords.roch23.mbr[,2])


coords.roch13.hull<-as.matrix(coords.roch13[,c("x.trans","y.trans")])
coords.roch13.mbr<-MBR(coords.roch13.hull)
limits <- apply(coords.roch13.mbr, 2, function(x) c(min(x),max(x))) # Plotting limits
plot(ashape(coords.roch13.hull, alpha=1000), col="Gray", pch=20, 
     xlim=limits[,1], ylim=limits[,2])                # The hull
lines(coords.roch13.mbr, col="Blue", lwd=3)                         # The MBR
points(coords.roch13.hull, pch=19)    
polyarea(coords.roch13.mbr[,1], coords.roch13.mbr[,2])


coords.roch.hull<-as.matrix(coords.roch[,c("x.trans","y.trans")])
coords.roch.mbr<-MBR(coords.roch.hull)
limits <- apply(coords.roch.mbr, 2, function(x) c(min(x),max(x))) # Plotting limits
plot(ashape(coords.roch.hull, alpha=1000), col="Gray", pch=20, 
     xlim=limits[,1], ylim=limits[,2])                # The hull
lines(coords.roch.mbr, col="Blue", lwd=3)                         # The MBR
points(coords.roch.hull, pch=19)    
polyarea(coords.roch.mbr[,1], coords.roch.mbr[,2])


coords.trans.dr12.hull<-as.matrix(coords.trans.dr12[,c("x.trans","y.trans")])
coords.trans.dr12.mbr<-MBR(coords.trans.dr12.hull)
limits <- apply(coords.trans.dr12.mbr, 2, function(x) c(min(x),max(x))) # Plotting limits
plot(ashape(coords.trans.dr12.hull, alpha=1000), col="Gray", pch=20, 
     xlim=limits[,1], ylim=limits[,2])                # The hull
lines(coords.trans.dr12.mbr, col="Blue", lwd=3)                         # The MBR
points(coords.trans.dr12.hull, pch=19)    
polyarea(coords.trans.dr12.mbr[,1], coords.trans.dr12.mbr[,2])


coords.trans.dr13.hull<-as.matrix(coords.trans.dr13[,c("x.trans","y.trans")])
coords.trans.dr13.mbr<-MBR(coords.trans.dr13.hull)
limits <- apply(coords.trans.dr13.mbr, 2, function(x) c(min(x),max(x))) # Plotting limits
plot(ashape(coords.trans.dr13.hull, alpha=1000), col="Gray", pch=20, 
     xlim=limits[,1], ylim=limits[,2])                # The hull
lines(coords.trans.dr13.mbr, col="Blue", lwd=3)                         # The MBR
points(coords.trans.dr13.hull, pch=19)    
polyarea(coords.trans.dr13.mbr[,1], coords.trans.dr13.mbr[,2])


coords.trans.dr14.hull<-as.matrix(coords.trans.dr14[,c("x.trans","y.trans")])
coords.trans.dr14.mbr<-MBR(coords.trans.dr14.hull)
limits <- apply(coords.trans.dr14.mbr, 2, function(x) c(min(x),max(x))) # Plotting limits
plot(ashape(coords.trans.dr14.hull, alpha=1000), col="Gray", pch=20, 
     xlim=limits[,1], ylim=limits[,2])                # The hull
lines(coords.trans.dr14.mbr, col="Blue", lwd=3)                         # The MBR
points(coords.trans.dr14.hull, pch=19)    
polyarea(coords.trans.dr14.mbr[,1], coords.trans.dr14.mbr[,2])


coords.trans.dr23.hull<-as.matrix(coords.trans.dr23[,c("x.trans","y.trans")])
coords.trans.dr23.mbr<-MBR(coords.trans.dr23.hull)
limits <- apply(coords.trans.dr23.mbr, 2, function(x) c(min(x),max(x))) # Plotting limits
plot(ashape(coords.trans.dr23.hull, alpha=1000), col="Gray", pch=20, 
     xlim=limits[,1], ylim=limits[,2])                # The hull
lines(coords.trans.dr23.mbr, col="Blue", lwd=3)                         # The MBR
points(coords.trans.dr23.hull, pch=19)    
polyarea(coords.trans.dr23.mbr[,1], coords.trans.dr23.mbr[,2])


coords.trans.dr24.hull<-as.matrix(coords.trans.dr24[,c("x.trans","y.trans")])
coords.trans.dr24.mbr<-MBR(coords.trans.dr24.hull)
limits <- apply(coords.trans.dr24.mbr, 2, function(x) c(min(x),max(x))) # Plotting limits
plot(ashape(coords.trans.dr24.hull, alpha=1000), col="Gray", pch=20, 
     xlim=limits[,1], ylim=limits[,2])                # The hull
lines(coords.trans.dr24.mbr, col="Blue", lwd=3)                         # The MBR
points(coords.trans.dr24.hull, pch=19)    
polyarea(coords.trans.dr24.mbr[,1], coords.trans.dr24.mbr[,2])


coords.trans.dr34.hull<-as.matrix(coords.trans.dr34[,c("x.trans","y.trans")])
coords.trans.dr34.mbr<-MBR(coords.trans.dr34.hull)
limits <- apply(coords.trans.dr34.mbr, 2, function(x) c(min(x),max(x))) # Plotting limits
plot(ashape(coords.trans.dr34.hull, alpha=1000), col="Gray", pch=20, 
     xlim=limits[,1], ylim=limits[,2])                # The hull
lines(coords.trans.dr34.mbr, col="Blue", lwd=3)                         # The MBR
points(coords.trans.dr34.hull, pch=19)    
polyarea(coords.trans.dr34.mbr[,1], coords.trans.dr34.mbr[,2])


coords.trans.dr.hull<-as.matrix(coords.trans.dr[,c("x.trans","y.trans")])
coords.trans.dr.mbr<-MBR(coords.trans.dr.hull)
limits <- apply(coords.trans.dr.mbr, 2, function(x) c(min(x),max(x))) # Plotting limits
plot(ashape(coords.trans.dr.hull, alpha=1000), col="Gray", pch=20, 
     xlim=limits[,1], ylim=limits[,2])                # The hull
lines(coords.trans.dr.mbr, col="Blue", lwd=3)                         # The MBR
points(coords.trans.dr.hull, pch=19)    
polyarea(coords.trans.dr.mbr[,1], coords.trans.dr.mbr[,2])





coords.trans.dr23.04.snp.hull<-as.matrix(coords.trans.dr23.04.snp[,c("x.trans","y.trans")])
coords.trans.dr23.04.snp.mbr<-MBR(coords.trans.dr23.04.snp.hull)
limits <- apply(coords.trans.dr23.04.snp.mbr, 2, function(x) c(min(x),max(x))) # Plotting limits
plot(ashape(coords.trans.dr23.04.snp.hull, alpha=1000), col="Gray", pch=20, 
     xlim=limits[,1], ylim=limits[,2])                # The hull
lines(coords.trans.dr23.04.snp.mbr, col="Blue", lwd=3)                         # The MBR
points(coords.trans.dr23.04.snp.hull, pch=19)    
polyarea(coords.trans.dr23.04.snp.mbr[,1], coords.trans.dr23.04.snp.mbr[,2])


coords.trans.dr23.14.snp.hull<-as.matrix(coords.trans.dr23.14.snp[,c("x.trans","y.trans")])
coords.trans.dr23.14.snp.mbr<-MBR(coords.trans.dr23.14.snp.hull)
limits <- apply(coords.trans.dr23.14.snp.mbr, 2, function(x) c(min(x),max(x))) # Plotting limits
plot(ashape(coords.trans.dr23.14.snp.hull, alpha=1000), col="Gray", pch=20, 
     xlim=limits[,1], ylim=limits[,2])                # The hull
lines(coords.trans.dr23.14.snp.mbr, col="Blue", lwd=3)                         # The MBR
points(coords.trans.dr23.14.snp.hull, pch=19)    
polyarea(coords.trans.dr23.14.snp.mbr[,1], coords.trans.dr23.14.snp.mbr[,2])



coords.trans.dr23.15.snp.hull<-as.matrix(coords.trans.dr23.15.snp[,c("x.trans","y.trans")])
coords.trans.dr23.15.snp.mbr<-MBR(coords.trans.dr23.15.snp.hull)
limits <- apply(coords.trans.dr23.15.snp.mbr, 2, function(x) c(min(x),max(x))) # Plotting limits
plot(ashape(coords.trans.dr23.15.snp.hull, alpha=1000), col="Gray", pch=20, 
     xlim=limits[,1], ylim=limits[,2])                # The hull
lines(coords.trans.dr23.15.snp.mbr, col="Blue", lwd=3)                         # The MBR
points(coords.trans.dr23.15.snp.hull, pch=19)    
polyarea(coords.trans.dr23.15.snp.mbr[,1], coords.trans.dr23.15.snp.mbr[,2])


###################################
library(ggplot2)
setwd("~/Dropbox/GENETS_finishing")
mantel<-read.csv("mantel.n.area.4R.csv", head=T, stringsAsFactors = F)

mantel.sig<-mantel[mantel$signif=="sig",]

#ggplot(mantel[c(2:6,8:23),], aes(x=area, y=r))+geom_point()+geom_smooth(method='lm')
#ggplot(mantel[c(2:23),], aes(x=mushies, y=r))+geom_point()+geom_smooth(method='lm')

#ggplot(mantel.sig[4:9,], aes(x=area, y=r, size=mushies))+geom_point()+geom_smooth(method='lm')


ggplot(mantel, aes(x=area, y=r, size=mushies, color=signif))+geom_point()+geom_smooth(method='lm')
cor.test(mantel$r, mantel$area, method="s")
cor.test(mantel[mantel$signif=="sig",]$r, mantel[mantel$signif=="sig",]$area, method="s")


#only within pops
ggplot(mantel[1:18,], aes(x=area, y=r, size=mushies, color=signif))+geom_point()+geom_smooth(method='lm')

summary(lm(r~area, mantel[1:18,]))

summary(lm(r~area, mantel[mantel$signif == "sig",][1:18,]))
summary(lm(r~area, mantel[mantel$signif == "in",][1:18,]))

#only between pops
ggplot(mantel[19:nrow(mantel),], aes(x=area, y=r, size=mushies, color=signif))+geom_point()+geom_smooth(method='lm')

summary(lm(r~area, mantel[22:39,]))
summary(lm(r~area, mantel[mantel$signif =="sig",][22:39,]))
summary(lm(r~area, mantel[mantel$signif =="in",][22:39,]))

#Serbia and JL are very large, so just see what it looks like with them removed
ggplot(mantel[c(2:5,7:18),], aes(x=area, y=r, size=mushies, color=signif))+geom_point()+geom_smooth(method='lm')
summary(lm(r~area, mantel[c(2:5,7:18),]))


library(ggrepel)
ggplot(mantel, aes(area, r, label = pop)) +
  geom_point(color = ifelse(mantel$area > 100000, "grey50", "red")) +
  geom_text_repel()




ggplot(mantel, aes(area, r, label = pop, shape=signif, size=mushies)) +
  geom_text_repel(
    data          = subset(mantel, area > 10000),
    nudge_x       = subset(mantel, area > 10000)$area,
    segment.size  = 0.2,
    segment.color = "grey50",
    direction     = "x"
  ) +
  geom_point(color = ifelse(mantel$area > 10000, "red", "black")) +
  scale_x_continuous(expand = c(0, 12500)) +
  scale_y_continuous(limits = c(-1, 1))+geom_smooth(method="lm")

#regardless of r mantel significance
summary(lm(r~area, mantel))

summary(lm(r~area, mantel[mantel$signif == "in",]))

summary(lm(r~area, mantel[mantel$signif == "sig",]))


mantel.noS<-mantel[-1,]
ggplot(mantel.noS, aes(x=area, y=r, size=mushies, color=signif))+geom_point()+geom_smooth(method='lm')

ggplot(mantel.noS, aes(area, r, label = pop)) +
  geom_point(color = ifelse(mantel.noS$area > 100000, "grey50", "red")) +
  geom_text_repel()




ggplot(mantel.noS, aes(area, r, label = pop, shape=signif, size=mushies)) +
  geom_text_repel(
    data          = subset(mantel.noS, area > 10000),
    nudge_x       = subset(mantel.noS, area > 10000)$area,
    segment.size  = 0.2,
    segment.color = "grey50",
    direction     = "x"
  ) +
  geom_point(color = ifelse(mantel.noS$area > 10000, "red", "black")) +
  scale_x_continuous(expand = c(0, 12500)) +
  scale_y_continuous(limits = c(-1, 1))+geom_smooth(method="lm")


summary(lm(r~area, mantel[mantel.noS$signif == "in",]))

summary(lm(r~area, mantel[mantel.noS$signif == "sig",]))


ggplot(mantel, aes(x=log(area), y=r, color=signif))+geom_point()+geom_smooth(method="lm")
summary(lm(r~log(area), mantel[mantel$signif == "sig",]))
summary(lm(r~log(area), mantel[mantel$signif == "in",]))

