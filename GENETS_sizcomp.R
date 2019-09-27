####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
#                                         GENET SIZE COMPARISON
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
setwd("~/Dropbox/GENETS_finishing")
coords<-read.csv("ALL_COORDS_geo_AMANITABASED.csv", head=T, stringsAsFactors = F) #just so that I had the right labels

coords.d1<-coords[coords$pop =="Drake1",]
coords.d1.aflp<-coords.d1

max(dist(coords.d1.aflp[coords.d1.aflp$genet.1 =="d1A",][,5:6]))
max(dist(coords.d1.aflp[coords.d1.aflp$genet.1 =="d1B",][,5:6]))
max(dist(coords.d1.aflp[coords.d1.aflp$genet.1 =="d1C",][,5:6]))
max(dist(coords.d1.aflp[coords.d1.aflp$genet.1 =="d1D",][,5:6]))
max(dist(coords.d1.aflp[coords.d1.aflp$genet.1 =="d1E",][,5:6]))
max(dist(coords.d1.aflp[coords.d1.aflp$genet.1 =="d1F",][,5:6]))




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

max(dist(coords.d3.04.aflp[coords.d3.04.aflp$shape =="tri",][,5:6]))


coords.d3.14<-coords.d3[coords.d3$year =="Sample_2014",]
coords.d3.14.snp<-coords.d3.14

coords.d3.15<-coords.d3[coords.d3$year =="Sample_2015",]
coords.d3.15.snp<-coords.d3.15


coords.d4<-coords[coords$pop =="Drake4",]
coords.d4.aflp<-coords.d4[-2,]

max(dist(coords.d4.aflp[coords.d4.aflp$genet.1 =="d4.A",][,5:6]))


coords.hd1<-coords[coords$pop =="HD1",]
coords.hd1.aflp<-coords.hd1
coords.hd2<-coords[coords$pop =="HD2",]
coords.hd2.aflp<-coords.hd2[-4,]
coords.hd3<-coords[coords$pop =="HD3",]
coords.hd3.aflp<-coords.hd3

max(dist(coords.hd2.aflp[coords.hd2.aflp$genet.1 =="A",][,5:6]))
max(dist(coords.hd2.aflp[coords.hd2.aflp$genet.1 =="B",][,5:6]))

max(dist(coords.hd3.aflp[coords.hd3.aflp$genet.1 =="A",][,5:6]))


coords.mira<-coords[coords$pop =="mira",]
coords.mira.snp<-coords.mira
coords.agraria<-coords[coords$pop =="agraria",]
coords.agraria.snp<-coords.agraria

coords.JL<-coords[coords$pop =="JL",]
coords.JL.aflp<-coords.JL

max(dist(coords.JL.aflp[coords.JL.aflp$genet.1 =="A",][,5:6]))
max(dist(coords.JL.aflp[coords.JL.aflp$genet.1 =="B",][,5:6]))
max(dist(coords.JL.aflp[coords.JL.aflp$genet.1 =="C",][,5:6]))
max(dist(coords.JL.aflp[coords.JL.aflp$genet.1 =="D",][,5:6]))
max(dist(coords.JL.aflp[coords.JL.aflp$genet.1 =="E",][,5:6]))
max(dist(coords.JL.aflp[coords.JL.aflp$genet.1 =="F",][,5:6]))
max(dist(coords.JL.aflp[coords.JL.aflp$genet.1 =="G",][,5:6]))
max(dist(coords.JL.aflp[coords.JL.aflp$genet.1 =="H",][,5:6]))


coords.CESAC02<-coords[coords$pop =="CESAC02",]
coords.CESAC02.aflp<-coords.CESAC02
coords.CESAC06<-coords[coords$pop =="CESAC06",]
coords.CESAC06.aflp<-coords.CESAC06

max(dist(coords.CESAC02.aflp[coords.CESAC02.aflp$genet.1 =="D",][,5:6]))

max(dist(coords.CESAC06.aflp[coords.CESAC06.aflp$genet.1 =="A",][,5:6]))
max(dist(coords.CESAC06.aflp[coords.CESAC06.aflp$genet.1 =="B",][,5:6]))
max(dist(coords.CESAC06.aflp[coords.CESAC06.aflp$genet.1 =="C",][,5:6]))
max(dist(coords.CESAC06.aflp[coords.CESAC06.aflp$genet.1 =="D",][,5:6]))
max(dist(coords.CESAC06.aflp[coords.CESAC06.aflp$genet.1 =="E",][,5:6]))

coords.S<-coords[coords$pop =="S",]
coords.S.aflp<-coords.S # note S 11-14 do not have coords, need to remove from genind

max(dist(coords.S.aflp[coords.S.aflp$genet.1 =="SA",][,5:6]))
max(dist(coords.S.aflp[coords.S.aflp$genet.1 =="SB",][,5:6]))
max(dist(coords.S.aflp[coords.S.aflp$genet.1 =="SC",][,5:6]))


coords.str<-coords[coords$pop =="str",]
coords.str.aflp<-coords.str[-c(2,6,8,12,15,16),]

max(dist(coords.str.aflp[coords.str.aflp$genet.1 =="A",][,5:6]))


coords.roch<-coords[coords$pop =="roch",]
coords.roch1<-coords.roch[1:13,]
coords.roch1.aflp<-coords.roch1
coords.roch2<-coords.roch[14:20,]
coords.roch2.aflp<-coords.roch2[-4,] #need to remove 2_6b from genind file
coords.roch3<-coords.roch[21:30,]
coords.roch3.aflp<-coords.roch3[3:4,]

max(dist(coords.roch1.aflp[coords.roch1.aflp$genet.1 =="A1",][,5:6]))
max(dist(coords.roch2.aflp[coords.roch2.aflp$genet.1 =="A2",][,5:6]))


aflp.max<-c(1.48,2.08,0.81,2.00,1.73,4.2,0.48,2.26,5.14,0.69,2.34,0.72,7.51,0.19,4.08,095,0.64,0.40,0.25,1.48,1.12,0,1.66,0.84,2.73,3.20,0.09,2.97,2.20,7.39)
median(aflp.max)
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
getmode(aflp.max)

max.EU<-c(0.25,1.48,1.12,0,1.66,0.84,2.73,3.20,0.09)
max.EC<-c(2.34,0.72,7.51,0.19,4.08,0.95,0.64,0.40,2.97,2.20,7.39)
max.C<-c(1.48,2.08,0.81,0.18,2.00,1.73,4.2,0.48,2.26,5.14,0.69)

getmode(max.EU)
getmode(max.EC)
getmode(max.C)



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


##################################


maxgen<-read.csv("maxgen.lengths.2018.csv", head=T, stringsAsFactors = F)

ggplot(maxgen, aes(x=range, y=max.dist, fill=range))+geom_boxplot()
ggplot(maxgen, aes(x=cont, y=max.dist, fill=cont))+geom_boxplot()

t.test(max.dist ~ cont, data = maxgen)

t.test(max.dist ~ range, data = maxgen[maxgen$range == c("CA", "EC"),])
t.test(max.dist ~ range, data = maxgen[maxgen$range == c("CA", "EU"),])
t.test(max.dist ~ range, data = maxgen[maxgen$range == c("EU", "EC"),])

ggplot(maxgen, aes(x=max.dist, fill=range))+geom_histogram()

wilcox.test(maxgen$max.dist[maxgen$range=="CA"],maxgen$max.dist[maxgen$range=="EU"])
wilcox.test(maxgen$max.dist[maxgen$range=="EU"],maxgen$max.dist[maxgen$range=="EC"])
wilcox.test(maxgen$max.dist[maxgen$range=="EU"],maxgen$max.dist[maxgen$range=="EC"|maxgen$range=="CA"])

#

upper.tri(as.matrix(dist(coords.d1.aflp[coords.d1.aflp$genet.1 =="d1A",][,5:6])), diag=T)

d1A[upper.tri(d1A), diag=T] <- "rem"

d1A<-as.matrix(dist(coords.d1.aflp[coords.d1.aflp$genet.1 =="d1A",][,5:6]))
d1B<-as.matrix(dist(coords.d1.aflp[coords.d1.aflp$genet.1 =="d1B",][,5:6]))
d1C<-as.matrix(dist(coords.d1.aflp[coords.d1.aflp$genet.1 =="d1C",][,5:6]))
d1D<-as.matrix(dist(coords.d1.aflp[coords.d1.aflp$genet.1 =="d1D",][,5:6]))
d1E<-as.matrix(dist(coords.d1.aflp[coords.d1.aflp$genet.1 =="d1E",][,5:6]))
d1F<-as.matrix(dist(coords.d1.aflp[coords.d1.aflp$genet.1 =="d1F",][,5:6]))


d3tri<-as.matrix(dist(coords.d3.04.aflp[coords.d3.04.aflp$shape =="tri",][,5:6]))

d4A<-as.matrix(dist(coords.d4.aflp[coords.d4.aflp$genet.1 =="d4.A",][,5:6]))

hd1A<-as.matrix(dist(coords.hd2.aflp[coords.hd2.aflp$genet.1 =="A",][,5:6]))
hd1B<-as.matrix(dist(coords.hd2.aflp[coords.hd2.aflp$genet.1 =="B",][,5:6]))
hd2A<-as.matrix(dist(coords.hd3.aflp[coords.hd3.aflp$genet.1 =="A",][,5:6]))



JLA<-as.matrix(dist(coords.JL.aflp[coords.JL.aflp$genet.1 =="A",][,5:6]))
JLB<-as.matrix(dist(coords.JL.aflp[coords.JL.aflp$genet.1 =="B",][,5:6]))
JLC<-as.matrix(dist(coords.JL.aflp[coords.JL.aflp$genet.1 =="C",][,5:6]))
JLD<-as.matrix(dist(coords.JL.aflp[coords.JL.aflp$genet.1 =="D",][,5:6]))
JLE<-as.matrix(dist(coords.JL.aflp[coords.JL.aflp$genet.1 =="E",][,5:6]))
JLF<-as.matrix(dist(coords.JL.aflp[coords.JL.aflp$genet.1 =="F",][,5:6]))
JLG<-as.matrix(dist(coords.JL.aflp[coords.JL.aflp$genet.1 =="G",][,5:6]))
JLH<-as.matrix(dist(coords.JL.aflp[coords.JL.aflp$genet.1 =="H",][,5:6]))




ce02D<-as.matrix(dist(coords.CESAC02.aflp[coords.CESAC02.aflp$genet.1 =="D",][,5:6]))

ce06A<-as.matrix(dist(coords.CESAC06.aflp[coords.CESAC06.aflp$genet.1 =="A",][,5:6]))
ce06B<-as.matrix(dist(coords.CESAC06.aflp[coords.CESAC06.aflp$genet.1 =="B",][,5:6]))
ce06C<-as.matrix(dist(coords.CESAC06.aflp[coords.CESAC06.aflp$genet.1 =="C",][,5:6]))
ce06D<-as.matrix(dist(coords.CESAC06.aflp[coords.CESAC06.aflp$genet.1 =="D",][,5:6]))
ce06E<-as.matrix(dist(coords.CESAC06.aflp[coords.CESAC06.aflp$genet.1 =="E",][,5:6]))



SA<-as.matrix(dist(coords.S.aflp[coords.S.aflp$genet.1 =="SA",][,5:6]))/100
SB<-as.matrix(dist(coords.S.aflp[coords.S.aflp$genet.1 =="SB",][,5:6]))/100
SC<-as.matrix(dist(coords.S.aflp[coords.S.aflp$genet.1 =="SC",][,5:6]))/100



strA<-as.matrix(dist(coords.str.aflp[coords.str.aflp$genet.1 =="A",][,5:6]))


roch1A1<-as.matrix(dist(coords.roch1.aflp[coords.roch1.aflp$genet.1 =="A1",][,5:6]))
roch2A2<-as.matrix(dist(coords.roch2.aflp[coords.roch2.aflp$genet.1 =="A2",][,5:6]))

library(reshape2)
d1A.mlt<-melt(d1A)
d1B.mlt<-melt(d1B)
d1C.mlt<-melt(d1C)
d1D.mlt<-melt(d1D)
d1E.mlt<-melt(d1E)
d1F.mlt<-melt(d1F)
d3tri.mlt<-melt(d3tri)
d4A.mlt<-melt(d4A)
hd1A.mlt<-melt(hd1A)
hd1B.mlt<-melt(hd1B)
hd2A.mlt<-melt(hd2A)
JLA.mlt<-melt(JLA)
JLB.mlt<-melt(JLB)
JLC.mlt<-melt(JLC)
JLD.mlt<-melt(JLD)
JLE.mlt<-melt(JLE)
JLF.mlt<-melt(JLF)
JLG.mlt<-melt(JLG)
JLH.mlt<-melt(JLH)
ce02D.mlt<-melt(ce02D)
ce06A.mlt<-melt(ce06A)
ce06B.mlt<-melt(ce06B)
ce06C.mlt<-melt(ce06C)
ce06D.mlt<-melt(ce06D)
ce06E.mlt<-melt(ce06E)
SA.mlt<-melt(SA)
SB.mlt<-melt(SB)
SC.mlt<-melt(SC)
strA.mlt<-melt(strA)
roch1A1.mlt<-melt(roch1A1)
roch2A2.mlt<-melt(roch2A2)


d1A.mlt$lab<-"d1A"
d1B.mlt$lab<-"d1B"
d1C.mlt$lab<-"d1C"
d1D.mlt$lab<-"d1D"
d1E.mlt$lab<-"d1E"
d1F.mlt$lab<-"d1F"
d3tri.mlt$lab<-"d3tri"
d4A.mlt$lab<-"d4A"
hd1A.mlt$lab<-"hd1A"
hd1B.mlt$lab<-"hd1B"
hd2A.mlt$lab<-"hd2A"
JLA.mlt$lab<-"JLA"
JLB.mlt$lab<-"JLB"
JLC.mlt$lab<-"JLC"
JLD.mlt$lab<-"JLD"
JLE.mlt$lab<-"JLE"
JLF.mlt$lab<-"JLF"
JLG.mlt$lab<-"JLG"
JLH.mlt$lab<-"JLH"
ce02D.mlt$lab<-"ce02D"
ce06A.mlt$lab<-"ce06A"
ce06B.mlt$lab<-"ce06B"
ce06C.mlt$lab<-"ce06C"
ce06D.mlt$lab<-"ce06D"
ce06E.mlt$lab<-"ce06E"
SA.mlt$lab<-"SA"
SB.mlt$lab<-"SB"
SC.mlt$lab<-"SC"
strA.mlt$lab<-"strA"
roch1A1.mlt$lab<-"roch1A1"
roch2A2.mlt$lab<-"roch2A2"



maxgendists<-rbind(d1A.mlt,
                   d1B.mlt,
                   d1C.mlt,
                   d1D.mlt,
                   d1E.mlt,
                   d1F.mlt,
                   d3tri.mlt,
                   d4A.mlt,
                   hd1A.mlt,
                   hd1B.mlt,
                   hd2A.mlt,
                   JLA.mlt,
                   JLB.mlt,
                   JLC.mlt,
                   JLD.mlt,
                   JLE.mlt,
                   JLF.mlt,
                   JLG.mlt,
                   JLH.mlt,
                   ce02D.mlt,
                   ce06A.mlt,
                   ce06B.mlt,
                   ce06C.mlt,
                   ce06D.mlt,
                   ce06E.mlt,
                   SA.mlt,
                   SB.mlt,
                   SC.mlt,
                   strA.mlt,
                   roch1A1.mlt,
                   roch2A2.mlt)


library("dplyr")
#removes duplicates conditional on column 3 and 4 matching together
test<-distinct(maxgendists[,3:4])

clean.gendists<-test[test$value !=0,]
clean.gendists$row<-1:nrow(clean.gendists)

conts<-c(
  replicate(75, "US"),
  replicate(31, "EU"),
  replicate(10, "US")
)


range<-c(
  replicate(37, "CA"),
  replicate(38, "EC"),
  replicate(31, "EU"),
  replicate(10, "EC")
)

gendists<-cbind(clean.gendists, conts, range)

ggplot(gendists, aes(x=conts, y=value, fill=conts))+geom_boxplot()
ggplot(gendists, aes(x=range, y=value, fill=range))+geom_boxplot()




t.test(value ~ conts, data = gendists)

t.test(value ~ range, data = gendists[gendists$range == c("CA", "EC"),])
t.test(value ~ range, data = gendists[gendists$range == c("CA", "EU"),])
t.test(value ~ range, data = gendists[gendists$range == c("EU", "EC"),])
