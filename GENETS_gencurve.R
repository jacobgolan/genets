
############################################################################################
############################################################################################
############################################################################################
##############################          GENCURVES           ################################
############################################################################################
############################################################################################
############################################################################################
setwd("~/Dropbox/finalizing_genets")

library(adegenet)
library(poppr)
library(phangorn)
library(dartR)

#create gencurve objects
m6.1<-as.genind(as.matrix(read.PLINK("HESS_fix/gencurve/vcf.MMDP60.maf0.006.minQ30.mx01.raw", 
                                     map.file = "HESS_fix/gencurve/vcf.MMDP60.maf0.006.minQ30.mx1.map", quiet = FALSE,
                                     parallel = require("parallel"), n.cores = NULL)))


m6.05<-as.genind(as.matrix(read.PLINK("HESS_fix/gencurve/vcf.MMDP60.maf0.006.minQ30.mx05.raw", 
                                      map.file = "HESS_fix/gencurve/vcf.MMDP60.maf0.006.minQ30.mx05.map", quiet = FALSE,
                                      parallel = require("parallel"), n.cores = NULL)))

m6.0<-as.genind(as.matrix(read.PLINK("HESS_fix/gencurve/vcf.MMDP60.maf0.006.minQ30.mx0.raw", 
                                     map.file = "HESS_fix/gencurve/vcf.MMDP60.maf0.006.minQ30.mx0.map", quiet = FALSE,
                                     parallel = require("parallel"), n.cores = NULL)))


m49.1<-as.genind(as.matrix(read.PLINK("HESS_fix/gencurve/vcf.MMDP60.maf0.49.minQ30.mx1.raw", 
                                      map.file = "HESS_fix/gencurve/vcf.MMDP60.maf0.49.minQ30.mx1.map", quiet = FALSE,
                                      parallel = require("parallel"), n.cores = NULL)))

m49.05<-as.genind(as.matrix(read.PLINK("HESS_fix/gencurve/vcf.MMDP60.maf0.49.minQ30.mx05.raw", 
                                       map.file = "HESS_fix/gencurve/vcf.MMDP60.maf0.49.minQ30.mx05.map", quiet = FALSE,
                                       parallel = require("parallel"), n.cores = NULL)))

m49.0<-as.genind(as.matrix(read.PLINK("HESS_fix/gencurve/vcf.MMDP60.maf0.49.minQ30.mx01.raw", 
                                      map.file = "HESS_fix/gencurve/vcf.MMDP60.maf0.49.minQ30.mx1.map", quiet = FALSE,
                                      parallel = require("parallel"), n.cores = NULL)))

# saveRDS(m6.1, "HESS_fix/gencurve/robjs/m6.1.rds")
# saveRDS(m6.05, "HESS_fix/gencurve/robjs/m6.05.rds")
# saveRDS(m6.0, "HESS_fix/gencurve/robjs/m6.0.rds")
# 
# saveRDS(m49.1, "HESS_fix/gencurve/robjs/m49.1.rds")
# saveRDS(m49.05, "HESS_fix/gencurve/robjs/m49.05.rds")
# saveRDS(m49.0, "HESS_fix/gencurve/robjs/m49.0.rds")





# now import them to cluster
# run script GENETS_GC.aux.R in cluster





######################################################################
# compile

############ 250loci

setwd("~/Dropbox/finalizing_genets/")

pm6.0<-read.csv("HESS_fix/gencurve/pink_out/m6.0.rds.genoc200_compile.csv", head=T, stringsAsFactors = F)
pm6.05<-read.csv("HESS_fix/gencurve/pink_out/m6.05.rds.genoc200_compile.csv", head=T, stringsAsFactors = F)
pm6.1<-read.csv("HESS_fix/gencurve/pink_out/m6.1.rds.genoc200_compile.csv", head=T, stringsAsFactors = F)

pm49.0<-read.csv("HESS_fix/gencurve/pink_out/m49.0.rds.genoc200_compile.csv", head=T, stringsAsFactors = F)
pm49.05<-read.csv("HESS_fix/gencurve/pink_out/m49.05.rds.genoc200_compile.csv", head=T, stringsAsFactors = F)
pm49.1<-read.csv("HESS_fix/gencurve/pink_out/m49.1.rds.genoc200_compile.csv", head=T, stringsAsFactors = F)

pm6.0$lab<-"pm6.0"
pm6.05$lab<-"pm6.05"
pm6.1$lab<-"pm6.1"

pm49.0$lab<-"pm49.0"
pm49.05$lab<-"pm49.05"
pm49.1$lab<-"pm49.1"

rep(1, 250)
numbers<-c(1:250)
iterat<-rep(numbers, each=1000)

pm6.0$iterat<-iterat
pm6.05$iterat<-iterat
pm6.1$iterat<-iterat

pm49.0$iterat<-iterat
pm49.05$iterat<-iterat
pm49.1$iterat<-iterat

gencurve.compile.novo<-rbind(pm6.0,
                             pm6.05,
                             pm6.1,
                             pm49.0,
                             pm49.05,
                             pm49.1)

# data_summary <- function(data, varname, groupnames){
#   require(plyr)
#   summary_func <- function(x, col){
#     c(mean = mean(x[[col]], na.rm=TRUE),
#       sd = sd(x[[col]], na.rm=TRUE))
#   }
#   data_sum<-ddply(data, groupnames, .fun=summary_func,
#                   varname)
#   data_sum <- rename(data_sum, c("mean" = varname))
#   return(data_sum)
# }
# 
# pm_sum <- data_summary(gencurve.compile.novo, varname="MLG", 
#                        groupnames=c("lab", "iterat"))


library(dplyr)
bubbles <- gencurve.compile.novo %>%
  group_by(lab, iterat) %>%
  summarise( 
    n=n(),
    mean=mean(MLG),
    sd=sd(MLG)
  ) %>%
  mutate (se=sd/sqrt(n)) 

pm_sum<-bubbles
colnames(pm_sum)[4]<-"MLG"

# Convert iterta to a factor variable
pm_sum$iterat=as.factor(pm_sum$iterat)
head(pm_sum)

library(ggplot2)

ggplot(pm_sum, aes(x=iterat, y=MLG, group=lab, color=factor(lab))) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=MLG-sd, ymax=MLG+sd), width=.2,
                position=position_dodge(0.05))+
  theme(axis.text.x = element_text(angle=45),
        panel.background = element_rect(fill = "NA", colour = "black"), legend.text=element_text(size=7),legend.key=element_blank(),legend.position = "right",axis.text=element_text(size=6),
        axis.title=element_text(size=10))+
  labs(x="Number of SNP Loci Randomly Sampled", y="Multilocus Genotypes Identified", color="VCF Filter")+
  geom_hline(yintercept=86, linetype="solid", color = "darkgrey", size=1)+
  annotate("text", x = 10, y = 87, label = "Total Sporocarps")+
  geom_hline(yintercept=86, linetype="dotted", color = "red", size=1)+
  annotate("text", x = 10, y = 85, label = "Maximum Multilocus Genotypes Identified")+
  scale_color_manual(values=c("#ffe119", "#4363d8", "#3cb44b", "#911eb4", "#000000", "#e6194B"))

#no legend
snpgc<-ggplot(pm_sum, aes(x=iterat, y=MLG, group=lab, color=factor(lab))) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=MLG-sd, ymax=MLG+sd), width=.2,
                position=position_dodge(0.05))+
  theme(axis.text.x = element_text(angle=45),
        panel.background = element_rect(fill = "NA", colour = "black"), legend.text=element_text(size=7),legend.key=element_blank(),legend.position = "none",axis.text=element_text(size=6),
        axis.title=element_text(size=10))+
  labs( x="Number of SNP Loci Randomly Sampled", y="Multilocus Genotypes Identified", color="VCF Filter")+
  geom_hline(yintercept=86, linetype="solid", color = "darkgrey", size=1)+
  annotate("text", x = 15, y = 85, label = "Total Sporocarps")+
  geom_hline(yintercept=86, linetype="dotted", color = "red", size=1)+
  annotate("text", x = 15, y = 87, label = "Maximum Multilocus Genotypes Identified")+
  scale_color_manual(values=c("#ffe119", "#4363d8", "#3cb44b", "#911eb4", "#000000", "#e6194B"))

# ############ 100loci
# setwd("~/Dropbox/finalizing_genets/")
# 
# pm6.0<-read.csv("HESS_fix/gencurve/pink_out/m6.0.rds.genoc_compile.csv", head=T, stringsAsFactors = F)
# pm6.05<-read.csv("HESS_fix/gencurve/pink_out/m6.05.rds.genoc_compile.csv", head=T, stringsAsFactors = F)
# pm6.1<-read.csv("HESS_fix/gencurve/pink_out/m6.1.rds.genoc_compile.csv", head=T, stringsAsFactors = F)
# 
# pm49.0<-read.csv("HESS_fix/gencurve/pink_out/m49.0.rds.genoc_compile.csv", head=T, stringsAsFactors = F)
# pm49.05<-read.csv("HESS_fix/gencurve/pink_out/m49.05.rds.genoc_compile.csv", head=T, stringsAsFactors = F)
# pm49.1<-read.csv("HESS_fix/gencurve/pink_out/m49.1.rds.genoc_compile.csv", head=T, stringsAsFactors = F)
# 
# pm6.0$lab<-"pm6.0"
# pm6.05$lab<-"pm6.05"
# pm6.1$lab<-"pm6.1"
# 
# pm49.0$lab<-"pm49.0"
# pm49.05$lab<-"pm49.05"
# pm49.1$lab<-"pm49.1"
# 
# rep(1, 100)
# numbers<-c(1:100)
# iterat<-rep(numbers, each=1000)
# 
# pm6.0$iterat<-iterat
# pm6.05$iterat<-iterat
# pm6.1$iterat<-iterat
# 
# pm49.0$iterat<-iterat
# pm49.05$iterat<-iterat
# pm49.1$iterat<-iterat
# 
# gencurve.compile.novo<-rbind(pm6.0,
#                              pm6.05,
#                              pm6.1,
#                              pm49.0,
#                              pm49.05,
#                              pm49.1)
# 
# data_summary <- function(data, varname, groupnames){
#   require(plyr)
#   summary_func <- function(x, col){
#     c(mean = mean(x[[col]], na.rm=TRUE),
#       sd = sd(x[[col]], na.rm=TRUE))
#   }
#   data_sum<-ddply(data, groupnames, .fun=summary_func,
#                   varname)
#   data_sum <- rename(data_sum, c("mean" = varname))
#   return(data_sum)
# }
# 
# pm_sum <- data_summary(gencurve.compile.novo, varname="MLG", 
#                     groupnames=c("lab", "iterat"))
# 
# # Convert iterta to a factor variable
# pm_sum$iterat=as.factor(pm_sum$iterat)
# head(pm_sum)
# 
# library(ggplot2)
# 
# ggplot(pm_sum, aes(x=iterat, y=MLG, group=lab, color=factor(lab))) + 
#   geom_line() +
#   geom_point()+
#   geom_errorbar(aes(ymin=MLG-sd, ymax=MLG+sd), width=.2,
#                 position=position_dodge(0.05))+
#   theme(axis.text.x = element_text(angle=45),
#         panel.background = element_rect(fill = "NA", colour = "black"), legend.text=element_text(size=7),legend.key=element_blank(),legend.position = "right",axis.text=element_text(size=6),
#         axis.title=element_text(size=10))+
#   labs(x="Number of SNP Loci Randomly Sampled", y="Multilocus Genotypes Identified", color="VCF Filter")+
#   geom_hline(yintercept=86, linetype="solid", color = "darkgrey", size=1)+
#   annotate("text", x = 10, y = 87, label = "Total Sporocarps")+
#   geom_hline(yintercept=86, linetype="dotted", color = "red", size=1)+
#   annotate("text", x = 10, y = 85, label = "Maximum Multilocus Genotypes Identified")+
#   scale_color_manual(values=c("#ffe119", "#4363d8", "#3cb44b", "#911eb4", "#000000", "#e6194B"))
# 
# #no legend
# ggplot(pm_sum, aes(x=iterat, y=MLG, group=lab, color=factor(lab))) + 
#   geom_line() +
#   geom_point()+
#   geom_errorbar(aes(ymin=MLG-sd, ymax=MLG+sd), width=.2,
#                 position=position_dodge(0.05))+
#   theme(axis.text.x = element_text(angle=45),
#         panel.background = element_rect(fill = "NA", colour = "black"), legend.text=element_text(size=7),legend.key=element_blank(),legend.position = "none",axis.text=element_text(size=6),
#         axis.title=element_text(size=10))+
#   labs( x="Number of SNP Loci Randomly Sampled", y="Multilocus Genotypes Identified", color="VCF Filter")+
#   geom_hline(yintercept=86, linetype="solid", color = "darkgrey", size=1)+
#   annotate("text", x = 15, y = 85, label = "Total Sporocarps")+
#   geom_hline(yintercept=86, linetype="dotted", color = "red", size=1)+
#   annotate("text", x = 15, y = 87, label = "Maximum Multilocus Genotypes Identified")+
#   scale_color_manual(values=c("#ffe119", "#4363d8", "#3cb44b", "#911eb4", "#000000", "#e6194B"))



# # # # # # # # # # # # # # # # # #
#AFLP
setwd("~/Dropbox/finalizing_genets/gencurves.raw")

hiergc.102<-read.csv("hier.gi.rds.y100.102.genoc_compile.csv", head=T, stringsAsFactors = F)
hiergc<-read.csv("hier.gi.rds.genoc_compile.csv", head=T, stringsAsFactors = F)

rep(1, 100)
numbers<-c(1:100)
iterat<-rep(numbers, each=1000)
iterat2<-rep(numbers, each=102)

hiergc$iterat<-iterat
hiergc.102$iterat<-iterat2

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}


df.hier<- data_summary(hiergc, varname="MLG", 
                       groupnames=c("lab", "iterat"))

df.hier.102<- data_summary(hiergc.102, varname="MLG", 
                           groupnames=c("lab", "iterat"))


# Convert iterta to a factor variable
df2$iterat=as.factor(df2$iterat)
head(df2)

df.hier$iterat=as.factor(df.hier$iterat)
df.hier.102$iterat=as.factor(df.hier.102$iterat)

library(ggplot2)

aflpgc<-ggplot(df.hier, aes(x=iterat, y=MLG, group=lab)) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=MLG-sd, ymax=MLG+sd), width=.2,
                position=position_dodge(0.05))+
  theme(axis.text.x = element_text(angle=45),
        panel.background = element_rect(fill = "NA", colour = "black"), legend.text=element_text(size=7),legend.key=element_blank(),legend.position = "bottom",axis.text=element_text(size=6),
        axis.title=element_text(size=10))+
  labs(x="Number of AFLP Loci Randomly Sampled", y="Multilocus Genotypes Identified", color="VCF Filter")+
  geom_hline(yintercept=221, linetype="solid", color = "darkgrey", size=1)+
  annotate("text", x = 10, y = 218, label = "Total Sporocarps")+
  geom_hline(yintercept=160, linetype="dotted", color = "red", size=1)+
  annotate("text", x = 10, y = 163, label = "Maximum Multilocus Genotypes Identified")

#  panel
library(gridExtra)
r <- rectGrob(gp=gpar(fill="white"))
grid.arrange(aflpgc, r, snpgc,r, nrow=2, ncol=2)










