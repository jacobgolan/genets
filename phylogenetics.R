#AFLP
setwd("~/Documents/2020_AphalloidesPopgen")
library(poppr)
library(ggtree)
library(ape)
library(ggplot2)
library(phangorn)

hier<-read.genalex("data/aflp/AFLP_genalex_AMANITABASED.csv")
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


cols<-c(Drake_1="#b2182b",
        Drake_2="#d6604d",
        Drake_3="#f4a582",
        Drake_4="#fddbc7",
        Hearts_Desire_1="#2166ac",
        Hearts_Desire_2="#67a9cf",
        Hearts_Desire_3="#d1e5f0",
        Monterrey="#888888",
        Rochester_NY_2007_1="#4d9221",
        Rochester_NY_2007_2="#a1d76a",
        Rochester_NY_2007_3="#e6f5d0",
        Jakes="#542788",
        New_Jersey="#d8daeb",
        Austria="#EF3340",
        CESAC="#002654",
        Arfons="#002654",
        Serbia="#EDB92E")


hier.gi<-as.genind(as.matrix(hier))

pop(hier.gi)<-ifelse(grepl("Rochester",pops.aflp),pops.aflp,
  ifelse(grepl("Drake",pops.aflp)|grepl("Hearts",pops.aflp),
                     sub("_[^_]*$","",pops.aflp),ifelse(pops.aflp=="New_Jersey_str_2006","New_Jersey"
                                      ,sub("_.*$","",pops.aflp))))


hier.gi<-clonecorrect(hier.gi)


tre<-upgma(diss.dist(hier.gi))

pop.col.df<-data.frame(
  ind=indNames(hier.gi),
  pop=pop(hier.gi)
)


tre_boot<-aboot(hier.gi)


p <- ggtree(tre_boot) %<+% pop.col.df
dir.create("results/phylogeny")
pdf("results/phylogeny/clonecorrectUPGMA.pdf",width=8,height=10)
p + geom_tiplab(offset = 0, hjust = 0,size=2) +
  geom_tippoint(aes(color = pop), size=1) + 
  theme(legend.position = "right")+
  scale_color_manual(values=cols)+
  geom_text(aes(label = c(NA,rep(NA,nrow(tre_boot$edge)-length(tre_boot$node.label)),
                        tre_boot$node.label)), hjust = 1, vjust = -0.4, size = 3)

dev.off()

pdf("results/phylogeny/clonecorrectUPGMA-notip.pdf",width=8,height=10)
p + geom_tree() + theme_tree() +
  geom_tippoint(aes(color = pop), size=1.5) + 
  theme(legend.position = "right")+
  scale_color_manual(values=cols)#+
  #geom_text(aes(label = c(NA,rep(NA,nrow(tre_boot$edge)-length(tre_boot$node.label)),
  #                        tre_boot$node.label)), hjust = 1, vjust = -0.4, size = 3)

dev.off()

### CA + CESAC

KEEP<-c("Aph275",
        indNames(popsub(hier.gi, sublist = c("Drake_1","Drake_2","Drake_3","Drake_4","Hearts_Desire_1","Hearts_Desire_2","Hearts_Desire_3","Monterrey"))))
CA.plus<-hier.gi[row.names(hier.gi@tab) %in% KEEP]
tre<-upgma(diss.dist(CA.plus))

pop.col.df<-data.frame(
  ind=indNames(CA.plus),
  pop=pop(CA.plus)
)




tre_boot<-aboot(CA.plus)

p <- ggtree(tre_boot) %<+% pop.col.df
pdf("results/phylogeny/CASerbiaUPGMA.pdf",width=8,height=10)
p + geom_tiplab(offset = 0, hjust = 0,size=2) +
  geom_tippoint(aes(color = pop), size=1) + 
  theme(legend.position = "right")+
  scale_color_manual(values=cols[names(cols)%in%pop(CA.plus)])+
  geom_text(aes(label = c(NA,rep(NA,nrow(tre_boot$edge)-length(tre_boot$node.label)),
                          tre_boot$node.label)), hjust = 1, vjust = -0.4, size = 3)
dev.off()

p <- ggtree(tre_boot) %<+% pop.col.df
pdf("results/phylogeny/CASerbiaUPGMA-notip.pdf",width=8,height=10)
p + geom_tree() + theme_tree() +
  geom_tippoint(aes(color = pop), size=5) + 
  theme(legend.position = "right")+
  scale_color_manual(values=cols[names(cols)%in%pop(CA.plus)])#+
  #geom_text(aes(label = c(NA,rep(NA,nrow(tre_boot$edge)-length(tre_boot$node.label)),
  #                        tre_boot$node.label)), hjust = 1, vjust = -0.4, size = 3)
dev.off()


pdf("results/phylogeny/CASerbianbNet.pdf",width=8,height=10)
nbNet<-neighborNet(diss.dist(CA.plus))
nbnetcols<-cols[as.character(pop.col.df$pop[match(nbNet$tip.label,pop.col.df$ind)])]
plot(nbNet,tip.color = nbnetcols)
dev.off()


pdf("results/phylogeny/CASerbianbNet-notip.pdf",width=8,height=10)
nbNet$tip.label<-rep("O",length(nbNet$tip.label))
plot(nbNet,tip.color = nbnetcols)
dev.off()



### East Coast+Serbia

KEEP<-c("Aph275",
        indNames(popsub(hier.gi, sublist = c("Rochester_NY_2007_1","Rochester_NY_2007_2","Rochester_NY_2007_3","Jakes","New_Jersey"))))
EC.plus<-hier.gi[row.names(hier.gi@tab) %in% KEEP]
tre<-upgma(diss.dist(EC.plus))

pop.col.df<-data.frame(
  ind=indNames(EC.plus),
  pop=pop(EC.plus)
)
tre<-ape::rotate(tre,getMRCA(tre,c("str1","10136")))



tre_boot<-aboot(EC.plus)
tre_boot<-phytools::rotateNodes(tre_boot,getMRCA(tre_boot,c("str1","10136")))


p <- ggtree(tre_boot) %<+% pop.col.df
pdf("results/phylogeny/ECSerbiaUPGMA.pdf",width=8,height=10)
p + geom_tiplab(offset = 0, hjust = 0,size=2) +
  geom_tippoint(aes(color = pop), size=1) + 
  theme(legend.position = "right")+
  scale_color_manual(values=cols[names(cols)%in%pop(EC.plus)])+
  geom_text(aes(label = c(NA,rep(NA,nrow(tre_boot$edge)-length(tre_boot$node.label)),
                          tre_boot$node.label)), hjust = 1, vjust = -0.4, size = 3)
dev.off()

p <- ggtree(tre_boot) %<+% pop.col.df
p + geom_tree() + theme_tree() +
  geom_tippoint(aes(color = pop), size=5) + 
  theme(legend.position = "right")+
  scale_color_manual(values=cols[names(cols)%in%pop(EC.plus)])

rotate(node=getMRCA(tre_boot,c("str1","10136")))
rotate(node=getMRCA(tre_boot,c("str12","10057")))
pdf("results/phylogeny/ECSerbiaUPGMA-notip.pdf",width=8,height=10)
rotate(node=getMRCA(tre_boot,c("10134","10128")))
dev.off()




 