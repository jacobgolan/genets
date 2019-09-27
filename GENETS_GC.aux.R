library(poppr)

GENINPUT<-readRDS("----.rds", refhook = NULL)
file<-"---.rds"

GENINPUT.gi<-GENINPUT


genoc<-vector("list")
y<-1
#how many loci to consider = 250
while(y<250+1){
  
  #create 1000 random integers w/in range of no. loci
  sub.list<-vector("list")
  x<-1
  #randomly sample 1-100 loci, 1000 times per y
  while(x<1000+1){
    sub.list[[x]]<-sample.int(ncol(GENINPUT.gi$tab), y) #try different no. (=range of y) of loci, loop 100 times
    x<-x+1
  }
  
  #create 1000 randomly subsampled genind objects
  mafsub.list<-vector("list")
  x<-1
  while(x<length(sub.list)+1){
    mafsub.list[x]<-GENINPUT.gi[loc=sub.list[[x]]]
    x<-x+1
  }
  
  pop.list<-vector("list")
  x<-1
  while(x<length(mafsub.list)+1){
    pop.list[[x]]<-as.data.frame(poppr(mafsub.list[[x]]))
    x<-x+1
  }
  
  #merge into single df
  library(plyr)
  genoc.pre <- ldply(pop.list, data.frame)
  genoc.pre$lab<-paste("it_",y, sep='')
  genoc[[y]]<-genoc.pre
  
  y<-y+1
}

genoc.compile <- ldply(genoc, data.frame)
write.csv(genoc.compile, paste(file,".genoc_compile.csv",sep=""))

          