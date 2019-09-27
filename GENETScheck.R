confui<-read.PLINK("HESS_fix/mx0.5/vcf.MMDP60.maf0.006.minQ30.recode.raw", 
                   map.file = "HESS_fix/mx0.5/vcf.MMDP60.maf0.006.minQ30.recode.map", quiet = FALSE,
                   parallel = require("parallel"), n.cores = NULL)
dat <- list(toto=c(1,1,0,0), titi=c(NA,1,1,0), tata=c(NA,0,3, NA))
x <- new("genlight", dat)
x


library(vcfR)
vcf <- read.vcfR( "HESS_fix/pops/vcf.MMDP60.maf0.006.minQ30.recode.D2.14.HWE.recode.vcf", verbose = FALSE )
x <- vcfR2genlight(vcf)

gt <- extract.gt(vcf, element = "GT")
gt[c(2,6,18), 1:3]







my_genind <- vcfR2genind(vcf)


vcf2genind <- function(filename,otherUsefullOption){
  my.loci <- read.vcfR(file)
  options(adegenet.check.ploidy = FALSE)
  my.genind <- vcfR2genind(my.loci)
  options(adegenet.check.ploidy = TRUE)
}


any(as.matrix(confui)==2)

as.data.frame(table(as.matrix(confui)))
