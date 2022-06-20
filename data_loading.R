## dissertation code work ##

##load the data
setwd("e:/rdata")
load("pheno_raw.Rdata") 
load("cross.Rdata")

length(pheno_raw$Cisplatin)
cross$pheno$Cadmium_Chloride

length(cross$geno$`1`$data)/1008
head(cross$geno$`1`$data)
cross$geno$`1`$data

geno_matrix<-as.matrix(cross$geno$`16`$data)
ncol(geno_matrix)
geno_matrix[,208]

for(i in 1:16){
a<-length(cross$geno[[i]]$data)+a
}
a/1008
cross$pheno$Caffeine

length(pheno_raw$Cadmium_Chloride)

yeast_geno<-matrix(NA,nrow=11623,ncol=1008)
num1<-c(seq(1:16))
yeast_geno<-rbind(as.matrix(cross$geno[[num1]]$data))
int1<-1
for(i in 1:16){
  num2<-length(cross$geno[[i]]$data)
  <-length(cross$geno[[i]]$data)+a
}

