## dissertation code work--load the data ##

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

cross$geno$`1`$data

##store the yeast genotype data in a matrix
yeast_geno<-matrix(NA,nrow=1008,ncol=11623)
int1<-1
int2<-0
for(i in 1:16){
  int2<-int2+length(cross$geno[[i]]$data)/1008
  yeast_geno[,int1:int2]<-as.matrix(cross$geno[[i]]$data)
  int1<-int1+length(cross$geno[[i]]$data)/1008
}




