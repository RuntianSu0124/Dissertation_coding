library(pls)
library(glmnet)
##store the yeast genotype data in a matrix
yeast_geno<-matrix(NA,nrow=1008,ncol=11623)
int1<-1
int2<-0
for(i in 1:16){
  int2<-int2+length(cross$geno[[i]]$data)/1008
  yeast_geno[,int1:int2]<-as.matrix(cross$geno[[i]]$data)
  int1<-int1+length(cross$geno[[i]]$data)/1008
}


##Ridge regression:Alpha=0
lambdas<-10^seq(3,-4,by=-.1)
missings<-c(NA,length=46)
R2_ridge<-c(NA,length=46)
for(i in 1:46){
Y_raw<-c(cross$pheno[[i]])
nn<-c(which(Y_raw>-1000))## Exclusion of missing values
missings[i]<-1008-length(nn)
Y<-Y_raw[nn] #pheotype data without missing values
X<-yeast_geno[nn,]
#Dividing the training and test set
num<-sample(length(Y),0.8*length(Y))
Y.train<-Y[num]
X.train<-X[num,]
Y.test<-Y[-num]
X.test<-X[-num,]
#train the model by using train set 
cv_fit<-cv.glmnet(X.train,Y.train,alpha=0,lambda=lambdas,nfolds=10)#10 fold Cross-validation
lambda<-cv_fit$lambda.min #the optimal lambda
RR_model<-glmnet(X.train,Y.train,alpha=0, lambda=lambda)
#calculate R-square
y_bar<-mean(Y.test)
predict_y<-predict(RR_model,s=lambda,newx=X.test)
R2_ridge[i]<-1-sum((predict_y-y_bar)^2)/sum((Y.test-y_bar)^2)
}
R2_ridge


##Lasso regression:Alpha=1
R2_lasso<-c(NA,length=46)
for(i in 1:46){
  Y_raw<-c(cross$pheno[[i]])
  nn<-c(which(Y_raw>-1000))## Exclusion of missing values
  Y<-Y_raw[nn] #pheotype data without missing values
  X<-yeast_geno[nn,]
  #Dividing the training and test set
  num<-sample(length(Y),0.8*length(Y))
  Y.train<-Y[num]
  X.train<-X[num,]
  Y.test<-Y[-num]
  X.test<-X[-num,]
  #train the model by using train set 
  cv_fit<-cv.glmnet(X.train,Y.train,alpha=1,lambda=lambdas,nfolds=10)#10 fold Cross-validation
  lambda<-cv_fit$lambda.min #the optimal lambda
  LA_model<-glmnet(X.train,Y.train,alpha=1, lambda=lambda)
  #calculate R-square on test set
  y_bar<-mean(Y.test)
  predict_y<-predict(LA_model,s=lambda,newx=X.test)
  R2_lasso[i]<-1-sum((predict_y-y_bar)^2)/sum((Y.test-y_bar)^2)
}
R2_lasso


##Elastic net: Alpha between 0 to 1
Alphas<-seq(0,1,by=0.05)
R2_elastic<-c(NA,length=46)
alphas<-c(NA,length=46)
for(i in 1:46){
  R2_Alphas<-c(NA,length=21)
  for(k in 1:21){
  Y_raw<-c(cross$pheno[[i]])
  nn<-c(which(Y_raw>-1000))## Exclusion of missing values
  Y<-Y_raw[nn] #pheotype data without missing values
  X<-yeast_geno[nn,]
  #Dividing the training and test set
  num<-sample(length(Y),0.8*length(Y))
  Y.train<-Y[num]
  X.train<-X[num,]
  Y.test<-Y[-num]
  X.test<-X[-num,]
  #
  cv_fit<-cv.glmnet(X.train,Y.train,alpha=Alphas[k],lambda=lambdas,nfolds=10)#10 fold Cross-validation
  lambda<-cv_fit$lambda.min #the optimal lambda
  EN_model<-glmnet(X,Y,alpha=Alphas[k],lambda=lambda)
  #calculate R-square
  y_bar<-mean(Y.test)
  predict_y<-predict(EN_model,s=lambda,newx=X.test)
  R2_Alphas[k]<-1-sum((predict_y-y_bar)^2)/sum((Y.test-y_bar)^2)
  }
  alphas[i]<-which.max(R2_Alphas)
  R2_elastic[i]<-max(R2_Alphas)
}
Alphas[alphas]# the value of alpha
R2_elastic


## Principal component regression:PCR

R2_PCR<-c(NA,length=46)
for(i in 1:46){
  
  
  
}
PCR_model<-pcr(Y.train~X.train,scale=FALSE,validation="CV")



