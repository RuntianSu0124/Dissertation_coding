install.packages("randomForest")
library("randomForest")
library(pls)
library(glmnet)

## EDA plots


## Ridge regression
# use the elbow method to find the optimal lambda
# yeast growing condition: Caffeine
library(glmnet)
Caffeine_y<-c(cross$pheno$Caffeine)
nn<-c(which(Caffeine_y>-1000))## Exclusion of missing values

lambdas<-10^seq(3,-4,by=-.1)
life.ridge<-glmnet(yeast_geno[nn,],Caffeine_y[nn],alpha=0,lambda=lambdas)
plot(life.ridge,xvar='lambda')
cv_fit<-cv.glmnet(yeast_geno[nn,],Caffeine_y[nn],alpha=0,lambda=lambdas)
plot(cv_fit) 
cv_fit$lambda.min # the optimal lambda
# fit a ridge regression model
ridge_m1<-glmnet(yeast_geno[nn,],Caffeine_y[nn],alpha=0, lambda=cv_fit$lambda.min)
# Calculate r-squared
y_bar<-mean(Caffeine_y[nn])
predict_y<-predict(ridge_m1,s=cv_fit$lambda.min,newx=yeast_geno[nn,])
R_Square =  1 - sum((predict_y-y_bar)^2)/sum((Caffeine_y[nn]-y_bar)^2)

# Ridge regression Cross-validation(Caffeine condition)
Caffeine_y_num<-Caffeine_y[nn]
yeast_geno_num<-yeast_geno[nn,]
num<-length(Caffeine_y_num)
R_Squares.RR<-c(NA,length=10)
for(i in 1:10){
  #Dividing the training and test set
  train_num<-sample(num,0.8*num)
  yeast_geno_train<-yeast_geno_num[train_num,]
  yeast_geno_test<-yeast_geno_num[-train_num,]
  Caffeine_y_train<-Caffeine_y_num[train_num]
  Caffeine_y_test<-Caffeine_y_num[-train_num]
  
  #ridge regression
  cv_fit<-cv.glmnet(yeast_geno_train,Caffeine_y_train,alpha=0,lambda=lambdas)
  ridge_model<-glmnet(yeast_geno_train,Caffeine_y_train,alpha=0,lambda=cv_fit$lambda.min)
  y_bar<-mean(Caffeine_y_train)
  predict_y<-predict(ridge_model,s=cv_fit$lambda.min,newx=yeast_geno_test)
  R_Squares.RR[i]<-1-sum((predict_y-y_bar)^2)/sum((Caffeine_y_test-y_bar)^2)
}
R_Squares.RR
CV_Rsqare.RR<-mean(R_Squares.RR)
CV_Rsqare.RR

# coefficient of Covariate >0 or <0?
which(ridge_model$beta<0)
which(ridge_model$beta>0)


## PCA : Dimension reduction and Visualization
## PCA based on covariance matrix
geno_pca<-prcomp(yeast_geno,scale=FALSE)
geno_pca
summary(geno_pca)
plot(geno_pca$sdev)#scree plot
cumvar<-100*cumsum(geno_pca$sdev^2)/sum(geno_pca$sdev^2)
plot(cumvar,ylab="Cumulative proportion of variance explained",
     xlab="Number of PCs used", ylim=c(0,100))
abline(h=90, lty=2)
abline(v=min(which(cumvar>90)), lty=2)
abline(h=95, lty=2)
abline(v=min(which(cumvar>95)), lty=2)
cumvar
# I’ve drawn on horizontal lines at 90% and 95% of variance explained, to help identify when we
# cross these thresholds. We need 146 components to explain 90% of the variance,
# and 237 components to explain 95% of the variance.

## PCR method, start from scratch
geno_146PC<-geno_pca$x[nn,1:146]#select the first 146 PCs
geno_237PC<-geno_pca$x[nn,1:237]#select the first 237 PCs
pheno_PCR1<-lm(Caffeine_y[nn]~geno_146PC)
pheno_PCR2<-lm(Caffeine_y[nn]~geno_237PC)
plot(Caffeine_y,type="l",col=1)
lines(Caffeine_y[nn])
# Calculate R-square of PCR model
summary(pheno_PCR2)

# PCR Cross-validation (under Caffeine condition)
# Start from scratch
R_Squares.PCR<-c(NA,10)
for(i in 1:10){
  #Dividing the training and test set
  train_num<-sample(num,0.8*num)
  yeast_geno_train<-yeast_geno_num[train_num,]
  yeast_geno_test<-yeast_geno_num[-train_num,]
  Caffeine_y_train<-Caffeine_y_num[train_num]
  Caffeine_y_test<-Caffeine_y_num[-train_num]
  
  #PCR
  geno_pca<-prcomp(yeast_geno_num,scale=FALSE)
  cumvar<-100*cumsum(geno_pca$sdev^2)/sum(geno_pca$sdev^2)
  npc<-which(cumvar>99)[1] #the amount of selected PCs
  PCS<-geno_pca$x[train_num,1:npc]
  pheno_PCR<-lm(Caffeine_y_train~PCS)
  PC.test<-geno_pca$x[-train_num,1:npc]
  #Calculate R-square of PCR model
  y_bar<-mean(Caffeine_y_train)
  predict_y<-predict.lm(pheno_PCR,newdata=data.frame(PC.test))
  R_Squares.PCR[i]<-1-sum((predict_y-y_bar)^2)/sum((Caffeine_y_test-y_bar)^2)
}
R_Squares.PCR
CV_Rsqare.RCR<-mean(R_Squares.PCR)
CV_Rsqare.RCR

## PCR Cross-validation:use prc function
R_Squares.PCR<-c(NA,10)
for(i in 1:10){
  #Dividing the training and test set
  train_num<-sample(num,0.8*num)
  yeast_geno_train<-yeast_geno_num[train_num,]
  yeast_geno_test<-yeast_geno_num[-train_num,]
  Caffeine_y_train<-Caffeine_y_num[train_num]
  Caffeine_y_test<-Caffeine_y_num[-train_num]
  
  #PCR
  pheno_PCR<-pcr(Caffeine_y_train~yeast_geno_train,scale=FALSE,validation="CV")
  predict_y<-predict(pheno_PCR,yeast_geno_test,ncomp=400)
  
  #Calculate R-square of PCR model
  y_bar<-mean(Caffeine_y_train)
  R_Squares.PCR[i]<-1-sum((predict_y-y_bar)^2)/sum((Caffeine_y_test-y_bar)^2)
}
R_Squares.PCR
CV_Rsqare.RCR<-mean(R_Squares.PCR)
CV_Rsqare.RCR
pheno_PCR$scores

## PCR Cross-validation:use prc function 2
  #Dividing the training and test set
  train_num<-sample(num,0.8*num)
  yeast_geno_train<-yeast_geno_num[train_num,]
  yeast_geno_test<-yeast_geno_num[-train_num,]
  Caffeine_y_train<-Caffeine_y_num[train_num]
  Caffeine_y_test<-Caffeine_y_num[-train_num]
  
  #PCR
  pheno_PCR<-pcr(Caffeine_y_train~yeast_geno_train,scale=FALSE,validation="CV")
  
  predict_y<-predict(pheno_PCR,yeast_geno_test,ncomp=500)
  
  #Calculate R-square of PCR model
  y_bar<-mean(Caffeine_y_train)
  R_Squares.PCR2<-1-sum((predict_y-y_bar)^2)/sum((Caffeine_y_test-y_bar)^2)
  sqrt(mean((predict_y - Caffeine_y_test)^2))
  
  
R_Squares.PCR2
CV_Rsqare.RCR<-mean(R_Squares.PCR)
CV_Rsqare.RCR
pheno_PCR$scores


## Random Forest regression 

