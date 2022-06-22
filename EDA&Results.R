## EDA plots


## PCA : Dimension reduction and Visualization
## PCA based on covariance matrix
geno_pca<-prcomp(t(geno_matrix),scale=FALSE)
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
## I’ve drawn on horizontal lines at 90% and 95% of variance explained, to help identify when we
## cross these thresholds. We need 13 components to explain 90% of the variance,
## and 25 components to explain 95% of the variance.




## Ridge regression
## use the elbow method to find the optimal lambda
## yeast growing condition: Caffeine
library(glmnet)
Caffeine_y<-c(cross$pheno$Caffeine)
nn<-c(which(Caffeine_y>-1000))## Exclusion of missing values

lambdas<-10^seq(3,-4,by=-.1)
life.ridge<-glmnet(yeast_geno[nn,],Caffeine_y[nn],alpha=0,lambda=lambdas)
plot(life.ridge,xvar='lambda')
cv_fit<-cv.glmnet(yeast_geno[nn,],Caffeine_y[nn],alpha=0,lambda=lambdas)
plot(cv_fit) 
cv_fit$lambda.min ## the optimal lambda
## fit a ridge regression model
ridge_m1<-glmnet(yeast_geno[nn,],Caffeine_y[nn],alpha=0, lambda=cv_fit$lambda.min)
## Calculate r-squared
y_bar<-mean(Caffeine_y[nn])
predict_y<-predict(ridge_m1,s=cv_fit$lambda.min,newx=yeast_geno[nn,])
R_Square =  1 - sum((predict_y-y_bar)^2)/sum((Caffeine_y[nn]-y_bar)^2)

## Cross-validation(Caffeine condition)
Caffeine_y_num<-Caffeine_y[nn]
yeast_geno_num<-yeast_geno[nn,]
num<-length(Caffeine_y_num)
R_Squares<-c(NA,length=10)
for(i in 1:10){
  ##Dividing the training and test sets
  train_num<-sample(num,0.8*num)
  yeast_geno_train<-yeast_geno_num[train_num,]
  yeast_geno_test<-yeast_geno_num[-train_num,]
  Caffeine_y_train<-Caffeine_y_num[train_num]
  Caffeine_y_test<-Caffeine_y_num[-train_num]
  
  ##ridge regression
  cv_fit<-cv.glmnet(yeast_geno_train,Caffeine_y_train,alpha=0,lambda=lambdas)
  ridge_model<-glmnet(yeast_geno_train,Caffeine_y_train,alpha=0,lambda=cv_fit$lambda.min)
  y_bar<-mean(Caffeine_y_train)
  predict_y<-predict(ridge_model,s=cv_fit$lambda.min,newx=yeast_geno_test)
  R_Squares[i]<-1-sum((predict_y-y_bar)^2)/sum((Caffeine_y_test-y_bar)^2)
}
R_Squares
CV_Rsqare<-mean(R_Squares)
CV_Rsqare










