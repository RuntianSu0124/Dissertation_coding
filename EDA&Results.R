## EDA plots





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
