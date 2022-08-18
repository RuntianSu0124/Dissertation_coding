#produce plots
set.seed(123)
yeast_geno<-matrix(NA,nrow=1008,ncol=11623)
int1<-1
int2<-0
for(i in 1:16){
  int2<-int2+length(cross$geno[[i]]$data)/1008
  yeast_geno[,int1:int2]<-as.matrix(cross$geno[[i]]$data)
  int1<-int1+length(cross$geno[[i]]$data)/1008
}

i=2 #use the dataset under Caffeine growth condition

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

##Ridge regression
cv_fit<-cv.glmnet(X.train,Y.train,alpha=0,lambda=lambdas,nfolds=10)
lambda<-cv_fit$lambda.min
RR_model<-glmnet(X.train,Y.train,alpha=0, lambda=lambda)
predict_y<-predict(RR_model,s=lambda,newx=X.test)
rr.model<-glmnet(X.train,Y.train,alpha=0,lambda=lambdas)
plot(rr.model,xvar='lambda')
plot(cv_fit)
plot(x=Y.test,y=predict_y,xlab="Actual",ylab="Prediction",
     main="Actual VS Prediction(Ridge regression)")
trainlinmod<-lm(predict_y~Y.test)
abline(trainlinmod,col=3,lwd=2.5,lty="solid")
abline(a=0,b=1,col=2,lwd=2.5,lty="dashed")
legend("topleft",
       legen=c("Model","Base"),
       col=c("green","red"),text.width = 0.5,bty="n",
       lwd=2.5,lty=c("solid","dashed"))

##Lasso regression
cv_fit<-cv.glmnet(X.train,Y.train,alpha=1,lambda=lambdas,nfolds=10)#Cross-validation
lambda<-cv_fit$lambda.min #the optimal lambda
LA_model<-glmnet(X.train,Y.train,alpha=1, lambda=lambda)
predict_y<-predict(LA_model,s=lambda,newx=X.test)
LA_models<-glmnet(X.train,Y.train,alpha=1, lambda=lambdas)
plot(LA_models,xvar='lambda')
plot(cv_fit)
plot(x=Y.test,y=predict_y,xlab="Actual",ylab="Prediction",
     main="Actual VS Prediction(Lasso regression)")
trainlinmod<-lm(predict_y~Y.test)
abline(trainlinmod,col=3,lwd=2.5,lty="solid")
abline(a=0,b=1,col=2,lwd=2.5,lty="dashed")
legend("topleft",
       legen=c("Model","Base"),
       col=c("green","red"),text.width = 0.5,bty="n",
       lwd=2.5,lty=c("solid","dashed"))

##Elastic net



##OLS


##PCR


##Ramdom Forest



##GBLUP