#produce plots

##genotype plots
plot(yeast_geno[2,],pch=19,ann = F, xaxt = "n", yaxt = "n",xlab="Position",ylab="Allele")
axis(2,seq(1, 2, 1), seq(1, 2, 1))
axis(1,seq(1, 11623, 2000), seq(1, 11623, 2000))
for(i in 1:20){
  png(filename = paste0(i, "_", ".jpg"),width = 2400,height = 1800,res = 200)
  plot(yeast_geno[i,],pch=19,ann = F, xaxt = "n", yaxt = "n",xlab="Position",ylab="Allele")
  axis(2,seq(1, 2, 1), seq(1, 2, 1))
  axis(1,seq(1, 11623, 2000), seq(1, 11623, 2000))
  dev.off()
  }
#proportion plot
genosum<-colSums(yeast_geno)
proportion<-(genosum-1008)/1008
plot(proportion,type="l",xlab="position",ylab="proportion")

#correlation
cor<-cor(yeast_geno)
kappa(cor,exact=TRUE)

##phenotype plots
i=2 
Y_raw<-c(cross$pheno[[i]])
plot(Y_raw,ylab="Population growth",xlab="",main="Genotype data under Caffeine
condition")

for(i in 1:15){
  Y_raw<-c(cross$pheno[[i]])
  png(filename = paste0(i, "phenotype", ".jpg"),width = 2400,height = 1800,res = 200)
  plot(Y_raw,ylab="Population growth",xlab="")
  dev.off()
}




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
EL_model<-glmnet(X.train,Y.train,alpha=0.05, lambda=0.794)
predict_y<-predict(LA_model,newx=X.test)
plot(x=Y.test,y=predict_y,xlab="Actual",ylab="Prediction",
     main="Actual VS Prediction(Elastic net)")
trainlinmod<-lm(predict_y~Y.test)
abline(trainlinmod,col=3,lwd=2.5,lty="solid")
abline(a=0,b=1,col=2,lwd=2.5,lty="dashed")
legend("topleft",
       legen=c("Model","Base"),
       col=c("green","red"),text.width = 0.5,bty="n",
       lwd=2.5,lty=c("solid","dashed"))



##OLS
OLS_model<-glmnet(X.train,Y.train,alpha=0, lambda=0)
predict_y<-predict(OLS_model,newx=X.test)
plot(x=Y.test,y=predict_y,xlab="Actual",ylab="Prediction",
     main="Actual VS Prediction(OLS)")
trainlinmod<-lm(predict_y~Y.test)
abline(trainlinmod,col=3,lwd=2.5,lty="solid")
abline(a=0,b=1,col=2,lwd=2.5,lty="dashed")
legend("topleft",
       legen=c("Model","Base"),
       col=c("green","red"),text.width = 0.5,bty="n",
       lwd=2.5,lty=c("solid","dashed"))

##PCR
PCR_model<-pcr(Y.train~X.train,scale=FALSE,validation="CV",segments=10)# Cross Validation
cverr<-RMSEP(PCR_model)$val[1,,]
predict_y<-predict(PCR_model,X.test,ncomp=which.min(cverr)-1)
plot(cverr,type="l",ylab="mean squared error of prediction",xlab="number of components")
abline(v=which.min(cverr)-1,lwd=1,col=2)
plot(x=Y.test,y=predict_y,xlab="Actual",ylab="Prediction",
     main="Actual VS Prediction(PCR)")
trainlinmod<-lm(predict_y~Y.test)
abline(trainlinmod,col=3,lwd=2.5,lty="solid")
abline(a=0,b=1,col=2,lwd=2.5,lty="dashed")
legend("topleft",
       legen=c("Model","Base"),
       col=c("green","red"),text.width = 0.5,bty="n",
       lwd=2.5,lty=c("solid","dashed"))

##Ramdom Forest
rf.train<-randomForest(Y.train~.,data=X.train,importance=TRUE,ntree=500)
plot(rf.train,main="ERROR & TREES")
which.min(rf.train$mse)
rf.pre<-predict(rf.train,newdata=X.test)
plot(x=Y.test,y=rf.pre,xlab="Actual",ylab="Prediction",
     main="Actual VS Prediction(RandomForest)")
trainlinmod<-lm(rf.pre~Y.test)
abline(trainlinmod,col=3,lwd=2.5,lty="solid")
abline(a=0,b=1,col=2,lwd=2.5,lty="dashed")
legend("topleft",
       legen=c("Model","Base"),
       col=c("green","red"),text.width = 0.5,bty="n",
       lwd=2.5,lty=c("solid","dashed"))
imp.geno<-rf.train$importance
head(imp.geno)
plot(imp.geno[,2],type="l")

##GBLUP
ans.GAUSS<-kinship.BLUP(y=Y.train,G.train=X.train,G.pred = X.test,K.method = "GAUSS")
blup.pre<-ans.GAUSS$g.pred
plot(x=Y.test,y=blup.pre,xlab="Actual",ylab="Prediction",
     main="Actual VS Prediction(GBLUP)")
trainlinmod<-lm(blup.pre~Y.test)
abline(trainlinmod,col=3,lwd=2.5,lty="solid")
abline(a=0,b=1,col=2,lwd=2.5,lty="dashed")
legend("topleft",
       legen=c("Model","Base"),
       col=c("green","red"),text.width = 0.5,bty="n",
       lwd=2.5,lty=c("solid","dashed"))