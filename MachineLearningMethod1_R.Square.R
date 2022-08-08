library(pls)
library(glmnet)
library(randomForest)
library(skimr)
library(DataExplorer)
library(e1071)
library(caret)
library(rrBLUP)
set.seed(123)
setwd("e:/rdata")
load("pheno_raw.Rdata") 
load("cross.Rdata")



##store the yeast genotype data in a matrix
yeast_geno<-matrix(NA,nrow=1008,ncol=11623)
int1<-1
int2<-0
for(i in 1:16){
  int2<-int2+length(cross$geno[[i]]$data)/1008
  yeast_geno[,int1:int2]<-as.matrix(cross$geno[[i]]$data)
  int1<-int1+length(cross$geno[[i]]$data)/1008
}

Y_raw<-c(cross$pheno[[i]])
skim(Y_raw)
plot_missing(Y_raw)



##Ridge regression:Alpha=0
lambdas<-10^seq(3,-4,by=-.1)
missings<-c(NA,length=46)
R2_ridge<-matrix(NA,nrow=46,ncol=10)
for(i in 1:46){
for(k in 1:10){
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
cv_fit<-cv.glmnet(X.train,Y.train,alpha=0,lambda=lambdas,nfolds=10)# Cross-validation
lambda<-cv_fit$lambda.min #the optimal lambda
RR_model<-glmnet(X.train,Y.train,alpha=0, lambda=lambda)
#calculate R-square
y_bar<-mean(Y.test)
predict_y<-predict(RR_model,s=lambda,newx=X.test)
R2_ridge[i,k]<-1-sum((predict_y-y_bar)^2)/sum((Y.test-y_bar)^2)
}}
plot(cv_fit)
rr.model<-glmnet(X.test,Y.test,alpha=0,lambda=lambdas)
plot(rr.model,xvar='lambda')




##Lasso regression:Alpha=1
R2_lasso<-matrix(NA,nrow=46,ncol=10)
for(i in 1:46){
  for(k in 1:10){
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
  cv_fit<-cv.glmnet(X.train,Y.train,alpha=1,lambda=lambdas,nfolds=10)#Cross-validation
  lambda<-cv_fit$lambda.min #the optimal lambda
  LA_model<-glmnet(X.train,Y.train,alpha=1, lambda=lambda)
  #calculate R-square on test set
  y_bar<-mean(Y.test)
  predict_y<-predict(LA_model,s=lambda,newx=X.test)
  R2_lasso[i,k]<-1-sum((predict_y-y_bar)^2)/sum((Y.test-y_bar)^2)
}}
R2_lasso
i=2
plot(cv_fit)
lasso.models<-glmnet(X.test,Y.test,alpha=1,lambda=lambdas)
plot(lasso.models,xvar='lambda')

##Elastic net: Alpha between 0 to 1
Alphas<-seq(0,1,by=0.05)
R2_elastic<-matrix(NA,nrow=46,ncol=10)
alphas<-matrix(NA,nrow=46,ncol=10)

for(g in 1:10){
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
  #CV
  cv_fit<-cv.glmnet(X.train,Y.train,alpha=Alphas[k],lambda=lambdas,nfolds=10)#Cross-validation
  lambda<-cv_fit$lambda.min #the optimal lambda
  EN_model<-glmnet(X,Y,alpha=Alphas[k],lambda=lambda)
  #calculate R-square
  y_bar<-mean(Y.test)
  predict_y<-predict(EN_model,s=lambda,newx=X.test)
  R2_Alphas[k]<-1-sum((predict_y-y_bar)^2)/sum((Y.test-y_bar)^2)
  }
  alphas[i,g]<-which.max(R2_Alphas)
  R2_elastic[i,g]<-max(R2_Alphas)
}}
plot(Alphas,R2_Alphas,type="l",xlab="alpha",ylab="R square")
abline(h=max(R2_Alphas),col=2,lty=2)
abline(v=Alphas[which.max(R2_Alphas)],col=2,lty=2)

cv_fit<-cv.glmnet(X.train,Y.train,alpha=Alphas[2],lambda=lambdas,nfolds=10)
plot(cv_fit)

en.models<-glmnet(X.test,Y.test,alpha=Alphas[2],lambda=lambdas)
plot(en.models,xvar='lambda')

Alphas[alphas]# the values of alpha, 需要修改
R2_elastic





## Principal component regression:PCR
npcs<-matrix(NA,nrow=46,ncol=10)#record the number of PCs 
R2_PCR<-matrix(NA,nrow=46,ncol=10)
for(g in 1:10){
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
  PCR_model<-pcr(Y.train~X.train,scale=FALSE,validation="CV",segments=10)# Cross Validation
  cverr<-RMSEP(PCR_model)$val[1,,]
  npcs[i,g]<-which.min(cverr)-1
  predict_y<-predict(PCR_model,X.test,ncomp=npcs[i,g])
  y_bar<-mean(Y.test)
  R2_PCR[i,g]<-1-sum((predict_y-y_bar)^2)/sum((Y.test-y_bar)^2)
}}

PCR_model<-pcr(Y.train~X.train,scale=FALSE,validation="CV",segments=10)
validationplot(PCR_model,val.type="MSEP") #使用 validatìonplot ()函数作出交叉验证得分的图像，MSE是均方误差
validationplot(PCR_model,val.type="R2")
validationplot(cv_fit,val.type="MSEP")
validationplot?
plot(cverr,type="l",ylab="mean squared error of prediction",xlab="number of components")
  abline(v=which.min(cverr)-1,lyt=2,col=2)
plot(RMSEP(PCR_model), legendpos = "topright")
plot(PCR_model, method="onesigma")
?pcr
predict_y<-predict(PCR_model,X.test,ncomp=110)
y_bar<-mean(Y.test)
1-sum((predict_y-y_bar)^2)/sum((Y.test-y_bar)^2)
cverr<-RMSEP(PCR_model)$val[1,,]
npcs<-which.min(cverr)-1
# true values and predictions on the test set
plot(Y.test,type="l",xlab="sample",ylab="yeast population growth")
lines(predict_y,col=2)
legend("topright",pch=c(15,15),legend=c("True values","Predictions"),col=c(1,2),bty="n")




##OLS
R2_OLS<-matrix(NA,nrow=46,ncol=10)
for(g in 1:10){
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
  OLS_model<-glmnet(X.train,Y.train,alpha=0, lambda=0)
  predict_y<-predict(OLS_model,newx=X.test)
  y_bar<-mean(Y.test)
  R2_OLS[i,g]<-1-sum((predict_y-y_bar)^2)/sum((Y.test-y_bar)^2)
}
}

##Random Forest

rf.train<-randomForest(Y.train~.,data=X.train,importance=TRUE,ntree=500)
#树的数量与Error的关系
plot(rf.train,main="ERROR & TREES")
which.min(rf.train$mse)#mse与ntrees的关系
#变量重要性
importance(rf.train)
varImpPlot(rf.train,main="Variable Importance Plot",type=1)
varImpPlot(rf.train,main="Variable Importance Plot",type=2)
varImpPlot(rf.train,main="Variable Importance Plot")
#偏依赖性 太多了 应该不看了
partialPlot(x=rf.train,
            pred.data=)
#测试集的预测
rf.pre<-predict(rf.train,newdata=X.test)
rf.pre
plot(Y.test,type="l")
#另一种检验预测结果的方法
plot(x=Y.test,y=rf.pre,
     xlab="Actual",ylab="Prediction",
     main="Actual VS Prediction")
trainlinmod<-lm(rf.pre~Y.test)
abline(trainlinmod,col=3,lwd=2.5,lty="solid")
abline(a=0,b=1,col=2,lwd=2.5,lty="dashed")
legend("topleft",
       legen=c("Model","Base"),
       col=c("green","red"),text.width = 0.5,bty="n",
       lwd=2.5,lty=c("solid","dashed")
       )

# R square on test set

y_bar<-mean(rf.pre)
1-sum((rf.pre-y_bar)^2)/sum((Y.test-y_bar)^2)



plot(rf.train$y,rf.pre,main="Test Set")
abline(1,1)
##得到的效果不理想 要剔除一部分因变量
##做筛选,查看importance
imp.geno<-rf.train$importance
head(imp.geno)
##作图查看前30个因变量
varImpPlot(rf.train, n.var = min(30, nrow(rf.train$importance)),
           main = 'Top 30 - variable importance')
# %IncMSE”即increase in mean squared error，通过对每一个预测变量随机赋值
#如果该预测变量更为重要，那么其值被随机替换后模型预测的误差会增大。
#因此，该值越大表示该变量的重要性越大

# IncNodePurity”即increase in node purity，通过残差平方和来度量，
#代表了每个变量对分类树每个节点上观测值的异质性的影响，
#从而比较变量的重要性。该值越大表示该变量的重要性越大。

#对于“%IncMSE”或“IncNodePurity”，二选一作为判断预测变量重要性的指标。
#需注意的是，二者的排名存在一定的差异。
## Cross Validation
R2_RF<-matrix(NA,nrow=46,ncol=10)
for(i in 1:46){
  for(k in 1:10){
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
    
    
    
    
    #calculate R-square
    y_bar<-mean(Y.test)
    predict_y<-predict(RR_model,s=lambda,newx=X.test)
    R2_RF[i,k]<-1-sum((predict_y-y_bar)^2)/sum((Y.test-y_bar)^2)
  }}
## 用函数交叉验证 mtrys& ntree
mtrys<-c(5811,3874,2905)
mses<-matrix(NA,nrow=500,ncol=3)
R2_RF<-matrix(NA,nrow=46,ncol=10)

  for(i in 1:46){
    for(j in 1:10){
      for(k in 1:3){
        mtry<-mtrys[k]
      #split dataset
      Y_raw<-c(cross$pheno[[i]])
      nn<-c(which(Y_raw>-1000))## Exclusion of missing values
      missings[i]<-1008-length(nn)
      Y<-Y_raw[nn] #pheotype data without missing values
      X<-yeast_geno[nn,]
      num<-sample(length(Y),0.8*length(Y))
      Y.train<-Y[num]
      X.train<-X[num,]
      Y.test<-Y[-num]
      X.test<-X[-num,]
      ##train the model by using train set 
      rf.train<-randomForest(Y.train~.,data=X.train,ntree=500,mtry=mtry)
      mses[,k]<-rf.train$mse
      }
     a<- which(mses==min(mses),arr.ind=TRUE)
     ntree<-a[1,1]
     mtry<-mtrys[a[1,2]]
     rf.train<-randomForest(Y.train~.,data=X.train,ntree=ntree,mtry=mtry)
     rf.pre<-predict(rf.train,newdata=X.test)
    
    }}

##rfcv减少特征的数量

## 不做参数选择
R2_RF<-matrix(NA,nrow=46,ncol=10)

for(i in 1:46){
  for(j in 1:10){
    #split dataset
    Y_raw<-c(cross$pheno[[i]])
    nn<-c(which(Y_raw>-1000))## Exclusion of missing values
    Y<-Y_raw[nn] #pheotype data without missing values
    X<-yeast_geno[nn,]
    num<-sample(length(Y),0.8*length(Y))
    Y.train<-Y[num]
    X.train<-X[num,]
    Y.test<-Y[-num]
    X.test<-X[-num,]
    ##train the model by using train set 
    rf.train<-randomForest(Y.train~.,data=X.train,ntree=500)
    rf.pre<-predict(rf.train,newdata=X.test)
    y_bar<-mean(rf.pre)
    R2_RF[i,j]<-1-sum((rf.pre-y_bar)^2)/sum((Y.test-y_bar)^2)
    }
}



##
R2s<-cbind(R2_ridge,R2_lasso,R2_elastic,R2_PCR)



## Save the results
save(R2_ridge,file="R2_ridge.RData")
save(R2_lasso,file="R2_lasso.RData")
save(R2_elastic,file="R2_elastic.RData")
Alpha<-Alphas[alphas]
save(Alpha,file="alphas.RData")
save(R2_PCR,file="R2_PCR.RData")
save(R2_OLS,file="R2_OLS.RData")
save(R2_RF,file="R2_RF.RData")
rowMeans(R2_PCR)
results<-cbind(OLS=c(rowMeans(R2_OLS)),PCR=c(rowMeans(R2_PCR)),Ridge=c(rowMeans(R2_ridge)),
               Lasso=c(rowMeans(R2_lasso)),Elastic=c(rowMeans(R2_elastic)))
write.csv(results,file="results.csv")
