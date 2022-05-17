##############Data road#########################################################
X <- read.csv("C:/Users/user/Box/thyrsocope/NGS 연구/data/21_03/X.csv")
y <- read.csv("C:/Users/user/Box/thyrsocope/NGS 연구/data/21_03/y.csv")

str(X) #9471 features
unique(y) #four classes
y_ordinal = factor(y$x,levels=c("normal prostate tissue","normal prostate adjacent to tumor","primary prostate tumor","metastatic prostate tumor"))

################################################################################

##############library & function################################################
source('C:/Users/user/Box/thyrsocope/NGS 연구/code/21_04/CFS.R', echo=TRUE)
library(ABCoptim)
library(caret) #createFolds
library(e1071) #svm
library(pROC)
assign("last.warning", NULL, envir = baseenv())
################################################################################


###############################practice#########################################
CFS_y = as.numeric(y_ordinal=="normal prostate tissue")
normal_features_CFS = CFS(50,X,CFS_y)
length(which(colnames(X) %in% normal_features_CFS))
X_CFS = X[,which(colnames(X) %in% normal_features_CFS)]

CV_list = createFolds(CFS_y,k=3)
fun <- function(x){
  if(sum(round(x)==1)==0){
    return(0)
  }
  
  auc <- c()
  for(i in 1:3){
    train_X = X_CFS[-(CV_list[[i]]),which(round(x)==1)]
    test_X = X_CFS[(CV_list[[i]]),which(round(x)==1)]
    
    train_y = CFS_y[-(CV_list[[i]])]
    test_y = CFS_y[(CV_list[[i]])]
    
    model <- svm(x=train_X,y=train_y,type="C-classification")
    predict_y = predict(model,test_X)
    if(sum(predict_y==1)==length(test_y) | sum(predict_y==0)==length(test_y)){
      auc <- c(auc,0.5)
    }else{
      auc <- c(auc,roc(as.numeric(predict_y),as.numeric(test_y),quiet = TRUE)$auc)
    }
  }
  
  print(mean(auc))
  return(-mean(auc))
}

optim = abc_optim(rep(1,dim(X_CFS)[2]), fun, lb=rep(0,dim(X_CFS)[2]), ub=rep(1,dim(X_CFS)[2]), criter=3, maxCycle = 8000)
plot(optim)
round(optim$par)

auc <- c()
for(i in 1:3){
  train_X = X_CFS[-(CV_list[[i]]),which(round(optim$par)==1)]
  test_X = X_CFS[(CV_list[[i]]),which(round(optim$par)==1)]
  
  train_y = CFS_y[-(CV_list[[i]])]
  test_y = CFS_y[(CV_list[[i]])]
  
  model <- svm(x=train_X,y=train_y,type="C-classification")
  predict_y = predict(model,test_X)
  if(sum(predict_y==1)==length(test_y) | sum(predict_y==0)==length(test_y)){
    auc <- c(auc,0.5)
  }else{
    auc <- c(auc,roc(as.numeric(predict_y),as.numeric(test_y),quiet = TRUE)$auc)
  }
}
auc
################################################################################

############################Real Application####################################
CV_list = createFolds(as.numeric(y_ordinal=="normal prostate tissue"),k=10)

auc <- c()

fun_ABC <- function(x){
  if(sum(round(x)==1)==0){
    return(0)
  }
  
  auc <- c()
  for(i in 1:3){
    train_X_ABC = X_CFS[-(CV_list_ABC[[i]]),which(round(x)==1)]
    cv_X = X_CFS[(CV_list_ABC[[i]]),which(round(x)==1)]
    
    train_y_ABC = train_y[-(CV_list_ABC[[i]])]
    cv_y = train_y[(CV_list_ABC[[i]])]
    
    model <- svm(x=train_X_ABC,y=train_y_ABC,type="C-classification")
    predict_y = predict(model,cv_X)
    if(sum(predict_y==1)==length(cv_y) | sum(predict_y==0)==length(cv_y)){
      auc <- c(auc,0.5)
    }else{
      auc <- c(auc,roc(as.numeric(predict_y),as.numeric(cv_y),quiet = TRUE)$auc)
    }
  }
  
  print(mean(auc))
  return(-mean(auc))
}

for(k in 1:10){
  train_X = X[-(CV_list[[i]]),]
  
  train_y = as.numeric(y_ordinal=="normal prostate tissue")[-(CV_list[[i]])]
  test_y = as.numeric(y_ordinal=="normal prostate tissue")[(CV_list[[i]])]
  
  normal_features_CFS = CFS(5000,train_X,train_y)
  X_CFS = train_X[,which(colnames(X) %in% normal_features_CFS)]
  test_X = X[(CV_list[[i]]),which(colnames(X) %in% normal_features_CFS)]
  
  CV_list_ABC = createFolds(train_y,k=3)
  
  optim = abc_optim(rep(1,dim(X_CFS)[2]), fun_ABC, lb=rep(0,dim(X_CFS)[2]), ub=rep(1,dim(X_CFS)[2]), criter=3)
  colnames(X_CFS)[which(round(optim$par)==1)]
    
  model <- svm(x=X_CFS[,which(round(optim$par)==1)],y=train_y,type="C-classification")
  predict_y = predict(model,test_X[,which(round(optim$par)==1)])
  
  if(sum(predict_y==1)==length(test_y) | sum(predict_y==0)==length(test_y)){
    auc <- c(auc,0.5)
  }else{
    auc <- c(auc,roc(as.numeric(predict_y),as.numeric(test_y),quiet = TRUE)$auc)
  }
}
