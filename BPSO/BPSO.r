##############Data road#########################################################
X <- read.csv("C:/Users/user/Box/thyrsocope/NGS 楷备/data/21_03/X.csv")
y <- read.csv("C:/Users/user/Box/thyrsocope/NGS 楷备/data/21_03/y.csv")

str(X) #9471 features
unique(y) #four classes
y_ordinal = factor(y$x,levels=c("normal prostate tissue","normal prostate adjacent to tumor","primary prostate tumor","metastatic prostate tumor"))
################################################################################

##############library & function################################################
source('C:/Users/user/Box/thyrsocope/NGS 楷备/code/21_06/CFS_0702.R', echo=TRUE)
source('C:/Users/user/Box/thyrsocope/NGS 楷备/code/21_06/iBPSO.R', echo=TRUE)
library('BPSO') #https://rdrr.io/github/DominikMueller64/BPSO/man/bpsoptim.html
library('splitTools')
library('klaR')
library('pROC')
################################################################################

###############################practice#########################################
#CFS
CFS_features <- CFS_0702(X,as.numeric(y_ordinal))
X_CFS = X[,which(colnames(X) %in% CFS_features)]

#BPSO objective function
fn <- function(position) {
  train_bpso <- train_X[,seq(dim(X_CFS)[2])[position]]
  valid_bpso <- valid_X[,seq(dim(X_CFS)[2])[position]]
  nb <- NaiveBayes(train_bpso,train_Y)
  pred_y <- predict(nb, valid_bpso)$class
  if(length(unique(pred_y))!=4){
    mauc <- 0.5
  }else{
    mauc <- multiclass.roc(valid_Y,as.numeric(pred_y),direction='<',levels=c("normal prostate tissue","normal prostate adjacent to tumor","primary prostate tumor","metastatic prostate tumor"))$auc
  }
  return(mauc)
}


#output
result <- c()
features <- vector(mode = "list", length = 10)

# 10 iterations (BPSO)
for(k in 1:10){
  
  #Data split
  inds <- partition(y_ordinal, p = c(train = 0.6, valid = 0.2, test = 0.2))
  train_X <- X_CFS[inds$train,]
  train_Y <- y_ordinal[inds$train]
  valid_X <- X_CFS[inds$valid,]
  valid_Y <- y_ordinal[inds$valid]
  test_X <- X_CFS[inds$test,]
  test_Y <- y_ordinal[inds$test]
  
  #BPSO
  fm <- bpsoptim(par = rep(FALSE,dim(X_CFS)[2]),fn=fn,control=list(s=60,maxit=100,c.p=2,c.g=2))
  bpso_par <- fm$par
  
  #test
  Train_X <- X_CFS[c(inds$train,inds$valid),seq(dim(X_CFS)[2])[bpso_par]]
  Train_Y <- y_ordinal[c(inds$train,inds$valid)]
  nb <- NaiveBayes(Train_X,Train_Y)
  pred_y <- predict(nb, test_X[,seq(dim(X_CFS)[2])[bpso_par]])$class
  result <- c(result,multiclass.roc(test_Y,as.numeric(pred_y),direction='<',levels=c("normal prostate tissue","normal prostate adjacent to tumor","primary prostate tumor","metastatic prostate tumor"))$auc)
  
  # save selected features
  features[[k]] <- colnames(X_CFS)[bpso_par]
}

print(c(mean(result),sd(result)))
sapply(features,length)


#output
result <- c()
features <- vector(mode = "list", length = 10)

# 10 iterations (iBPSO)
for(k in 1:10){
  
  #Data split
  inds <- partition(y_ordinal, p = c(train = 0.6, valid = 0.2, test = 0.2))
  train_X <- X_CFS[inds$train,]
  train_Y <- y_ordinal[inds$train]
  valid_X <- X_CFS[inds$valid,]
  valid_Y <- y_ordinal[inds$valid]
  test_X <- X_CFS[inds$test,]
  test_Y <- y_ordinal[inds$test]
  
  #BPSO
  fm <- ibpso(par = rep(FALSE,dim(X_CFS)[2]),fn=fn,control=list(s=60,maxit=100,c.p=2,c.g=2))
  bpso_par <- fm$par
  
  #test
  Train_X <- X_CFS[c(inds$train,inds$valid),seq(dim(X_CFS)[2])[bpso_par]]
  Train_Y <- y_ordinal[c(inds$train,inds$valid)]
  nb <- NaiveBayes(Train_X,Train_Y)
  pred_y <- predict(nb, test_X[,seq(dim(X_CFS)[2])[bpso_par]])$class
  result <- c(result,multiclass.roc(test_Y,as.numeric(pred_y),direction='<',levels=c("normal prostate tissue","normal prostate adjacent to tumor","primary prostate tumor","metastatic prostate tumor"))$auc)
  
  # save selected features
  features[[k]] <- colnames(X_CFS)[bpso_par]
}

print(c(mean(result),sd(result)))
sapply(features,length)