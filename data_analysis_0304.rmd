##############Data road#########################################################
X <- read.csv("C:/Users/user/Box/thyrsocope/NGS 楷备/data/21_03/X.csv")
y <- read.csv("C:/Users/user/Box/thyrsocope/NGS 楷备/data/21_03/y.csv")

str(X) #9471 features
unique(y) #four classes
y_ordinal = factor(y$x,levels=c("normal prostate tissue","normal prostate adjacent to tumor","primary prostate tumor","metastatic prostate tumor"))

################################################################################

##############mRMR##############################################################
#mRMR
library(praznik)
library(randomForest)
library(pROC)

#practice
MRMR(X,y_ordinal,k=34)

library(caret)
CV_list = createFolds(y_ordinal,k=10)

k_decision <- c()
for(k in 10:50){
  auc <-c()
  for(i in 1:length(CV_list)){
    train_X = X[-(CV_list[[i]]),]
    test_X = X[(CV_list[[i]]),]
    
    train_y = y_ordinal[-(CV_list[[i]])]
    test_y = y_ordinal[(CV_list[[i]])]
    
    mRMR <- MRMR(as.data.frame(train_X),train_y,k=k)
    selected_features <- mRMR$selection
    
    #colnames(X_normalization) %in% selected_features
    
    selected_train_X <- train_X[,selected_features]
    selected_test_X <- test_X[,selected_features]
    
    #library(randomForest)
    rf <- randomForest(selected_train_X,train_y)
    pred.y <- predict(rf,selected_test_X)
    
    #library(pROC)
    auc <- c(auc,multiclass.roc(test_y,as.numeric(pred.y),direction='<',levels=c("normal prostate tissue","normal prostate adjacent to tumor","primary prostate tumor","metastatic prostate tumor"))$auc)
    
  }
  
  k_decision <- rbind(k_decision,data.frame(auc.avg=mean(auc),k=k))
}

plot(k_decision$k,k_decision$auc.avg,type='b',lty=2,lwd=3,xlab='# features',ylab="accuracy")

########################################
#LARS

library(caret)
library(dplyr)

k_decision <- c()
for(lars in 5:50){
  auc <-c()
  for(i in 1:length(CV_list)){
    train_X = X[-(CV_list[[i]]),]
    test_X = X[(CV_list[[i]]),]
    
    train_y = y_ordinal[-(CV_list[[i]])]
    test_y = y_ordinal[(CV_list[[i]])]
    
    mRMR <- MRMR(as.data.frame(train_X),train_y,k=300)
    selected_features <- mRMR$selection
    
    selected_train_X <- train_X[,selected_features]
    
    df = data.frame(cbind(selected_train_X, class = as.numeric(train_y)))
    rpartmod = train(class~.,data=df, method="lars2")
    rrfImp <- varImp(rpartmod)$importance %>% arrange(desc(Overall))
    
    #library(randomForest)
    rf <- randomForest(selected_train_X[colnames(selected_train_X) %in% rownames(rrfImp)[1:lars]],train_y)
    pred.y <- predict(rf,test_X[colnames(test_X) %in% rownames(rrfImp)[1:lars]])
    
    #library(pROC)
    auc <- c(auc,multiclass.roc(test_y,as.numeric(pred.y),direction='<',levels=c("normal prostate tissue","normal prostate adjacent to tumor","primary prostate tumor","metastatic prostate tumor"))$auc)
    
  }
  
  k_decision <- rbind(k_decision,data.frame(auc.avg=mean(auc),k=lars))
}

plot(k_decision$k,k_decision$auc.avg,type='b',lty=2,lwd=3,xlab='# features',ylab="accuracy")

##########################################
#lasso

library(glmnet)
library(coefplot)

binary_y = as.numeric(y_ordinal=="normal prostate tissue")
# Fit the LASSO model (Lasso: Alpha = 1)
# The best lambda value is stored inside 'cv.lasso$lambda.min'.
cv.lasso <- cv.glmnet(as.matrix(X), binary_y, family='binomial', alpha=1, parallel=TRUE, nfolds=10,standardize=FALSE, type.measure='auc')

coef_lasso <- extract.coef(cv.lasso)



###############################################
#Fast OSFS
source('C:/Users/user/Box/thyrsocope/NGS 楷备/code/货 弃歹/partial_corr_coef.R', echo=TRUE)
source('C:/Users/user/Box/thyrsocope/NGS 楷备/code/货 弃歹/my_cond_indep_fisher_z.R', echo=TRUE)
source('C:/Users/user/Box/thyrsocope/NGS 楷备/code/货 弃歹/computer_dep_2.R', echo=TRUE)
source('C:/Users/user/Box/thyrsocope/NGS 楷备/code/货 弃歹/optimal_compter_dep_2.R', echo=TRUE)
source('C:/Users/user/Box/thyrsocope/NGS 楷备/code/货 弃歹/fast_osfs_z.R', echo=TRUE)

#practice
#data1 = cbind(X_normalization,(as.numeric(y)))
#class_index = ncol(data1)
#alpha = 0.05

#fast_osfs_results2 <- fast_osfs_z(data1,ncol(data1),0.05)

library(caret)
library(randomForest)
library(pROC)

CV_list = createFolds(y_ordinal,k=10)

auc <-c()
selected_feature_set <- c()
for(i in 1:length(CV_list)){
  train_X <- X[-(CV_list[[i]]),]
  train_y <- y_ordinal[-(CV_list[[i]])]
  train_df = cbind(train_X,train_y)
  test_X = X[(CV_list[[i]]),]
  test_y = y_ordinal[(CV_list[[i]])]
  
  class_index = ncol(train_df)
  alpha = 0.05
  
  selected_features <- fast_osfs_z(train_df,class_index,0.05)
  
  #colnames(X_normalization) %in% selected_features
  
  selected_train_X <- train_X[,selected_features]
  selected_test_X <- test_X[,selected_features]
  
  #library(randomForest)
  rf <- randomForest(selected_train_X,train_y)
  pred.y <- predict(rf,selected_test_X)
  table(pred.y,test_y)
  
  #library(pROC)
  auc <- c(auc,multiclass.roc(as.numeric(pred.y),as.numeric(test_y))$auc)
  selected_feature_set[[i]] <- selected_features
  
}
mean(auc)
sd(auc)
feature_len <- sapply(selected_feature_set,length)
mean(feature_len)
sd(feature_len)

####################################################
#DRPT
library(corpcor) #pesudoinvers
library(pracma) #rands
library(prospectr) #savitzkyGolay
library(randomForest)
library(pROC)
library(caret)

allF = ncol(X_normalization)
t= 50
clusters= 50
runIter = 10
featuresPicked = matrix(0,nrow=clusters,runIter)
max_acc = matrix(0,nrow=runIter,ncol=2)

df = X_normalization
y_ordinal = factor(y,levels=c("normal prostate tissue","normal prostate adjacent to tumor","primary prostate tumor","metastatic prostate tumor"))
CV_list = createFolds(y_ordinal,k=10)

for(counter in 1:runIter){
  #data partition
  #data
  A = df[-CV_list[[counter]],]
  C = A
  r = nrow(A)
  B = as.numeric(y_ordinal[-CV_list[[counter]]])
  
  # A= train_data
  iA = pseudoinverse(A) # A should be a matrix
  X = iA%*%B
  outliersLen = allF
  listX = matrix(0,nrow=5,ncol=length(X))
  cleanedF = matrix(0,nrow=5,ncol=length(X))
  listX[1,] = X
  ii=1
  
  while(outliersLen > (allF*0.21)){
    outliers = isoutlier(abs(listX[ii,1:outliersLen]))
    tmp = which(outliers)
    outliersLen = length(tmp)
    cleanedF[ii,1:outliersLen] =tmp
    ii = ii+1
    listX[ii,1:outliersLen] = listX[(ii-1),tmp]
  }
  
  if(outliersLen<10){
    TF = isoutlier(abs(X))
    outliersLen = sum(TF)
  }
  
  uniqueX = sort(unique(abs(X)),decreasing = TRUE)
  threshold = mean(uniqueX[1:outliersLen*10])
  cleanedF =1:allF
  
  while(length(cleanedF) > (outliersLen*(2/(ii-1)))){
    irrF = which(abs(X) < threshold)
    irrF = cbind(irrF,X[irrF])
    irrF = irrF[order(-irrF[,2]),]
    irrF = irrF[,1]
    threshold = threshold*1.03
    cleanedF = setxor((1:allF),irrF)
  }
  
  A = A[,cleanedF]
  C = C[,cleanedF]
  allF = length(cleanedF)
  
  #========================Perturbantion matrix======================
  svdA = svd(A)$d
  smallestAan = min(svdA)
  iA = pseudoinverse(A)
  X = iA%*%B
  minPer = apply(A,2,min)*1e-3*smallestAan
  maxPer = apply(A,2,max)*1e-2*smallestAan
  mError = 1e-3*smallestAan
  perVal = matrix(0,nrow=r,ncol=allF)
  px = matrix(0,nrow=r,ncol=allF)
  for(z in 1:t){
    perVal = mError*rands(n=r,N=allF-1)
    nr = norm(perVal,type="2")
    pA = A + perVal
    piA = pseudoinverse(pA)
    Xtilda = piA%*%B
    DX = abs(Xtilda - X)
    px[z,]=DX
  }
  
  pX = apply(px,2,mean)
  ent = -apply(C*log(C),2,function(x) sum(x,na.rm=TRUE))##error
  
  #====================Sorting PX and Finding top ranked features====================
  pX = savitzkyGolay(pX,m=2,p=3,w=5)
  uniquePX = length(unique(c(pX)))
  roundMetric = 20
  roundedPX = pX 
  
  while(uniquePX >50){
    roundedPX=round(pX,roundMetric)
    roundMetric = roundMetric -1
    uniquePX = length(unique(c(roundedPX)))
  }
  
  out=c()
  uniquekeys = unique(c(roundedPX))
  indexOut = 1
  
  for(key in 1:uniquePX){
    selectedPX = which(roundedPX==uniquekeys[key])
    filteredEnt = ent[selectedPX]
    uniqueEnt = length(unique(filteredEnt))
    roundMetric = 5
    roundedEnt = filteredEnt
    while(uniqueEnt > 20){
      roundedEnt = round(filteredEnt,roundMetric)
      roundMetric = roundMetric -1
      uniqueEnt = length(unique(roundedEnt))
    }
    lenEnt= length(unique(roundedEnt))
    uniqueEnt = unique(roundedEnt)
    for(index in 1:lenEnt){
      selectedEnt = which(roundedEnt==uniqueEnt[index])
      filteredX = X[selectedPX[selectedEnt]]
      rankFilteredX = cbind(selectedPX[selectedEnt],filteredX)
      sortedFilteredX = rankFilteredX[order(-rankFilteredX[,2]),]
      out[indexOut] = sortedFilteredX[1]
      indexOut = indexOut + 1
    }
  }
  outEnt = ent[out]
  rankOutEnt = cbind(out,outEnt)
  sortedRankOutEnt = rankOutEnt[order(rankOutEnt[,2]),]
  sortedRankOutEnt = cbind((1:length(sortedRankOutEnt[,1])),sortedRankOutEnt)
  outX = X[out,1]
  rankOutX = cbind(out,outX)
  sortedRankOutX = rankOutX[order(-rankOutX[,2]),]
  sortedRankOutX = cbind((1:length(sortedRankOutX[,1])),sortedRankOutX)
  sortedRankOutEnt = cbind(sortedRankOutEnt,rep(NA,nrow(sortedRankOutEnt)))
  for(i in 1:length(sortedRankOutX[,1])){
    indexX = which(sortedRankOutEnt[i,2] == sortedRankOutX[,2])
    sortedRankOutEnt[i,4] = indexX+sortedRankOutEnt[i,1]
  }
  sortedRankOutEnt = sortedRankOutEnt[order(sortedRankOutEnt[,4]),]
  out = cleanedF[sortedRankOutEnt[,2]]
  upperBand = min(clusters,length(out))
  
  # classification
  best_result = matrix(0,nrow=3,ncol=2)
  m_acc= 0
  o_f = 0
  for(k in 2:upperBand){
    centres = out[1:k]
    rf <- randomForest(df[-CV_list[[counter]],centres],y_ordinal[-CV_list[[counter]]])
    pred.y <- predict(rf,df[CV_list[[counter]],centres])
    acc= multiclass.roc(as.numeric(pred.y),as.numeric(y_ordinal[CV_list[[counter]]]))$auc
    if(acc > m_acc){
      m_acc = acc
      o_f = k
      featuresPicked[1:k,counter] = out[1:k]
    }
    if(m_acc==1){
      break
    }
  }
  max_acc[counter,1] = m_acc
  max_acc[counter,2] = o_f
  cat(" SF = ",o_f," , CA =", m_acc)
}

apply(max_acc,2,mean)
apply(max_acc,2,sd)

