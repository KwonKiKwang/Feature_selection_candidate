##############Data road#########################################################
X <- read.csv("C:/Users/user/Box/thyrsocope/NGS 楷备/data/21_03/X.csv")
y <- read.csv("C:/Users/user/Box/thyrsocope/NGS 楷备/data/21_03/y.csv")

str(X) #9471 features
unique(y) #four classes
y_ordinal = factor(y$x,levels=c("normal prostate tissue","normal prostate adjacent to tumor","primary prostate tumor","metastatic prostate tumor"))
################################################################################


##############library & function################################################
source('C:/Users/user/Box/thyrsocope/NGS 楷备/code/21_07/fisher_score.R')
source('C:/Users/user/Box/thyrsocope/NGS 楷备/code/21_07/GSFSJNE_fn.R')
library(dplyr)
library(stats)
library('splitTools')
library('klaR')
library('pROC')
################################################################################

###############################practice#########################################
# essential to standardization for distance based learning
X_sdd = scale(X)

#fisher score
fs <- fisher_score(X_sdd,y_ordinal)
fs_col <- data.frame(colnames(X_sdd),fs) %>% arrange(desc(fs))
fs_selected_cols <- fs_col[1:200,1]# 200 : hyperparameter

#GSFSJNE
GSFSJNE(X=X_sdd,y=y_ordinal,S=fs_selected_cols,delta=0.3) #0.3:hyperparameter

#output
result <- c()
features <- vector(mode = "list", length = 10)

# 10 iterations (BPSO)
for(k in 1:10){
  
  #Data split
  inds <- partition(y_ordinal, p = c(train = 0.6, valid = 0.2, test = 0.2))
  train_X <- X_sdd[inds$train,]
  train_Y <- y_ordinal[inds$train]
  valid_X <- X_sdd[inds$valid,]
  valid_Y <- y_ordinal[inds$valid]
  test_X <- X_sdd[inds$test,]
  test_Y <- y_ordinal[inds$test]
  
  #tuning hyperparameter
  tuning <- data.frame(hyper1=integer(),
                       hyper2 = double(),
                       accuracy = double())
  
  fs <- fisher_score(train_X,train_Y)
  fs_col <- data.frame(colnames(train_X),fs) %>% arrange(desc(fs))
  
  for(l in c(10,50,100,200,300)){
    for(d in (1:10)*0.1){
      selected_features <- GSFSJNE(X=train_X,y=train_Y,S=fs_col[1:l,1],delta=d)
      if(length(selected_features) > 1){
        Train_X <- train_X[,colnames(train_X) %in% selected_features]
        Valid_X <- valid_X[,colnames(valid_X) %in% selected_features]
        nb <- NaiveBayes(Train_X,train_Y)
        pred_y <- predict(nb, Valid_X)$class
        result <- multiclass.roc(valid_Y,as.numeric(pred_y),direction='<',levels=c("normal prostate tissue","normal prostate adjacent to tumor","primary prostate tumor","metastatic prostate tumor"))$auc
        
        tuning <- rbind(tuning,data.frame(hyper1=l,
                                          hyper2 = d,
                                          accuracy = result)) 
      }
    }
  }
  
  hyperparameter <- tuning[which.max(tuning$accuracy),]
  
  #test
  selected_features <- GSFSJNE(X=train_X,y=train_Y,S=fs_col[1:hyperparameter$hyper1,1],delta=hyperparameter$hyper2)
  Train_X <- train_X[,colnames(train_X) %in% selected_features]
  Test_X <- test_X[,colnames(test_X) %in% selected_features]
  nb <- NaiveBayes(Train_X,train_Y)
  pred_y <- predict(nb, Test_X)$class
  result <- c(result,multiclass.roc(test_Y,as.numeric(pred_y),direction='<',levels=c("normal prostate tissue","normal prostate adjacent to tumor","primary prostate tumor","metastatic prostate tumor"))$auc)
  
  # save selected features
  features[[k]] <- selected_features
}

print(c(mean(result),sd(result)))
sapply(features,length)
