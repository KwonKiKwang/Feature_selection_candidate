CFS_0702 <- function(X,CFS_y){
  S = c()
  
  #first iteration
  S <- colnames(X)[which.max(cor(X,CFS_y))]
  X_cand = X[,-which.max(cor(X,CFS_y))]
  
  #initial value
  i=0
  old_merits =0
  
  while(i < 5){
    k = length(S)
    merits = k*cor(X_cand,CFS_y)/sqrt(k+k*(k-1)*apply(abs(cor(X_cand)),2,mean))
    S <- c(S,colnames(X_cand)[which.max(abs(merits))])
    X_cand = X_cand[,-which.max(abs(merits))]
    
    #Finished when no improvement for five consecutive iterations
    if(max(abs(merits)) <= old_merits){ 
      i = i+1
    }
    old_merits = max(abs(merits))
    
    
    if(length(S) %% 10 == 0){
      print(length(S))
    }
  }
  
  
  return(S)
}

