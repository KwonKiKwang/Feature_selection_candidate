CFS <- function(m,X,CFS_y){
  S = c()
  
  S <- colnames(X)[which.max(cor(X,CFS_y))]
  X_cand = X[,-which.max(cor(X,CFS_y))]
  
  while(length(S)>0 & length(S) < m){
    k = length(S)
    merits = k*cor(X_cand,CFS_y)/sqrt(k+k*(k-1)*apply(abs(cor(X_cand)),2,mean))
    S <- c(S,colnames(X_cand)[which.max(abs(merits))])
    X_cand = X_cand[,-which.max(abs(merits))]
    
    if(length(S) %% 10 == 0){
      print(length(S))
    }
  }
  
  
  return(S)
}

