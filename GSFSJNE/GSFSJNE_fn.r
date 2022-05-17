joint_neighborhood_entropy <- function(X_jne, y_jne, red, gene_a, delta){
  B = c(red,gene_a)
  U = dim(X_jne)[1]
  x_jne <- matrix(X_jne[,colnames(X_jne) %in% B],nrow=U)
  
  C= levels(y_jne)
  n <- rep(0, length(C))
  nn = 1
  for(c in C){
    n[nn] <- sum(c == y_jne)
    nn <- nn + 1
  }
  
  xiD <- rep(0,length(y_jne))
  nn = 1
  for(c in C){
    xiD[y_jne == c] <- n[nn]
    nn <- nn + 1
  }
  
  d <- as.matrix(dist(x_jne))
  
  jne <- 0
  for(i in 1:dim(X_jne)[1]){
    jne <- jne + log((sum(d[i,] <= delta & y_jne == y_jne[i]))^2/(U*xiD[i]))
  }
  
  return(-jne/U)
}

GSFSJNE <- function(X,y,S,delta){
  
  red <- c()
  temp <- c()
  
  sig <- 100
  while(sig > 0){
    
    red <- temp
    
    ini_sig <- rep(0,length(S)-length(red))
    is <- 1
    for(a in S[!(S %in% red)]){
      ini_sig[is] <- joint_neighborhood_entropy(X,y,red,a,delta)
      is <- is+1
    }
    
    significane_measure <- data.frame(S[!(S %in% red)],ini_sig) %>% arrange(desc(ini_sig))
    temp <- c(red,significane_measure[1,1])
    
    sig <- joint_neighborhood_entropy(X,y,S,c(),delta) - joint_neighborhood_entropy(X,y,temp,c(),delta) 
  }
  
  return(red)
}