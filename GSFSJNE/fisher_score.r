fisher_score <- function(X_fs,Y_fs){
  fs_s <- rep(0,dim(X_fs)[2])
  
  C= levels(Y_fs)
  n <- rep(0, length(C))
  nn = 1
  for(c in C){
    n[nn] <- sum(c == Y_fs)
    nn <- nn + 1
  }
  
  for(i in 1:dim(X_fs)[2]){
    mu = mean(X_fs[,i])
    
    nn <- 1
    bcs <- 0
    for(c in C){
      bcs <- bcs + n[nn]*(mean(X_fs[c == Y_fs,i])-mu)^2
      nn <- nn + 1
    }
    
    nn <- 1
    wcs <- 0 
    for(c in C){
      wcs <- wcs + n[nn]*(sd(X_fs[c == Y_fs,i]))^2
      nn <- nn + 1
    }
    
    fs_s[i] = bcs/wcs
  }
  
  return(fs_s)
}