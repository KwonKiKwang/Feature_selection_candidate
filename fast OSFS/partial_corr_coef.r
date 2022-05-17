#partial_corr_coef
partial_corr_coef <- function(S,i,j,Y){
  #S = C
  #i = X1
  #j = Y1
  #Y = S1
  
  # S is the covariance (or correlation) matrix for X, Y, Z
  # where X=[i j], Y is conditioned on, and Z is marginalized out.
  # Let S2 = Cov[X | Y] be the partial covariance matrix.
  
  X = c(i,j)
  i2 = 1
  j2 = 2
  S2 = as.matrix(S[X,X],nrow=length(X)) - as.matrix(S[X,Y],nrow=length(X))%*%solve(S[Y,Y])%*%S[Y,X]
  c = S2[i2,j2]
  if(S2[i2,i2]*S2[j2,j2]==0){
    r = Inf
  } else {
    r = c/ sqrt(S2[i2,i2]*S2[j2,j2]) 
  }
  
  return(r)
}
