#input

my_cond_indep_fisher_z <- function(data,X,Y,S,N,alpha){
  
  C= cov(data[,c(X,Y,S)])
  size_C = ncol(C)
  X1 =1
  Y1 =2
  S1 = c(3:size_C)
  
  r = partial_corr_coef(C,X1,Y1,S1)
  z = 0.5*log((1+r)/(1-r))
  z0 = 0
  W = sqrt(N-length(S1)-3)*(z-z0) # W ~ N(0,1)
  cutoff = qnorm(1-0.5*alpha) # P(|W| <= cutoff) = 0.95
  
  if(abs(W) < cutoff){
    CI = 1
  } else{
    CI = 0
  }
  
  p=pnorm(W)
  r = abs(r)
  
  return(c(CI,r,p))
}


