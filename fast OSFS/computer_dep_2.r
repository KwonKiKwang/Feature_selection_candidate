#input 

computer_dep_2 <- function(bcf,var,target,max_k,alpha,data){
  #content
  dep1=0
  x=0
  n_bcf=length(bcf)
  code=bcf
  N = nrow(data)
  max_cond_size=max_k
  CI=0
  p_value=1
  if(max_cond_size>n_bcf){
    max_cond_size=n_bcf
  }
  cond = c()
  cond_size = 1
  while(cond_size <= max_cond_size){
    
    cond_index = rep(0,cond_size)
    for(i in 1:cond_size){
      cond_index[i] = i
    }
    stop = 0
    
    while(stop==0){
      cond = c()
      
      for(i in 1:cond_size){
        cond = c(cond,code[cond_index[i]])
      }
      
      #Only consider continuout data
      result1 <- my_cond_indep_fisher_z(data,var,target,cond,N,alpha)
      CI <- result1[1]
      r <- result1[2]
      p_value <- result1[3]
      x=r
      
      if(CI==1|is.na(x)){
        stop = 1
        cond_size = max_cond_size+1
      }
      
      if(stop==0){
        result2 <- next_cond_index(n_bcf,cond_size,cond_index)
        cond_index <- result2[1]
        stop <- result2[2]
      }
    }
    cond_size = cond_size+1
  }
  dep1=x
  
  return(c(CI,dep1,p_value))
}


# next_cond_index
#input
#n_bcf = n_bcf
#cond_size = cond_size
#cond_index1 = cond_index

next_cond_index <- function(n_bcf,cond_size,cond_index1){
  
  stop=1
  i=cond_size
  
  while(i>=1){
    if(cond_index1[i]<(n_bcf+i-cond_size)){
      cond_index1[i]=cond_index1[i]+1
      if(i < cond_size){
        for(j in (i+1):cond_size){
          cond_index1[j] = cond_index1[j-1]+1
        }
      }
      stop=0
      i=-1
    }
    i=i-1
  }
  cond_index= cond_index1
  
  return(c(cond_index,stop))
}
