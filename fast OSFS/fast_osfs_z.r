#fast_osfs_z input Á¤¸®
#for continuous data

fast_osfs_z <- function(data1,class_index,alpha){
  n = nrow(data1)
  p = ncol(data1)
  ns = apply(data1,2,max)
  selected_features = c()
  selected_features1 = c()
  b= c()
  
  for(i in 1:(p-1)){ #p-1
    
    print(i)
    #for very sparse data
    n1 = sum(data1[,i])
    if(n1 == 0){
      next
    }
    
    stop = 0
    CI = 1
    
    indep_fisher_z = cor.test(data1[,i],data1[,class_index],alternative = "two.sided",method="pearson")
    if(indep_fisher_z$p.value >= alpha){ #cannot reject the null hypothesis that true correlation is 0.
      next
    } else{
      stop = 1
      CI=0
    }
    
    if(stop){
      
      if(length(selected_features)>0){
        result1 <- computer_dep_2(selected_features,i,class_index,3,alpha,data1)
        CI <- as.integer(result1[1])
        dep <- result1[2]
      }
      
      if(CI==0){
        selected_features = c(selected_features,i)
        p2= length(selected_features)
        selected_features1 = selected_features
        
        for(j in 1:p2){#p2=4
          b = setdiff(selected_features1,selected_features[j])
          
          if(length(b)>0){
            result2 <- optimal_compter_dep_2(b,selected_features[j],class_index,3,alpha,data1)
            CI <- as.integer(result2[1])
            dep <- result2[2]
            
            if(CI==1|is.na(dep)){
              selected_features1=b
            }
          }
        }
      }
    }
    selected_features=selected_features1
  }
  return(selected_features)
}
