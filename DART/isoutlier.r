isoutlier <- function(data){
  (data > mean(data)+3*sd(data))|(data < mean(data)-3*sd(data))
}

setxor <- function(x, y) {
  setdiff(union(x, y), intersect(x, y))
}
