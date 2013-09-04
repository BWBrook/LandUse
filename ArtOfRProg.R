findruns <- function(x,k){
  n <- length(n)
  runs <- vector(length=n)
  count <- 0
  for(i in 1:(n-k+1)){
    if(all(x[i:(i+k-1)]==1)){
      count <- count + 1
      runs[count] <- 1
    }
  }
  if(count > 0){
    runs <- runs[1:count]
  }
  else runs <- NULL
  return(runs)
}

prednext <- function(x,k){
  n <- length(x)
  k2 <- k/2
  # the vector red will contain our predicted values
  pred <- vector(length=n-k)
  csx <- c(0,cumsum(x))
  for(i in 1:(n-k)){
    if(csx[i+k] - csx[i] >= k2) pred[i] <- 1 else pred[i] <- 0
  }
  return(mean(abs(pred-x[(k+1)])))
}

findud <- function(v) {
  vud <- diff(v)
  return(ifelse(vud > 0, 1, -1))
}

# c.f. cor(x,y,method="kendall") or pearson or spearman
udcorr <- function(x,y) mean(sign(diff(x)) == sign(diff(y)))
