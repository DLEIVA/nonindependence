
# Testing dyadic covariance by means of a F ratio #

SRM.dyadic.covariance.Ftest <- function(X){
  if (class(X) != "srmRRMatrix")
    stop("Data input need to be a srmRRMatrix object")
  N <- X$individuals
  T <- X$ntimes
  G <- X$ngroups
  names <- X$names
  times <- X$times
  groups <- X$groups
  variances <- SRM.variances(X)
  MSbetween <- variances[[1]]
  MSwithin <- variances[[2]]
  F.value <- array(dim=c(T,G),dimnames=list(times,groups))
  for (k in 1:T){
     for (l in 1:G){
        F.value[k,l] <- MSbetween[k,l]/MSwithin[k,l]}}
  df1 <- (N-1)*(N-2)/2-1
  df2 <- (N-1)*(N-2)/2
  F.p.value <- array(dim=c(T,G),dimnames=list(times,groups))
  for (k in 1:T){
     for (l in 1:G){
        if (sign(F.value[k,l])>0)
          F.p.value[k,l] <- pf(F.value[k,l],df1,df2,lower.tail=FALSE)
          else F.p.value[k,l] <- pf(F.value[k,l],df1,df2)}}
  res=list(MSb=MSbetween,MSw=MSwithin,F.value=F.value,df1=df1,df2=df2,p.value=F.p.value)
  class(res) <- "srmCovFTest"
  res
}
