
# Compute dyadic.reciprocity #

SRM.dyadic.reciprocity <- function(X){
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
  dyadic.reciprocity <- array(dim=c(T,G),dimnames=list(times,groups),0.)
  for (k in 1:T){
     for (l in 1:G){
        dyadic.reciprocity[k,l] <- (MSbetween[k,l] - MSwithin[k,l])/(MSbetween[k,l] + MSwithin[k,l])}}
  return(dyadic.reciprocity)
}
