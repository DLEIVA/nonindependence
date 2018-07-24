
# Function to estimate generalized reciprocity in block designs #

block.SRM.generalized.reciprocity <- function(X){
  if (class(X) != "srmBlockMatrix")
    stop("Data input need to be a srmBlockMatrix object")
  N <- X$individuals
  T <- X$ntimes
  G <- X$ngroups
  G1 <- X$individuals1
  G2 <- X$individuals2
  names1 <- X$names1
  names2 <- X$names2
  subgroups<- X$subgroups
  groups <- X$groups  
  times <- X$times  
  variances <- block.SRM.variances(X)
  variance.actor <- variances$actor.variance
  variance.partner <- variances$partner.variance
  covariance.relationship <- variances$dyadic.covariance
  generalized.reciprocity <- array(dim=c(2,T,G),dimnames=list(subgroups,times,groups),0.)
  for (m in 1:length(subgroups)){ 
     for (k in 1:T){
        for (l in 1:G){
           if ((variance.actor[m,k,l] == 0.) | (variance.partner[m,k,l] == 0.)) generalized.reciprocity[m,k,l] <- 0.
           else 
           generalized.reciprocity[m,k,l] <- covariance.relationship[m,k,l]/sqrt(variance.actor[m,k,l]*variance.partner[m,k,l])}}
    if ((generalized.reciprocity[m,k,l] > 1.0) | (generalized.reciprocity[m,k,l] < -1.0))
    generalized.reciprocity[m,k,l] <- sign(generalized.reciprocity[m,k,l])*1.0}
  return(generalized.reciprocity)
}