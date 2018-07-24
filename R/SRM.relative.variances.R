
# Function to compute relative variance due to actor, partner and relationship effects #

SRM.relative.variances <- function(X){
  if (class(X) != "srmRRMatrix")
    stop("Data input need to be a srmRRMatrix object")
  N <- X$individuals
  T <- X$ntimes
  G <- X$ngroups
  names <- X$names
  times <- X$times
  groups <- X$groups
  variances <- SRM.variances(X)
  variance.actor <- variances[[3]]
  variance.partner <- variances[[4]]
  variance.relationship <- variances[[5]]
  prop.variance.actor <- array(dim=c(T,G),dimnames=list(times,groups),0.)
  prop.variance.partner <- array(dim=c(T,G),dimnames=list(times,groups),0.)
  prop.variance.relationship <- array(dim=c(T,G),dimnames=list(times,groups),0.)
  total.variance <- sum(unlist(variances[3:5]),na.rm=TRUE)
  prop.variance.actor <- variance.actor/total.variance
  prop.variance.partner <- variance.partner/total.variance
  prop.variance.relationship <- variance.relationship/total.variance
  res<-list(actor.variance=prop.variance.actor,partner.variance=prop.variance.partner,
       relationship.variance=prop.variance.relationship)
  class(res) <- "srmRRrelVars"
  res
}