
# Function to compute relative variance due to actor, partner and relationship effects i block designs #

block.SRM.relative.variances <- function(X){ 
  if (class(X) != "srmBlockMatrix")
    stop("Data input need to be a srmBlockMatrix object")
  N <- X$individuals
  T <- X$ntimes
  G <- X$ngroups
  G1 <- X$individuals1
  G2 <- X$individuals2
  names1 <- X$names1
  names2 <- X$names2
  groups <- X$groups
  subgroups <- X$subgroups
  times <- X$times
  variances <- block.SRM.variances(X)
  variance.actor <- variances$actor.variance
  variance.partner <- variances$partner.variance
  variance.relationship <- variances$relationship.variance
  prop.variance.actor <- array(dim=c(2,T,G),dimnames=list(subgroups,times,groups),0.)
  prop.variance.partner <- array(dim=c(2,T,G),dimnames=list(subgroups,times,groups),0.)
  prop.variance.relationship <- array(dim=c(2,T,G),dimnames=list(subgroups,times,groups),0.)
  total.variance <- variance.actor + variance.partner + variance.relationship
  prop.variance.actor <- variance.actor/total.variance
  prop.variance.partner <- variance.partner/total.variance
  prop.variance.relationship <- variance.relationship/total.variance
  list(relative.actor.variance=prop.variance.actor,relative.partner.variance=prop.variance.partner,
       relative.relationship.variance=prop.variance.relationship)
}