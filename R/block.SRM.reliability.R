
# Function to compute reliability for the estimates of the actor and partner variances in block designs #

block.SRM.reliability <- function(X){
  if (class(X) != "srmBlockMatrix")
    stop("Data input need to be a srmBlockMatrix object")
  N <- X$individuals
  T <- X$ntimes
  G <- X$ngroups
  G1 <- X$individuals1
  G2 <- X$individuals2
  names1 <- X$names1
  names2 <- X$names2
  subgroups <- X$subgroups
  groups <- X$sgroups
  times <- X$times  
  variances <- block.SRM.variances(X)
  variance.actor <- variances$actor.variance
  variance.partner <- variances$partner.variance
  variance.relationship <- variances$relationship.variance
  actor.reliability <- array(dim=c(2,T,G),dimnames=list(subgroups,times,groups),0)
  partner.reliability <- array(dim=c(2,T,G),dimnames=list(subgroups,times,groups),0.)
  for (k in 1:T){
     for (l in 1:G){
        if ( (variance.actor[1,k,l] != 0) | (variance.relationship[1,k,l] != 0) )
          actor.reliability[1,k,l] <- variance.actor[1,k,l]/(variance.actor[1,k,l] + 
          variance.relationship[1,k,l]/G1)
        if ( (variance.actor[2,k,l] != 0) | (variance.relationship[2,k,l] != 0) )
          actor.reliability[2,k,l] <- variance.actor[2,k,l]/(variance.actor[2,k,l] + 
          variance.relationship[2,k,l]/G2)
        if ( (variance.partner[1,k,l] != 0) | (variance.relationship[1,k,l] != 0) )
          partner.reliability[1,k,l] <- variance.partner[1,k,l]/(variance.partner[1,k,l] + 
          variance.relationship[1,k,l]/G1)
        if ( (variance.partner[2,k,l] != 0) | (variance.relationship[2,k,l] != 0) )
          partner.reliability[2,k,l] <- variance.partner[2,k,l]/(variance.partner[2,k,l] + 
          variance.relationship[1,k,l]/G2)}}
  list(actor.reliability=actor.reliability,partner.reliability=partner.reliability)
}