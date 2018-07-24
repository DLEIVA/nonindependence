
# Function to compute reliability for the estimates of the actor and partner variances #

SRM.reliability <- function(X){
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
  covariance.relationship <- variances[[7]]
  actor.reliability <- array(dim=c(T,G),dimnames=list(times,groups),0.)
  actor.reliability <- variance.actor/(variance.actor + variance.relationship/(N-1)-
                       covariance.relationship/(N-1)**2)
  partner.reliability <- array(dim=c(T,G),dimnames=list(times,groups),0.)
  partner.reliability <- variance.partner/(variance.partner + variance.relationship/(N-1)-
                         covariance.relationship/(N-1)**2)
  list(actor.reliability=actor.reliability,partner.reliability=partner.reliability)
}