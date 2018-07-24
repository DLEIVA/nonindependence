
# Function to estimate generalized reciprocity #

SRM.generalized.reciprocity <- function(X){
  if (class(X) != "srmRRMatrix")
    stop("Data input need to be a srmRRMatrix object")
  N <- X$individuals
  T <- X$ntimes
  G <- X$ngroups
  names <- X$names
  times <- X$times
  groups <- X$groups
  effects <- SRM.effects(X)
  actor <- effects$actor.effects
  partner <- effects$partner.effects
  variances <- SRM.variances(X)
  variance.actor <- variances[[3]]
  variance.partner <- variances[[4]]
  MSbetween <- variances[[1]]
  MSwithin <- variances[[2]]
  generalized.reciprocity <- array(dim=c(T,G),dimnames=list(times,groups),0.)
  for (k in 1:T){
     for (l in 1:G){
        if (is.na(variance.actor[k,l]) | is.na(variance.partner[k,l])) generalized.reciprocity[k,l] <- NA
          else{
            
            generalized.reciprocity[k,l] <- (sum(actor[,k,l]*partner[,k,l])/(N-1) - MSbetween[k,l]/(2*(N-2)) +
                                            MSwithin[k,l]/(2*N))/sqrt(variance.actor[k,l]*variance.partner[k,l])
            if ((generalized.reciprocity[k,l] > 1.0) | (generalized.reciprocity[k,l] < -1.0))
            generalized.reciprocity[k,l] <- sign(generalized.reciprocity[k,l])*1.0}}}
  return(generalized.reciprocity)
}