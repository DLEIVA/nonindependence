
# Compute dyadic.reciprocity in block designs #

block.SRM.dyadic.reciprocity <- function(X,symmetric=FALSE){
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
  effects <- block.SRM.effects(X)
  relationship <- effects$relationship.effects 
  variances <- block.SRM.variances(X)
  variance.relationship <- variances$relationship.variance
  # Estimate covariance for actor-partner effect #
   dyads1 <- array(dim=c(G1*G2,T,G),0.)
   dyads2 <- array(dim=c(G1*G2,T,G),0.)
   for (k in 1:T){
      for (l in 1:G){
         count1 <- 1
         count2 <- 1
         for (i in 1:N){
            for (j in 1:N){
               if ((i <= G1) && (j > G1))
                {dyads1[count1,k,l] <- relationship[i,j,k,l]
                 count1 <- count1 +1}
               if ( (i > G1) && (j <= G1))
                {dyads2[count2,k,l] <- relationship[i,j,k,l]
                 count2 <- count2 +1}}}}}
  dyadic.reciprocity <- array(dim=c(T,G),dimnames=list(times,groups),0.)
  if (symmetric == FALSE) {
  for (k in 1:T){
     for (l in 1:G){
        dyadic.reciprocity[k,l] <- cor(dyads1[,k,l],dyads2[,k,l])}}}
  else {
    for (k in 1:T){
       for (l in 1:G){
          dyadic.reciprocity[k,l] <- cor(dyads1[,k,l],dyads2[,k,l])*(sqrt(variance.relationship[1,k,l]*
          variance.relationship[2,k,l])/((variance.relationship[1,k,l]+variance.relationship[2,k,l])/2))}}}
  return(dyadic.reciprocity)
}