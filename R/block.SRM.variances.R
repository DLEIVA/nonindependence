
# Function to estimate SRM variances in block designs #

block.SRM.variances <- function(X){ 
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
  effects <- block.SRM.effects(X)
  actor <- effects$actor.effects
  partner <- effects$partner.effects
  relationship <- effects$relationship.effects 
  # Estimate variance for relationship effect #
  variance.relationship <- array(dim=c(2,T,G),dimnames=list(subgroups,times,groups),0.)
  for (i in 1:N){
     for (j in 1:N){
        for (k in 1:T){
           for (l in 1:G){
              if ((i <= G1) && (j > G1))
                {variance.relationship[1,k,l] <- variance.relationship[1,k,l]+relationship[i,j,k,l]**2}
              if ( (i > G1) && (j <= G1)) 
                {variance.relationship[2,k,l] <- variance.relationship[2,k,l]+relationship[i,j,k,l]**2}}}}}
  variance.relationship <- variance.relationship/((G1-1)*(G2-1))
  # Estimate variance for actor effect #
  variance.actor <- array(dim=c(2,T,G),dimnames=list(subgroups,times,groups),0.)
  for (i in 1:N){
        for (k in 1:T){
           for (l in 1:G){
              if ( i <= G1 )
                {variance.actor[1,k,l] <- variance.actor[1,k,l]+actor[i,k,l]**2}
              if ( i > G1 ) 
                {variance.actor[2,k,l] <- variance.actor[2,k,l]+actor[i,k,l]**2}}}}
  variance.actor <- (variance.actor/(G1-1))-(variance.relationship/G2)
  for (k in 1:T){
     for (l in 1:G){
        for (m in 1:length(subgroups)){
           if ( sign(variance.actor[m,k,l]) < 0. ) variance.actor[m,k,l] <- 0.}}} 
  # Estimate variance for partner effect #
  variance.partner <- array(dim=c(2,T,G),dimnames=list(subgroups,times,groups),0.)
  for (j in 1:N){
        for (k in 1:T){
           for (l in 1:G){
              if ( j <= G2 )
                {variance.partner[1,k,l] <- variance.partner[1,k,l]+partner[j,k,l]**2}
              if ( j > G1 ) 
                {variance.partner[2,k,l] <- variance.partner[2,k,l]+partner[j,k,l]**2}}}}
  variance.partner[1,,] <- (variance.partner[1,,]/(G2-1))-(variance.relationship[2]/G2)
  variance.partner[2,,] <- (variance.partner[2,,]/(G2-1))-(variance.relationship[1]/G2)
  for (k in 1:T){
     for (l in 1:G){
        for (m in 1:length(subgroups)){
           if ( sign(variance.partner[m,k,l]) < 0. ) variance.partner[m,k,l] <- 0.}}} 
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
  covariance.actorpartner <- array(dim=c(T,G),dimnames=list(times,groups),0.)
  for (k in 1:T){
     for (l in 1:G){
        covariance.actorpartner[k,l] <- cov(dyads1[,k,l],dyads2[,k,l])}}
  # Estimate covariance for relationship effect #
  covariance.relationship <- array(dim=c(2,T,G),dimnames=list(subgroups,times,groups),0.)
  crosspod.relationship <- array(dim=c(N,N,T,G),0.)
  for (k in 1:T){
     for (l in 1:G){
        crosspod.relationship[,,k,l] <- relationship[,,k,l] * t(relationship[,,k,l])}}
  for (i in 1:N){
        for (k in 1:T){
           for (l in 1:G){
              if ( i <= G1 )
                {covariance.relationship[1,k,l] <- covariance.relationship[1,k,l]+actor[i,k,l]*partner[i,k,l]}
              if ( i > (N-G2) ) 
                {covariance.relationship[2,k,l] <- covariance.relationship[2,k,l]+actor[i,k,l]*partner[i,k,l]}}}}
  for (k in 1:T){
     for (l in 1:G){
        covariance.relationship[1,k,l] <- (covariance.relationship[1,k,l]/(G2-1))-(sum(crosspod.relationship[,,k,l])/
                                          2/((G2-1)*(G1-1)))/G2
        covariance.relationship[2,k,l] <- (covariance.relationship[2,k,l]/(G1-1))-(sum(crosspod.relationship[,,k,l])/
                                          2/((G2-1)*(G1-1)))/G2}}
  list(actor.variance=variance.actor,partner.variance=variance.partner,relationship.variance=variance.relationship,
       actorpartner.covariance=covariance.actorpartner,dyadic.covariance=covariance.relationship)
}
