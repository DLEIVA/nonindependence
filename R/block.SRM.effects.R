
# Function for obtaining SRM effects in Block designs #

block.SRM.effects <- function (X){
  if (class(X) != "srmBlockMatrix")
    stop("Data input need to be a srmBlockMatrix object")
  N <- X$individuals
  T <- X$ntimes
  G <- X$ngroups
  G1 <- X$individuals1
  G2 <- X$individuals2
  names1 <- X$names1
  names2 <- X$names2
  groups<- X$groups
  times <- X$times
  # Compute row means #
  time.rowmeans <- array(dim=c(N,T,G),0.)
  for (k in 1:T){
     for (l in 1:G){
        for (i in 1:N){
           if (i <= G1) time.rowmeans[i,k,l] <- sum(X$data[i,,k,l])/G2
           else time.rowmeans[i,k,l] <- sum(X$data[i,,k,l])/G1}}}
  # Compute column means #
  time.colmeans <- array(dim=c(N,T,G),0.)
  for (k in 1:T){
     for (l in 1:G){
         for (j in 1:N){
            if (j <= G1) time.colmeans[j,k,l] <- sum(X$data[,j,k,l])/G2
              else time.colmeans[j,k,l] <- sum(X$data[,j,k,l])/G1}}}
  # Compute grand mean #
  time.grandmean <- array(dim=c(2,T,G),0.)
  for (i in 1:N){
     for (j in 1:N){
        for (k in 1:T){
           for (l in 1:G){
              if ((i <= G1) && (j > G1)){
                time.grandmean[1,k,l] <- time.grandmean[1,k,l] + X$data[i,j,k,l]}
              if ( (i > G1) && (j <= G1)){
                time.grandmean[2,k,l] <- time.grandmean[2,k,l] + X$data[i,j,k,l]}}}}}
  time.grandmean <- time.grandmean/(G1*G2)
  # Estimate actor effects for all individuals #
  actor <- array(dim=c(N,T,G),dimnames=list(c(names1,names2),times,groups),0.)
  for (k in 1:T){
     for (l in 1:G){
        for (i in 1:N){
           if (i <= G1){
             actor[i,k,l] <- time.rowmeans[i,k,l]-time.grandmean[1,k,l]}
           else actor[i,k,l] <- time.rowmeans[i,k,l]-time.grandmean[2,k,l]}}}
  # Estimate partner effects for all individuals #
  partner <- array(dim=c(N,T,G),dimnames=list(c(names1,names2),times,groups),0.)
  for (j in 1:N){
     for (k in 1:T){
        for (l in 1:G){
           if (j > G1){
             partner[j,k,l] <- time.colmeans[j,k,l]-time.grandmean[1,k,l]}
           if (j <= G1){
             partner[j,k,l] <- time.colmeans[j,k,l]-time.grandmean[2,k,l]}}}}
  # Estimate relationships effects for all individuals #
  relationship <- array(dim=c(N,N,T,G),
                  dimnames=list(c(names1,names2),c(names1,names2),times,groups),0.)
  for (i in 1:N){
     for (j in 1:N){
        for (k in 1:T){
           for (l in 1:G){
              if ((i <= G1) && (j > G1))
                {relationship[i,j,k,l] <- X$data[i,j,k,l]-actor[i,k,l]-partner[j,k,l]-time.grandmean[1,k,l]}
              if ( (i > G1) && (j <= G1))
                {relationship[i,j,k,l] <- X$data[i,j,k,l]-actor[i,k,l]-partner[j,k,l]-time.grandmean[2,k,l]}}}}}
  list(actor.effects=actor,partner.effects=partner,relationship.effects=relationship)
}