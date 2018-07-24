
# Function for obtaining SRM effects #

SRM.effects <- function (X){
  if (class(X) != "srmRRMatrix")
  stop("Data input need to be a srmRRMatrix object")
  N <- X$individuals
  T <- X$ntimes
  G <- X$ngroups
  # Compute row means #
  row.means <- array(dim=c(N,G),0.)
  for (l in 1:G){
     row.means[,l] <- rowSums(X$data[,,,l])/(G*(N-1));}
  time.rowmeans <- array(dim=c(N,T,G),0.)
  for (k in 1:T){
     for (l in 1:G){
        time.rowmeans[,k,l] <- rowSums(X$data[,,k,l])/(N-1)}}
  # Compute column means #
  column.means <- array(dim=c(N,G),0.)
  if (T <= 1) {
    for (l in 1:G){
       column.means[,l] <- colSums(X$data[,,,l])/(T*(N-1))}}
  if (T > 1) {
    for (l in 1:G){
       column.means[,l] <- rowSums(colSums(X$data[,,,l]))/(T*(N-1))}}
  time.colmeans <- array(dim=c(N,T,G),0.)
  for (k in 1:T){
     for (l in 1:G){
        time.colmeans[,k,l] <- colSums(X$data[,,k,l])/(N-1)}}
  # Compute grand mean #
  grand.mean <- array(dim=c(G,1),0.)
  for (l in 1:G){
     grand.mean[l] <- sum(X$data[,,,l])/(T*(N*(N-1)))}
  time.grandmean <- array(dim=c(T,G),0.)
  for (k in 1:T){
     for (l in 1:G){
        time.grandmean[k,l] <- sum(X$data[,,k,l])/((N*(N-1)))}}
  # Estimate actor effects for all individuals #
  actor <- array(dim=c(N,T,G),dimnames=list(X$names,X$times,X$groups),0.)
  for (k in 1:T){
     for (l in 1:G){
        for (i in 1:N){
           actor[i,k,l] <- (N-1)*((N-1)*time.rowmeans[i,k,l]+time.colmeans[i,k,l]-
                           N*time.grandmean[k,l])/(N*(N-2))}}}
  # Estimate partner effects for all individuals #
  partner <- array(dim=c(N,T,G),dimnames=list(X$names,X$times,X$groups),0.)
  for (k in 1:T){
     for (l in 1:G){
        for (i in 1:N){
           partner[i,k,l] <- (N-1)*((N-1)*time.colmeans[i,k,l]+time.rowmeans[i,k,l]-
                             N*time.grandmean[k,l])/(N*(N-2))}}}
  # Estimate relationships effects for all individuals #
  relationship <- array(dim=c(N,N,T,G),
                  dimnames=list(X$names,X$names,X$times,X$groups),0.)
  for (i in 1:N){
     for (j in 1:N){
        for (k in 1:T){
           for (l in 1:G){
              if (i != j)
                {relationship[i,j,k,l] <- X$data[i,j,k,l]-actor[i,k,l]-partner[j,k,l]-time.grandmean[k,l]}}}}}
  res<-list(actor.effects=actor,partner.effects=partner,relationship.effects=relationship)
  class(res)="srmRREffects"
  res
}