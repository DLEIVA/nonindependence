
# Function to estimate SRM variances #

SRM.variances <- function(X){
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
  relationship <- effects$relationship.effects
  # Estimate mean squares for actor, partner and cross products actor-partner #
  MSactor <- array(dim=c(T,G),0.)
  for (k in 1:T){
     for (l in 1:G){
        MSactor[k,l] <- sum(actor[,k,l]**2)/(N-1)}}
  MSpartner <- array(dim=c(T,G),0.)
  for (k in 1:T){
     for (l in 1:G){
        MSpartner[k,l] <- sum(partner[,k,l]**2)/(N-1)}}
  MCP <- array(dim=c(T,G),0.)
  for (k in 1:T){
     for (l in 1:G){
        MCP[k,l] <- sum(actor[,k,l]*partner[,k,l])/(N-1)}}
  # Estimate mean squares between and within dyads #
  sum.relationship <- array(dim=c(N,N,T,G),0.)
  for (k in 1:T){
     for (l in 1:G){
        sum.relationship[,,k,l] <- relationship[,,k,l] + t(relationship[,,k,l])}}
  MSbetween <- array(dim=c(T,G),dimnames=list(times,groups),0.)
  for (k in 1:T){
     for (l in 1:G){
        MSbetween[k,l] <- (sum(sum.relationship[,,k,l]**2)/2)/((N-1)*(N-2)-2)}}
  subs.relationship <- array(dim=c(N,N,T,G),0.)
  for (k in 1:T){
     for (l in 1:G){
        subs.relationship[,,k,l] <- relationship[,,k,l] - t(relationship[,,k,l]);}}
  MSwithin <- array(dim=c(T,G),dimnames=list(times,groups),0.)
  for (k in 1:T){
     for (l in 1:G){
        MSwithin[k,l] <- (sum(subs.relationship[,,k,l]**2)/2)/((N-1)*(N-2))}}
  # Estimate variance for actor effect #
  variance.actor <- array(dim=c(T,G),dimnames=list(times,groups),0.)
  for (k in 1:T){
     for (l in 1:G){
        variance.actor[k,l] <- MSactor[k,l] - 0.5*((MSbetween[k,l]/(N-2))+(MSwithin[k,l]/N))
        if (variance.actor[k,l]<0) variance.actor[k,l]<- NA}}
  # Estimate variance for partner effect #
  variance.partner <- array(dim=c(T,G),dimnames=list(times,groups),0.)
  for (k in 1:T){
     for (l in 1:G){
        variance.partner[k,l] <- MSpartner[k,l] - 0.5*((MSbetween[k,l]/(N-2))+(MSwithin[k,l]/N))
        if (variance.partner[k,l]<0) variance.partner[k,l]<- NA}}
  # Estimate covariance for actor-partner effect #
  covariance.actorpartner <- array(dim=c(T,G),dimnames=list(times,groups),0.)
  for (k in 1:T){
     for (l in 1:G){
        covariance.actorpartner[k,l] <- MCP[k,l] - 0.5*((MSbetween[k,l]/(N-2))-(MSwithin[k,l]/N))}}
  # Estimate variance for relationship effect #
  variance.relationship <- array(dim=c(T,G),dimnames=list(times,groups),0.)
  for (k in 1:T){
     for (l in 1:G){
        variance.relationship[k,l] <- 0.5*(MSbetween[k,l]+MSwithin[k,l])
        if (variance.relationship[k,l]<0) variance.relationship[k,l]<- NA}}
  # Estimate covariance for relationship effect #
  covariance.relationship <- array(dim=c(T,G),dimnames=list(times,groups),0.)
  for (k in 1:T){
     for (l in 1:G){
        covariance.relationship[k,l] <- 0.5*(MSbetween[k,l]-MSwithin[k,l])}}
  res <- list(MSbetween=MSbetween,MSwithin=MSwithin,actor.variance=variance.actor,partner.variance=variance.partner,
       relationship.variance=variance.relationship,actorpartner.covariance=covariance.actorpartner,
       dyadic.covariance=covariance.relationship)

  class(res) <- "srmRRVars"
  res
}