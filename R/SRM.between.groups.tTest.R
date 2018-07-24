
# Testing SRM effects by means of a between-group t Test #

SRM.between.groups.tTest <- function(X){
  if (class(X) != "srmRRMatrix")
    stop("Data input need to be a srmRRMatrix object")
  N <- X$individuals
  T <- X$ntimes
  G <- X$ngroups
  names <- X$names
  times <- X$times
  # Are there more than 1 group? #
  if (G < 2)
    return("Error: Number of groups analyzed should be greater than 1.")

  variances <- SRM.variances(X)
  variance.actor <- variances[[3]]
  variance.partner <- variances[[4]]
  variance.relationship <- variances[[5]]
  actorpartner.covariance <- variances[[6]]
  dyadic.covariance <- variances[[7]]
  source <- c("Actor Var","Partner Var","Relationship Var","Actor-Partner Cov","Dyadic Cov")
  between.group.t.statistic <- array(dim=c(T,5),dimnames=list(times,source))
  between.group.p.value <- array(dim=c(T,5),dimnames=list(times,source))
  actor.variance.mean <- apply(variance.actor,1,mean)
  partner.variance.mean <- apply(variance.partner,1,mean)
  relationship.variance.mean <- apply(variance.relationship,1,mean)
  actorpartner.covariance.mean <- apply(actorpartner.covariance,1,mean)
  dyadic.covariance.mean <- apply(dyadic.covariance,1,mean)
  actor.variance.var <- apply(variance.actor,1,var)
  partner.variance.var <- apply(variance.partner,1,var)
  relationship.variance.var <- apply(variance.relationship,1,var)
  actorpartner.covariance.var <- apply(actorpartner.covariance,1,var)
  dyadic.covariance.var <- apply(dyadic.covariance,1,var)
  for (k in 1:T){
  between.group.t.statistic[k,1] <- actor.variance.mean[k]/sqrt(actor.variance.var[k]/G)
  between.group.t.statistic[k,2] <- partner.variance.mean[k]/sqrt(partner.variance.var[k]/G)
  between.group.t.statistic[k,3] <- relationship.variance.mean[k]/sqrt(relationship.variance.var[k]/G)
  between.group.t.statistic[k,4] <- dyadic.covariance.mean[k]/sqrt(dyadic.covariance.var[k]/G)
  between.group.t.statistic[k,5] <- actorpartner.covariance.mean[k]/sqrt(actorpartner.covariance.var[k]/G)}
  df <- G-1
  for (k in 1:T){ 
     for (m in 1:(dim(between.group.t.statistic)[2])){
        if ( !is.na(between.group.t.statistic[k,m])){
        if (sign(between.group.t.statistic[k,m]) > 0)
          between.group.p.value[k,m] <- pt(between.group.t.statistic[k,m],df,lower.tail=FALSE)
        else between.group.p.value[k,m] <- pt(between.group.t.statistic[k,m],df)}}}
  list(t.statistic=between.group.t.statistic,df=df,two.tailed.p.value=2*between.group.p.value)
}