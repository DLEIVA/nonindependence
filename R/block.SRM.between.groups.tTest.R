
# Testing SRM effects in block designs by means of a between-group t Test #

block.SRM.between.groups.tTest <-function(X){
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
  # Are there more than 1 group? #
  if (G < 2)
    return("Error: Number of groups analyzed should be greater than 1.")
  variances <- block.SRM.variances(X)
  variance.actor <- variances$actor.variance
  variance.partner <- variances$partner.variance
  variance.relationship <- variances$relationship.variance
  actorpartner.covariance <- variances$actorpartner.covariance
  dyadic.covariance <- variances$dyadic.covariance
  source <- c("Actor Var","Partner Var","Relationship Var","Actor-Partner Cov","Dyadic Cov")
  between.group.t.statistic <- array(dim=c(2,T,5),dimnames=list(subgroups,times,source))
  between.group.p.value <- array(dim=c(2,T,5),dimnames=list(subgroups,times,source))
  actor.variance.mean <- array(dim=c(2,T))
  partner.variance.mean <- array(dim=c(2,T))
  relationship.variance.mean <- array(dim=c(2,T)) 
  actorpartner.covariance.mean <- array(dim=c(T))
  dyadic.covariance.mean <- array(dim=c(2,T))
  actor.variance.var <- array(dim=c(2,T))
  partner.variance.var <- array(dim=c(2,T))
  relationship.variance.var <- array(dim=c(2,T)) 
  actorpartner.covariance.var <- array(dim=c(T))
  dyadic.covariance.var <- array(dim=c(2,T))
  for (k in 1:T){
  actor.variance.mean[,k] <- apply(variance.actor[,k,],1,mean)
  partner.variance.mean[,k] <-  apply(variance.partner[,k,],1,mean)
  relationship.variance.mean[,k] <-  apply(variance.relationship[,k,],1,mean)
  actorpartner.covariance.mean[k] <-  mean(actorpartner.covariance[k,])
  dyadic.covariance.mean[,k] <-  apply(dyadic.covariance[,k,],1,mean)
  actor.variance.var[,k] <- apply(variance.actor[,k,],1,var)
  partner.variance.var[,k] <-  apply(variance.partner[,k,],1,var)
  relationship.variance.var[,k] <-  apply(variance.relationship[,k,],1,var)
  actorpartner.covariance.var[k] <-  var(actorpartner.covariance[k,])
  dyadic.covariance.var[,k] <-  apply(dyadic.covariance[,k,],1,var)}
  for (m in 1:length(subgroups)){
     for (k in 1:T){
    between.group.t.statistic[m,k,1] <- actor.variance.mean[m,k]/sqrt(actor.variance.var[m,k]/G)
  	between.group.t.statistic[m,k,2] <- partner.variance.mean[m,k]/sqrt(partner.variance.var[m,k]/G)
  	between.group.t.statistic[m,k,3] <- relationship.variance.mean[m,k]/sqrt(relationship.variance.var[m,k]/G)
  	between.group.t.statistic[m,k,4] <- dyadic.covariance.mean[m,k]/sqrt(dyadic.covariance.var[m,k]/G)
  	between.group.t.statistic[m,k,5] <- actorpartner.covariance.mean[k]/sqrt(actorpartner.covariance.var[k]/G)}}
  df <- G-1
  for (m in 1:(dim(between.group.t.statistic)[1])){ 
  for (k in 1:T){ 
     for (n in 1:(dim(between.group.t.statistic)[3])){
        if ( !is.na(between.group.t.statistic[m,k,n]) ){
        if ( sign(between.group.t.statistic[m,k,n]) > 0 )
          between.group.p.value[m,k,n] <- pt(between.group.t.statistic[m,k,n],df,lower.tail=FALSE)
        else between.group.p.value[m,k,n] <- pt(between.group.t.statistic[m,k,n],df)}}}}
  list(t.statistic=between.group.t.statistic,df=df,two.tailed.p.value=2*between.group.p.value)
}