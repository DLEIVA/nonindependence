
# Testing SRM effects by means of Jackknife tests #

SRM.jackknife <- function(X){
  if (class(X) != "srmRRMatrix")
    stop("Data input need to be a srmRRMatrix object")
  N <- X$individuals
  T <- X$ntimes
  G <- X$ngroups
  names <- X$names
  times <- X$times
  groups <- X$groups
  # Is group(s) size greater than or equal to 5 individuals? #
  if (N < 5)
    return("Error: Group(s) size must be greater than or equal to 5 individuals for carrying out the jackknife.")
  source <- c("Actor Var","Partner Var","Relationship Var","Actor-Partner Cov","Dyadic Cov")
  jackknife <- array(dim=c(N,5,T,G),dimnames=list(names,source,times,groups),0.)
  jackknife.t.statistic <- array(dim=c(T,dim(jackknife)[2],G),dimnames=list(times,source,groups))
  jackknife.p.value <- array(dim=c(T,dim(jackknife)[2],G),dimnames=list(times,source,groups))
  jackknife.mean <- array(dim=c(T,dim(jackknife)[2],G),0.)
  jackknife.variance <- array(dim=c(T,dim(jackknife)[2],G),0.)
  Xtemp <- array(dim=c(N-1,N-1,T,G),0.)
  variances <- SRM.variances(X)
  variance.actor <- variances[[3]]
  variance.partner <- variances[[4]]
  variance.relationship <- variances[[5]]
  actorpartner.covariance <- variances[[6]]
  dyadic.covariance <- variances[[7]]
  for (i in 1:N){
     Xtemp <- as.numeric(X$data[-i,-i,,])
     Xtemp <- prepare.SRM.matrix(Xtemp,(N-1),T,G)
     variancestemp <- SRM.variances(Xtemp)
        for (k in 1:T){
           for (l in 1:G){
              jackknife[i,1,k,l] <- N*variance.actor[k,l]-variancestemp$actor.variance[k,l]*(N-1)
              jackknife[i,2,k,l] <- N*variance.partner[k,l]-variancestemp$partner.variance[k,l]*(N-1)
              jackknife[i,3,k,l] <- N*variance.relationship[k,l]-variancestemp$relationship.variance[k,l]*(N-1)
              jackknife[i,4,k,l] <- N*actorpartner.covariance[k,l]-variancestemp$actorpartner.covariance[k,l]*(N-1)
              jackknife[i,5,k,l] <- N*dyadic.covariance[k,l]-variancestemp$dyadic.covariance[k,l]*(N-1)}}}
  for (k in 1:T){
     for (l in 1:G){
        jackknife.mean[k,,l] <- apply(jackknife[,,k,l],2,mean)
        jackknife.variance[k,,l] <- apply(jackknife[,,k,l],2,var)/N
        jackknife.t.statistic[k,,l] <- jackknife.mean[k,,l]/sqrt(jackknife.variance[k,,l])}}
  df <- N-1
    for (k in 1:T){
       for (m in 1:(dim(jackknife)[2])){
          for (l in 1:G){
             if ( !is.na(jackknife.t.statistic[k,m,l])){
             if (sign(jackknife.t.statistic[k,m,l])>0)
               jackknife.p.value[k,m,l] <- pt(jackknife.t.statistic[k,m,l],df,lower.tail=FALSE)
             else jackknife.p.value[k,m,l] <- pt(jackknife.t.statistic[k,m,l],df)}}}}
  list(t.statistic=jackknife.t.statistic,df=df,two.tailed.p.value=2*jackknife.p.value)
}