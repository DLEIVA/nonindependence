
# Testing SRM effects in block designs by means of Jackknife tests #

block.SRM.jackknife <-  function(X){
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
  # Is group(s) size greater than or equal to 6 individuals? #
  if (N < 6)
    return("Error: Group(s) size must be greater than or equal to 6 individuals for carrying out the jackknife for a block design.")
  if (G1 < 3 | G2 < 3)
    return("Error: Subgroups must be greater than or equal to 3 individuals for carrying out the jackknife for a block design.")
  source <- c("Actor Var","Partner Var","Relationship Var","Actor-Partner Cov","Dyadic Cov")
  jackknife <- array(dim=c(length(subgroups),N,5,T,G),dimnames=list(subgroups,c(names1,names2),
                     source,times,groups),0.)
  jackknife.t.statistic <- array(dim=c(length(subgroups),T,dim(jackknife)[3],
                                 G),dimnames=list(subgroups,times,source,groups),0.)
  jackknife.p.value <- array(dim=c(length(subgroups),T,dim(jackknife)[3],G),
                             dimnames=list(subgroups,times,source,groups),NA)
  jackknife.mean <- array(dim=c(length(subgroups),T,dim(jackknife)[3],G),0.)
  jackknife.variance <- array(dim=c(length(subgroups),T,dim(jackknife)[3],G),0.)
  Xtemp <- array(dim=c(N-1,N-1,T,G),0.)
  variances <- block.SRM.variances(X)
  variance.actor <- variances$actor.variance
  variance.partner <- variances$partner.variance
  variance.relationship <- variances$relationship.variance
  actorpartner.covariance <- variances$actorpartner.covariance
  dyadic.covariance <- variances$dyadic.covariance
  for (m in 1:length(subgroups)){
     count <- 1
     for (i in 1:N){
        if (count <= G1) {
          ind.temp.1 <- G1-1
          ind.temp.2 <- G2}
        else {
          ind.temp.1 <- G1
          ind.temp.2 <- G2-1}
        if (count <= G1){ 
          term1 <- G1 
          term2 <- ind.temp.1}
        else {term1 <- G2 
            term2 <- ind.temp.2}
        count <- count + 1
        Xtemp <- prepare.block.matrix(as.numeric(X$data[-i,-i,,]),N-1,N-1,T,G)
        variancestemp <- block.SRM.variances(Xtemp)
        for (k in 1:T){
           for (l in 1:G){
              jackknife[m,i,1,k,l] <- term1*variance.actor[m,k,l]-variancestemp$actor.variance[m,k,l]*term2
              jackknife[m,i,2,k,l] <- term1*variance.partner[m,k,l]-variancestemp$partner.variance[m,k,l]*term2
              jackknife[m,i,3,k,l] <- term1*variance.relationship[m,k,l]-variancestemp$relationship.variance[m,k,l]*term2
              jackknife[m,i,4,k,l] <- term1*dyadic.covariance[m,k,l]-variancestemp$dyadic.covariance[m,k,l]*term2
              jackknife[m,i,5,k,l] <- term1*actorpartner.covariance[k,l]-variancestemp$actorpartner.covariance[k,l]*term2}}}}
  for (m in 1:length(subgroups)){
    for (k in 1:T){
       for (l in 1:G){
          jackknife.mean[m,k,,l] <- apply(jackknife[m,,,k,l],2,mean)
          jackknife.variance[m,k,,l] <- apply(jackknife[m,,,k,l],2,var)/(N)
          jackknife.t.statistic[m,k,,l] <- jackknife.mean[m,k,,l]/sqrt(jackknife.variance[m,k,,l])}}}
  df1 <- G1-1
  df2 <- G2-2
  df <- array(dim=c(2,1),dimnames=list(subgroups),c(df1,df2))
  for (m in 1:length(subgroups)){    
     for (k in 1:T){
        for (n in 1:(dim(jackknife)[3])){
           for (l in 1:G){
              if (sign(jackknife.t.statistic[m,k,n,l])>0){
                if (m == 1)
                jackknife.p.value[m,k,n,l] <- pt(jackknife.t.statistic[m,k,n,l],df1,lower.tail=FALSE)
                if (m == 2)
                jackknife.p.value[m,k,n,l] <- pt(jackknife.t.statistic[m,k,n,l],df2,lower.tail=FALSE)}
              else{
                if (m == 1)
                jackknife.p.value[m,k,n,l] <- pt(jackknife.t.statistic[m,k,n,l],df1)
                if (m == 2)
                jackknife.p.value[m,k,n,l] <- pt(jackknife.t.statistic[m,k,n,l],df2)}}}}}
  list(t.statistic=jackknife.t.statistic,df=df,two.tailed.p.value=2*jackknife.p.value)
}