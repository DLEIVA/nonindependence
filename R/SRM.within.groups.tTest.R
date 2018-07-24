
# Testing SRM effects by means of a within-group t Test #  

SRM.within.groups.tTest <- function(X){
  if (class(X) != "srmRRMatrix")
    stop("Data input need to be a srmRRMatrix object")
  N <- X$individuals
  T <- X$ntimes
  G <- X$ngroups
  names <- X$names
  times <- X$times
  groups <- X$groups
  # Computing some terms that will be used below #
  h1 <- 2*(G**2*N**6-7*G**2*N**5+18*G**2*N**4+10*G*N**4-20*G**2*N**3-
        46*G*N**3+8*G**2*N**2+70*G*N**2+24*N**2-36*G*N-48*N+32)
  h2 <- 2*(G**2*N**5-5*G**2*N**4+2*G*N**4+8*G**2*N**3-6*G*N**3-
        4*G**2*N**2+14*G*N**2+8*N**2-20*G*N-16*N+32)
  h3 <- 4*(G**2*N**5-5*G**2*N**4+8*G**2*N**3+12*G*N**3-
        4*G**2*N**2-38*G*N**2+28*G*N+32*N-32)
  h4 <- (N**2)*((N-2)**2)*(G*N-G+2)*(G*N**2-3*G*N+4)*(G*N**2-3*G*N+
        2*G+4)
  h5 <- G**3*N**7-7*G**3*N**6+19*G**3*N**5+10*G**2*N**5-25*G**3*N**4-
        54*G**2*N**4-4*G*N**4+16*G**3*N**3+120*G**2*N**3+44*G*N**3-
        4*G**3*N**2-132*G**2*N**2-124*G*N**2-16*N**2+56*G**2*N+
        168*G*N+32*N-64*G-64
  h6 <- G**3*N**7-7*G**3*N**6-2*G**2*N**6+19*G**3*N**5+26*G**2*N**5-
        25*G**3*N**4-100*G**2*N**4-20*G*N**4+16*G**3*N**3+176*G**2*N**3+
        124*G*N**3-4*G**3*N**2-156*G**2*N**2-236*G*N**2-48*N**2+56*G**2*N+
        200*G*N+96*N-64*G-64
  # Estimating mean and variance of round robin parameters #
  variances <- SRM.variances(X)
  variance.actor <- unlist(variances[[3]])
  variance.partner <- unlist(variances[[4]])
  variance.relationship <- unlist(variances[[5]])
  actorpartner.covariance <- unlist(variances[[6]])
  apcor <- unlist(SRM.generalized.reciprocity(X))
  dyadcor <- unlist(SRM.dyadic.reciprocity(X))
  dyadic.covariance <- unlist(variances[[7]])
  actor.variance.mean <- apply(variance.actor,1,function(x) mean(x,na.rm=TRUE))
  partner.variance.mean <- apply(variance.partner,1,function(x) mean(x,na.rm=TRUE))
  relationship.variance.mean <- apply(variance.relationship,1,function(x) mean(x,na.rm=TRUE))
  actorpartner.covariance.mean <- apply(actorpartner.covariance,1,function(x) mean(x,na.rm=TRUE))
  dyadic.covariance.mean <- apply(dyadic.covariance,1,function(x) mean(x,na.rm=TRUE))
  source <- c("Actor Var","Partner Var","Relationship Var","Actor-Partner Cov","Dyadic Cov")
  SRM.variance.parameters <- array(dim=c(T,5),dimnames=list(times,source),0.)
  within.group.t.statistic <- array(dim=c(T,5),dimnames=list(times,source),0.)
  within.group.p.value <- array(dim=c(T,5),dimnames=list(times,source),0.)
  for (k in 1:T){
  SRM.variance.parameters[k,1] <- ((2*(actor.variance.mean[k])**2)/(G*N-G+2))+(((4*actor.variance.mean[k])*
  ((N-1)*relationship.variance.mean[k]+dyadic.covariance.mean[k]))/(N*(N-2)*(G*N-G+2)))+
  ((h1*(relationship.variance.mean[k])**2+h2*(dyadic.covariance.mean[k])**2+h3*
  (relationship.variance.mean[k]*dyadic.covariance.mean[k]))/h4)
  SRM.variance.parameters[k,2] <- ((2*(partner.variance.mean[k])**2)/(G*N-G+2))+(((4*partner.variance.mean[k])*
  ((N-1)*relationship.variance.mean[k]+dyadic.covariance.mean[k]))/(N*(N-2)*(G*N-G+2)))+
  ((h1*(relationship.variance.mean[k])**2+h2*(dyadic.covariance.mean[k])**2+h3*
  (relationship.variance.mean[k]*dyadic.covariance.mean[k]))/h4)
  SRM.variance.parameters[k,3] <- SRM.variance.parameters[k,5] <- ((2*((G*N**2)-(3*G*N)+G+4)*
  (relationship.variance.mean[k]**2+dyadic.covariance.mean[k]**2))/(((G*N**2)-(3*G*N)+4)*
  ((G*N**2)-(3*G*N)+2*G+4)))+((4*G*relationship.variance.mean[k]*dyadic.covariance.mean[k])/(((G*N**2)-(3*G*N)+4)*
  ((G*N**2)-(3*G*N)+2*G+4)))
  SRM.variance.parameters[k,4] <- ((G*N-G-2)*actorpartner.covariance.mean[k]**2+G*(N-1)*
  actor.variance.mean[k]*partner.variance.mean[k])/((G*N-G+2)*(G*N-G-1))+(h5*relationship.variance.mean[k]**2+
  h6*dyadic.covariance.mean[k]**2+h3*relationship.variance.mean[k]*dyadic.covariance.mean[k])/
  ((G*N-G-1)*h4)+(G*(N-1)*(actor.variance.mean[k]+partner.variance.mean[k])*((N-1)*
  relationship.variance.mean[k]+dyadic.covariance.mean[k]))/(N*(N-2)*(G*N-G+2)*(G*N-G-1))+
  (2*(G*N-G-2)*actorpartner.covariance.mean[k]*(relationship.variance.mean[k]+(N-1)*dyadic.covariance.mean[k]))/
  (N*(N-2)*(G*N-G+2)*(G*N-G-1))}
  SRM.relative.parameters <- array(dim=c(T,5),dimnames=list(times,source),0.)
  SRM.mean.parameters <- array(dim=c(T,5),dimnames=list(times,source),c(actor.variance.mean,
  partner.variance.mean,relationship.variance.mean,actorpartner.covariance.mean,dyadic.covariance.mean))
  SRM.relative.parameters[1:3] <- SRM.mean.parameters[1:3]/sum(unlist(SRM.mean.parameters[1:3]),na.rm=TRUE)
  SRM.relative.parameters[4] <- apply(apcor,1,function(x) mean(x,na.rm=TRUE))
  SRM.relative.parameters[5] <- apply(dyadcor,1,function(x) mean(x,na.rm=TRUE))  
  SRM.stderror.parameters <- sqrt(SRM.variance.parameters)
  for (k in 1:T){
  within.group.t.statistic[k,1] <- actor.variance.mean[k]/sqrt(SRM.variance.parameters[k,1])
  within.group.t.statistic[k,2] <- partner.variance.mean[k]/sqrt(SRM.variance.parameters[k,2])
  within.group.t.statistic[k,3] <- relationship.variance.mean[k]/sqrt(SRM.variance.parameters[k,3])
  within.group.t.statistic[k,4] <- actorpartner.covariance.mean[k]/sqrt(SRM.variance.parameters[k,4])
  within.group.t.statistic[k,5] <- dyadic.covariance.mean[k]/sqrt(SRM.variance.parameters[k,5])}
  df <- G*(N-1)
  for (k in 1:T){ 
     for (m in 1:(dim(within.group.t.statistic)[2])){
        if ( !is.na(within.group.t.statistic[k,m])){
        if (sign(within.group.t.statistic[k,m]) > 0)
          within.group.p.value[k,m] <- pt(within.group.t.statistic[k,m],df,lower.tail=FALSE)
        else within.group.p.value[k,m] <- pt(within.group.t.statistic[k,m],df)}}}
  res <- list(estimates=SRM.mean.parameters,standard.values=SRM.relative.parameters,standard.error=SRM.stderror.parameters,
  t.value=within.group.t.statistic,df=df,p.value=within.group.p.value)
   class(res) <-"srmWithintTest"
  res
}
