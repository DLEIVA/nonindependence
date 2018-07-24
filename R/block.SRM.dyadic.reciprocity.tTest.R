
# Testing dyadic covariance in block designs by means of r Pearson correlation test #

block.SRM.dyadic.reciprocity.tTest <- function(X,symmetric=FALSE){
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
  dyadic.reciprocity <- block.SRM.dyadic.reciprocity(X,symmetric)
  t.value <- array(dim=c(T,G),dimnames=list(times,groups),0.)
  for (k in 1:T){
     for (l in 1:G){
        t.value[k,l] <- dyadic.reciprocity[k,l]*sqrt(G1*(G2-1)-1)/
                        sqrt(1-dyadic.reciprocity[k,l]**2)}}
  df <- (G1-1)*(G2-1)-1
  t.p.value <- array(dim=c(T,G),dimnames=list(times,groups),NA)
  for (k in 1:T){
     for (l in 1:G){
        if (sign(t.value[k,l])>0)
          t.p.value[k,l] <- pt(t.value[k,l],df,lower.tail=FALSE)
          else t.p.value[k,l] <- pt(t.value[k,l],df)}}
  list(t.value=t.value,two.tailed.p.value=2*t.p.value)
}