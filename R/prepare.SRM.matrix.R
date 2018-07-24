
prepare.SRM.matrix <- function(X,numb.individuals,numb.times,numb.groups,names=NULL,times=NULL,groups=NULL){
  if (is.null(names)) names <- paste('Ind.',c(1:numb.individuals))
  if (is.null(times)) times <- paste('Time',c(1:numb.times))
  if (is.null(groups)) groups <- paste('Group',c(1:numb.groups))
  X <- array(X,c(numb.individuals,numb.individuals,numb.times,numb.groups))
  
  temp <- array(0,c(numb.individuals,numb.individuals,numb.times,numb.groups))
  for (i in 1:(dim(X)[1])){
     for (j in 1:(dim(X)[1])){
        for (k in 1:(dim(X)[3])){
           for (l in 1:(dim(X)[4])){
           temp[i,j,k,l] <- X[j,i,k,l]}}}}
  count <- 1;
  count2 <- 1;
  for (k in 1:(dim(X)[3])){
        for (i in 1:(dim(X)[1])){
           if (count <= (dim(X)[3])){ 
             X[count2,,count,] <- temp[i,,k,]
             count <- count + 1}
           else {
             count <- 1
             count2 <- count2 + 1
             X[count2,,count,] <- temp[i,,k,]
             count <- count + 1}}}
  res<-list(data=X,individuals=numb.individuals,ntimes=numb.times,ngroups=numb.groups,
            names=names,times=times,groups=groups)
  class(res) <- 'srmRRMatrix'
  res
}