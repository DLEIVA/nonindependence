# Function for preparing data in order to be analized as a SRM block design #

prepare.block.matrix <- function (X,subgroup1,subgroup2,numb.times,numb.groups,
                                 names1=NULL,names2=NULL,times=NULL,groups=NULL,
                                 subgroups=NULL){
  numb.individuals <- subgroup1 + subgroup2
  if (is.null(names1)) names1 <- paste('Ind.',c(1:(subgroup1)))
    else names1 <- names1
  if (is.null(names2)) names2 <- paste('Ind.',c((numb.individuals-subgroup2+1):numb.individuals))
    else names2 <- names2
  if (is.null(times)) times <- paste('Time',c(1:numb.times))
    else times <- times
  if (is.null(groups)) groups <- paste('Group',c(1:numb.groups))
    else groups <- groups
  if (is.null(subgroups)) subgroups <- paste('Subgroup',c(1,2))
  else subgroups <- subgroups
    X <- array(X,c(numb.individuals,numb.individuals,numb.times,numb.groups),
         dimnames=list(c(names1,names2),c(names1,names2),times,groups))
    temp <- array(0,c(numb.individuals,numb.individuals,numb.times,numb.groups))
    for (i in 1:(dim(X)[1])){
       for (j in 1:(dim(X)[2])){
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
    res<-list(data=X,individuals1=subgroup1,individuals2=subgroup2,
              individuals=numb.individuals,ntimes=numb.times,ngroups=numb.groups,
              names1=names1,names2=names2,times=times,groups=groups,
              subgroups=subgroups)
    class(res) <- 'srmBlockMatrix'
    res
}