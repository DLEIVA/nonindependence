print.srmRREffects <- function(x,digits=max(4,getOption("digits")-4),...)
{
  cat("\n")
  cat("SRM Effects: Round Robin Design")
  cat("\n\n")
  N <- dim(x$actor.effects)[1]
  T <- dim(x$actor.effects)[2]
  G <- dim(x$actor.effects)[3]
  names <- unlist(dimnames(x$actor.effects)[1])
  times <- unlist(dimnames(x$actor.effects)[2])
  groups <- unlist(dimnames(x$actor.effects)[3])
  group <- rep(groups,2*N*T)
  time <- rep(times,2*N*G)
  individual <- rep(names,2)
  type <- c(rep("actor",N*T*G),rep("partner",N*T*G))
  effect <- c (x$actor.effects,x$partner.effects)
  results<-data.frame(group,time,individual,type,effect)
  cat("Actor and Partner Effects: \n")
  print(results)
  cat("\n\n")
  results <- vector()
 
  for (j in 1:G)
  {  
    for (i in 1:T)
    {  
      mat <- x$relationship.effects[,,i,j]
      actornames <- c(rownames(mat)[(row(mat)*(row(mat)<col(mat)))[row(mat)<col(mat)]],
                    rownames(mat)[(row(mat)*(row(mat)>col(mat)))[row(mat)>col(mat)]])
      partnernames <- c(colnames(mat)[(col(mat)*(row(mat)<col(mat)))[row(mat)<col(mat)]],
                      colnames(mat)[(col(mat)*(row(mat)>col(mat)))[row(mat)>col(mat)]])
      effects <- mat[cbind(c((row(mat)*(row(mat)<col(mat)))[row(mat)<col(mat)],
                 (row(mat)*(row(mat)>col(mat)))[row(mat)>col(mat)]),
                 c((col(mat)*(row(mat)<col(mat)))[row(mat)<col(mat)],
                 (col(mat)*(row(mat)>col(mat)))[row(mat)>col(mat)]))]
      group <- rep(groups[j],N*(N-1)/2)
      time <- rep(times[i],N*(N-1)/2)
      res <- data.frame(group,time,actor=actornames,partner=partnernames,relationship=effects)
      results<- rbind(results,res)
    }
  }
  results<-results[order(results$actor,results$partner),]
  rownames(results) <- NULL
  cat("Relationship Effects: \n")
  print(results)
  cat("\n\n")
  invisible(x)
}