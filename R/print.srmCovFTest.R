print.srmCovFTest <- function(x,digits=max(4,getOption("digits")-4),...)
{
  cat("\n")
  cat("SRM Dyadic Covariance F Test: Round Robin Design")
  cat("\n\n")
  T <- dim(x$MSb)[1]
  G <- dim(x$MSb)[2]
  times <- unlist(dimnames(x$MSb)[1])
  groups <- unlist(dimnames(x$MSb)[2])
  group <- rep(groups,T)
  time <- rep(times,G)
  MSb <- x[[1]]
  MSw <- x[[2]]
  df1 <- x[[4]]
  df2 <- x[[5]]
  F <- x[[3]]
  pval <- x[[6]]
  results<-data.frame(group,time,MSb,MSw,df1,df2,F, pval)
  colnames(results)<-c("Group","Time","MSb","MSw","df1","df2","F value",if (sign(F)>0) "Pr(> F)" else "Pr(< F)")
  rownames(results) <- ''
  print(results)
  cat("\n\n")
  invisible(x)
}