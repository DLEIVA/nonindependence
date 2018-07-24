print.srmWithintTest <- function(x,digits=max(4,getOption("digits")-4),...)
{
  cat("\n")
  cat("Within group t-Test: Round Robin Design")
  cat("\n\n")
  estimates <- as.numeric(x$estimates)
  std.value <- as.numeric(x$standard.values)
  std.error <- as.numeric(x$standard.error)
  t <- as.numeric(x$t.value)
  df <- as.numeric(x$df)
  p.value<- as.numeric(x$p.value)
  results<-data.frame(estimates,std.value,std.error,t,df,p.value)
  results[is.na(results)]<-NA
  results[is.na(results[,2]),3:6]<-NA
  rownames(results) <- c('actor var','partner var','relationship var','actor-partner cov','dyadic cov')
  print(results)
  cat("\n\n")
  invisible(x)
}