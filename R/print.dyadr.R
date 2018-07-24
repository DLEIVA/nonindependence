print.dyadr <- function(x,digits=max(4,getOption("digits")-4),...)
{
  cat("\n")
  cat("Distinguishable Members: Pearson's Correlation")
  cat("\n\n Call: \n")
  cat("",deparse(x$call), "\n\n")
  cat(" Data",  if (length(x$data[,1]) <= 5) ": " else " (5 first rows shown): ", "\n")
  print( if (length(x$data[,1]) >= 5) x$data[1:5,] else x$data[1:length(x$data),])
  cat("\n\n")
  cat("Descriptive Statistics","\n")
  print.table(x$stats, digits=digits)
  cat("\n\n")
  cat("t Test: ","\n")
  results <- cbind(x$t,as.integer(x$df),x$pval)
  colnames(results) <- c("t Value", "Df", "Pr(>|t|)")
  rownames(results) <- ""
  print.table(results,digits=digits)
  cat("\n\n")
  cat(paste((1-x$alpha)*100,"% Confidence Interval: ",sep=''),"\n")
  results <- t(x$rCI)
  colnames(results) <- c("Lower","Upper")
  rownames(results) <- ''
  print.table(results,digits=digits)
  cat("\n\n")
  cat("Sample estimate","\n")
  results <- cbind(x$r)
  colnames(results) <- "cor"
  rownames(results) <- ''
  print.table(results,digits=digits)
  cat("\n\n")
  invisible(x)
}