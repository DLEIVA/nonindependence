print.dyadkappa <- function(x,digits=max(4,getOption("digits")-4),...)
{
  cat("\n")
  cat("Distinguishable Members with categorical data: Kappa Index")
  cat("\n\n Call: \n")
  cat("",deparse(x$call), "\n\n")
  cat(" Contingency table: \n")
  print.table(x$contingencyTab)
  cat("\n\n")
  cat("Agreement Test: ","\n")
  results <- cbind(x$kappa,x$stdError,x$zval,x$pval)
  colnames(results) <- c("Kappa Value", "Asymp. Std Error", "z Value", "Approx. Sig.")
  rownames(results) <- ""
  print.table(results,digits=digits,scientific=FALSE)
  cat("\n\n")
  cat(paste((1-x$alpha)*100,"% Confidence Interval: ",sep=''),"\n")
  results <- t(x$kappaCI)
  colnames(results) <- c("Lower","Upper")
  rownames(results) <- ''
  print.table(results,digits=digits)
  cat("\n\n")
  invisible(x)
}