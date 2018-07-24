print.srmRRVars <- function(x,digits=max(4,getOption("digits")-4),...)
{
  cat("\n")
  cat("SRM Variances: Round Robin Design")
  cat("\n\n")
  T <- dim(x$actor.variance)[1]
  G <- dim(x$actor.variance)[2]
  times <- unlist(dimnames(x$actor.variance)[1])
  groups <- unlist(dimnames(x$actor.variance)[2])
  group <- rep(groups,T)
  time <- rep(times,G)
  actor <- array(x$actor.variance)
  partner <- array(x$partner.variance)
  relationship <- array(x$relationship.variance)
  actor.partner <- array(x$actorpartner.covariance)
  dyadic <- array(x$dyadic.covariance)
  results<-data.frame(group,time,actor,partner,relationship,actor.partner,dyadic)
  print(results)
  cat("\n\n")
  invisible(x)
}