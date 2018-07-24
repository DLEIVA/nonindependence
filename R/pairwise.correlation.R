# A function for estimating dyadic interdependence in standard dyadic designs with indistinguishable dyad members and interval responses by means of the pairwise correlation #

pairwise.correlation <- function(dataset,conf.int=0.95){
  alpha <- 1 - conf.int
  double.entry.dataset <- array(c(c(1:(2*length(dataset[,1]))),c(dataset[,2],dataset[,3]),c(dataset[,3],dataset[,2])),dim=c((2*length(dataset[,1])),
                          length(dataset)))
  colnames(double.entry.dataset) <- colnames(dataset)
  dataset <- double.entry.dataset
  N<-length(na.omit(dataset)[,2])
  na1<- sum(is.na(dataset[,2]))
  na2<- sum(is.na(dataset[,3]))
  dataset<-na.omit(dataset)
  summary.statistics <- array(c(N,na1,min(dataset[,2]),max(dataset[,2]),mean(dataset[,2]),
                        sd(dataset[,2]),quantile(dataset[,2],.25,names=F),quantile(dataset[,2],.50,names=F),
                        quantile(dataset[,2],.75,names=F),N,na2,min(dataset[,3]),max(dataset[,3]),
                        mean(dataset[,3]),sd(dataset[,3]),quantile(dataset[,3],.25,names=F),quantile(dataset[,3],.50,names=F),
                        quantile(dataset[,3],.75,names=F)),dim=c(9,2),dimnames=list(c("N:", "NAs:","Min:","Max:","Mean:","Sd:",
                        "25th Pctl:","50th Pctl:", "75th Pctl:"),c(colnames(dataset[,2:3]))))
  rp <- cor(dataset[,2],dataset[,3])
  standard.error <- 1/sqrt(N)/2
  z <- rp/standard.error
  p.value <- 1-pnorm(abs(z),0,1)
  rp.low <- rp - (qnorm((1-alpha/2),0,1)/sqrt(N)/2)
  rp.upper <- rp + (qnorm((1-alpha/2),0,1)/sqrt(N)/2)
  res <- list(call=match.call(),data=dataset,stats=summary.statistics,rp=rp,alpha=alpha,rpCI=c(rp.low,rp.upper),
       zVal=z,stdErr=standard.error,pval=2*p.value)
  class(res) <- "rp"
  res
}
