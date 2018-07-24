# A function for estimating dyadic interdependence in standard dyadic designs with distinguishable dyad members and interval responses #

distinguishable.dyad <- function(dataset,conf.int=0.95){
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
  
  r <- cor(dataset[,2],dataset[,3])
  df <- N-2
  alpha <- 1 - conf.int
  t.statistic <- r*sqrt(N-2)/(sqrt(1-r**2))
  p.value <- 1 - pt(t.statistic,df)
  zr.value <- 0.5*log((1+r)/(1-r))
  zr.low <- zr.value - qnorm((1-alpha/2),0,1)/sqrt(N-3)
  zr.upper <- zr.value + qnorm((1-alpha/2),0,1)/sqrt(N-3)
  r.low <- (exp(1)**(2*zr.low)-1)/(exp(1)**(2*zr.low)+1)
  r.upper <- (exp(1)**(2*zr.upper)-1)/(exp(1)**(2*zr.upper)+1)
  res <- list(call=match.call(),data=dataset,stats=summary.statistics,r=r,rCI=c(r.low,r.upper),
         t=t.statistic,df=df,pval=2*p.value)
  class(res) <- 'dyadr'
  res
}
