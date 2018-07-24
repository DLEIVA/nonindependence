# A function for estimating dyadic interdependence in standard dyadic designs with indistinguishable dyad members and interval responses #

indistinguishable.dyad <- function(dataset,conf.int=0.95){
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
  MSb <- 2*var(apply(dataset[,2:3],1,mean))
  MSw <- sum((dataset[,2]-dataset[,3])**2)/(2*N)
  ICC <- (MSb-MSw)/(MSb+MSw)
  alpha <- 1 - conf.int
  if ( MSb > MSw ) {
    F.statistic <- MSb/MSw
    df1 <- N-1
    df2 <- N}
    if (MSw > MSb ){
    F.statistic <- MSw/MSb
    df1 <- N
    df2 <- N-1}
    p.value <- pf(F.statistic,df1,df2,lower.tail=FALSE)
  fl <- MSb/MSw/qf(1-alpha/2,df1,df2)
  fu <- MSb/MSw*qf(1-alpha/2,df1,df2)
  icc.low <- (fl-1)/(fl+1)
  icc.upper <- (fu-1)/(fu+1)
  res=list(call=match.call(),data=dataset,stats=summary.statistics,
           intracor=ICC,MSb=MSb,MSw=MSw,Fstat=F.statistic,df1=df1,df2=df2,pval=2*p.value,
           alpha=alpha,iccCI=c(icc.low,icc.upper))
  class(res)="dyadicc"
  res
}
