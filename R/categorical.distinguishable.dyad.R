# A function for estimating dyadic interdependence in standard dyadic designs with distinguishable dyad members and categorical responses #

categorical.distinguishable.dyad <- function(dataset,conf.int=0.95){
  tabulated.data <- table(dataset[,2],dataset[,3])
  alpha <- 1 - conf.int
  no <- sum(diag(tabulated.data))
  ne <- sum(margin.table(tabulated.data,1)*margin.table(tabulated.data,2)/sum(tabulated.data))
  kappa <- (no - ne)/(sum(tabulated.data)-ne)
  table.props <- prop.table(tabulated.data)
  prop.margins <-array(c(margin.table(table.props,1),margin.table(prop.table(tabulated.data),2)),dim=c(nrow(tabulated.data),2))
  sumterm1 <- 0.
  sumterm2 <- 0.
  sumterm3 <- 0.
  pe <- ne/sum(tabulated.data)
  for (i in 1:nrow(table.props)){
     sumterm1 <- sumterm1+(table.props[i,i]*(1-sum(prop.margins[i,])*(1-kappa))**2+(1-kappa)**2)
     sumterm3 <- sumterm3+(prop.margins[i,1]*prop.margins[i,2]*sum(prop.margins[i,]))} 
  for (i in 1:nrow(table.props)){
     for (j in 1:ncol(table.props)){
        sumterm2 <- sumterm2+(table.props[i,j]*(prop.margins[i,2]+prop.margins[j,1])**2-(kappa-pe*(1-kappa))**2)}} 
  se.kappa2 <- sqrt((sumterm1*sumterm2)/(sum(tabulated.data)*(1-pe)**2))
  se.kappa <- sqrt((pe+pe**2-sumterm3)/(sum(tabulated.data)*(1-pe)**2))
  tabulated.data <- addmargins(tabulated.data)
  z.value <- kappa/se.kappa
  p.value <- 1-pnorm(abs(z.value),0,1)
  kappa.low <- kappa-qnorm((1-alpha/2),0,1)*se.kappa
  kappa.upper <- kappa+qnorm((1-alpha/2),0,1)*se.kappa
  res <-list(call=match.call(),contingencyTab=tabulated.data,kappa=kappa,alpha=alpha,kappaCI=c(kappa.low,kappa.upper),
       zval=z.value,stdError=se.kappa,pval=2*p.value)
  class(res) <- "dyadkappa"
  res
}