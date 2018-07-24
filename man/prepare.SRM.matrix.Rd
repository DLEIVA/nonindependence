\name{prepare.SRM.matrix}

\alias{prepare.SRM.matrix}

\title{Preparing SRM round robin data}

\description{
A function for preparing data in order to be analized as a SRM round robin design.
}

\usage{
   prepare.SRM.matrix(X,numb.individuals,numb.times,numb.groups,names=NULL,times=NULL,groups=NULL)
}

\arguments{
  \item{X}{Original round robin sociomatrix.}
  \item{numb.individuals}{Number of individuals in the round robin data.}
  \item{numb.times}{Number of occasions in which the dyads are analyzed.}
  \item{numb.groups}{Number of groups in the round robin data.}
  \item{names}{ (optional) Character array with the names of individuals in the group.}
  \item{times}{ (optional) Character array with the names of the occasions of data collection.}
  \item{groups}{ (optional) Character array with the names of groups in the round robin data.}
}

\details{
\code{prepare.SRM.matrix} organizes data according a SRM round robin design to be used by the several functions for SRM round robin designs included in the package \code{nonindependence}.
}

\value{
\item{X}{ Sociomatrix data organized as an SRM round robin design.}
}

\references{

Kenny, D. A., Kashy, D. A., & Cook, W. L. (2006). \emph{Dyadic data analysis}. New York: Guilford Press.

}

\author{
David Leiva <dleivaur@ub.edu>, Antonio Solanas <antonio.solanas@ub.edu> & Lluís Salafranca <lsalafranca@ub.edu>.
}

\seealso{ \code{\link{SRM.effects}},\code{\link{SRM.variances}},\code{\link{SRM.relative.variances}},\code{\link{SRM.reliability}},\code{\link{SRM.dyadic.reciprocity}},\code{\link{SRM.generalized.reciprocity}},\code{\link{SRM.dyadic.covariance.Ftest}},\code{\link{SRM.between.groups.tTest}},\code{\link{SRM.jackknife}}
}

\keyword{htest}
