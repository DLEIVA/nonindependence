\name{prepare.block.matrix}

\alias{prepare.block.matrix}

\title{Preparing SRM block data}

\description{
A function for preparing data in order to be analized as a SRM block design.
}

\usage{
   prepare.block.matrix(X,subgroup1,subgroup2,numb.times,numb.groups,names1=NULL,names2=NULL,times=NULL,groups=NULL)
}

\arguments{
  \item{X}{Original block sociomatrix.}
  \item{subgroup1}{Number of individuals of the first subgroup of the block data.}
  \item{subgroup2}{Number of individuals of the second subgroup of the block data.}
  \item{numb.times}{Number of occasions in which the dyads are analyzed.}
  \item{numb.groups}{Number of groups in the block data.}
  \item{names1}{ (optional) Character array with the names of individuals in the first subgroup.}
  \item{names2}{ (optional) Character array with the names of individuals in the second subgroup.}
  \item{times}{ (optional) Character array with the names of the occasions of data collection.}
  \item{groups}{ (optional) Character array with the names of groups in the block data.}
}

\details{
\code{prepare.block.matrix} organizes data according a SRM block design to be used by the several functions for SRM block designs included in the package \code{nonindependence}.
}

\value{
\item{X}{ Sociomatrix data organized as an SRM block design.}
}

\references{

Kenny, D. A., Kashy, D. A., & Cook, W. L. (2006). \emph{Dyadic data analysis}. New York: Guilford Press.

}

\author{
David Leiva <dleivaur@ub.edu>, Antonio Solanas <antonio.solanas@ub.edu> & Lluís Salafranca <lsalafranca@ub.edu>.
}

\seealso{ \code{\link{block.SRM.effects}},\code{\link{block.SRM.variances}},\code{\link{block.SRM.relative.variances}},\code{\link{block.SRM.reliability}},\code{\link{block.SRM.dyadic.reciprocity}},\code{\link{block.SRM.generalized.reciprocity}},\code{\link{block.SRM.dyadic.reciprocity.tTest}},\code{\link{block.SRM.between.groups.tTest}},\code{\link{block.SRM.jackknife}}
}

\keyword{htest}
