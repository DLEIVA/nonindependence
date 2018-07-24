\name{pairwise.correlation}

\alias{pairwise.correlation}

\title{Pairwise correlation for estimating nonindependence in indistinguishable standard dyadic designs with interval outcomes}

\description{
A function for estimating dyadic interdependence in standard dyadic designs with indistinguishable dyad members and interval responses by means of the pairwise correlation.
}

\usage{
   pairwise.correlation(dataset, conf.int=0.95)
}

\arguments{
  \item{dataset}{ Original data with dyadic responses.}
  \item{conf.int}{ Confidence level in order to estimate the confidence interval for the rp parameter. It is equal to 0.95 by default}
}

\details{
\code{pairwise.correlation} estimates dyadic interdependence for any standard dyadic design in which indistinguishable members and interval data are taken into account by means of a pairwise correlation.
}

\value{
\code{pairwise.correlation} returns a list containing the following components:

\item{descriptive.statistics}{ Summary of some descriptive statistics for the dyadic variables.}
\item{pairwise.correlation}{ rp coefficient for the original data.}
\item{alpha}{ Probability of the parameter not included in the interval estimation.}
\item{rp.confidence.interval}{ Lower and upper bounds of the rp parameter estimated at 1-alpha confidence level.}
\item{z.value}{ Statistic value for carrying out the statistical test regarding dyadic nonindependence.}
\item{standard.error}{ Standard error for the rp parameter under the null hypothesis of no dyadic correlation.}
\item{two.tailed.p.value}{ Statistical significance of the z statistic value under the null hypothesis of no dyadic correlation}
}

\references{

Kenny, D. A., Kashy, D. A., & Cook, W. L. (2006). \emph{Dyadic data analysis}. New York: Guilford Press.

}

\author{
David Leiva <dleivaur@ub.edu>, Antonio Solanas <antonio.solanas@ub.edu> & Lluís Salafranca <lsalafranca@ub.edu>.
}

\seealso{ \code{\link{indistinguishable.dyad}}
}

\keyword{htest}

