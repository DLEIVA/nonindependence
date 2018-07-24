\name{SRM.jackknife}

\alias{SRM.jackknife}

\title{Testing SRM effects in round robin designs by means of Jackknife tests}

\description{
A function for testing SRM effects in round robin designs by means of a Jackknife test.
}


\usage{
   SRM.jackknife(X)
}

\arguments{
  \item{X}{Original round robin sociomatrix.}
}

\details{
\code{SRM.jackknife} estimates statistical significance for the variance and covariance components obtained for SRM round robin designs.
}

\value{
\code{SRM.jackknife} returns a list containing the following components:

\item{t.statistic}{ t statistic for each variance component of the SRM block data estimated by means of a Jackknife procedure.}
\item{df}{Degrees of freedom for the statistical test.}
\item{two.tailed.p.value}{Statistical significance for the variance components of the SRM round robin data.}
}

\references{

Kenny, D. A., Kashy, D. A., & Cook, W. L. (2006). \emph{Dyadic data analysis}. New York: Guilford Press.
}

\author{
David Leiva <dleivaur@ub.edu>, Antonio Solanas <antonio.solanas@ub.edu> & Lluís Salafranca <lsalafranca@ub.edu>.
}

\seealso{ \code{\link{SRM.variances}}
}

\keyword{htest}

