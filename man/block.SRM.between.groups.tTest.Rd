\name{block.SRM.between.groups.tTest}

\alias{block.SRM.between.groups.tTest}

\title{Testing SRM effects in block designs by means of a between-group t Test}

\description{
A function for testing SRM effects in block designs by means of a between-group t Test.
}


\usage{
   block.SRM.between.groups.tTest(X)
}

\arguments{
  \item{X}{Original block sociomatrix.}
}

\details{
\code{block.SRM.between.groups.tTest} estimates statistical significance for the variance and covariance components obtained for a SRM block design.
}

\value{
\code{block.SRM.between.groups.tTest} returns a list containing the following components:

\item{t.statistic}{ t statistic for each variance component of the SRM block data.}
\item{df}{Degrees of freedom for the statistical test.}
\item{two.tailed.p.value}{Statistical significance for the variance components of the SRM block data.}
}

\references{

Kenny, D. A., Kashy, D. A., & Cook, W. L. (2006). \emph{Dyadic data analysis}. New York: Guilford Press.
}

\author{
David Leiva <dleivaur@ub.edu>, Antonio Solanas <antonio.solanas@ub.edu> & Lluís Salafranca <lsalafranca@ub.edu>.
}

\seealso{ \code{\link{block.SRM.variances}}
}

\keyword{htest}
