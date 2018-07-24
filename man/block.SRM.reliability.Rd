\name{block.SRM.reliability}

\alias{block.SRM.reliability}

\title{Estimating reliability for the SRM estimates in block designs}

\description{
A function to compute reliability for the estimates of the actor and partner variances in block designs.
}

\usage{
   block.SRM.reliability(X)
}

\arguments{
  \item{X}{Original block sociomatrix.}
  }

\details{
\code{block.SRM.reliability} estimates reliability for the actor and partner variances in SRM block designs.
}

\value{
\code{block.SRM.reliability} returns a list containing the following components:

\item{actor.reliability}{ Estimation of reliability of the actor variance for the SRM block data.}
\item{partner.reliability}{ Estimation of reliability of the partner variance for the SRM block data.}
}

\references{

Kenny, D. A., Kashy, D. A., & Cook, W. L. (2006). \emph{Dyadic data analysis}. New York: Guilford Press.

}

\author{
David Leiva <dleivaur@ub.edu>, Antonio Solanas <antonio.solanas@ub.edu> & Lluís Salafranca <lsalafranca@ub.edu>.
}

\seealso{ \code{\link{block.SRM.relative.variances}}
}

\keyword{htest}
