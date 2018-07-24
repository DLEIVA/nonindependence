\name{SRM.reliability}

\alias{SRM.reliability}

\title{Estimating reliability for the SRM estimates in round robin designs}

\description{
A function to compute reliability for the estimates of the actor and partner variances in round robin designs.
}

\usage{
   SRM.reliability(X)
}

\arguments{
  \item{X}{Original round robin sociomatrix.}
}

\details{
\code{SRM.reliability} estimates reliability for the actor and partner variances in SRM round robin designs.
}

\value{
\code{SRM.reliability} returns a list containing the following components:

\item{actor.reliability}{ Estimation of reliability of the actor variance for the SRM round robin data.}
\item{partner.reliability}{ Estimation of reliability of the partner variance for the SRM round robin data.}
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
