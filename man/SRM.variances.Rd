\name{SRM.variances}

\alias{SRM.variances}

\title{Estimating SRM absolute variances in round robin designs}

\description{
A function for obtaining SRM absolute variances in round robin designs.
}

\usage{
   SRM.variances(X)
}

\arguments{
  \item{X}{Original round robin sociomatrix.}
}

\details{
\code{SRM.variances} estimates actor,partner and relationship variances, as well as actor-partner and relationship covariances, for SRM round robin designs.
}

\value{
\code{SRM.variances} returns a list containing the following components:

\item{actor.variance}{ Estimation of the actor variance for the SRM round robin data.}
\item{partner.variance}{ Estimation of the partner variance for the SRM round robin data.}
\item{relationship.variance}{ Estimation of the relationship variance for the SRM robin data.}
\item{actorpartner.covariance}{ Estimation of the actor-partner covariance for the SRM round robin data.}
\item{dyadic.covariance}{ Estimation of the relationship covariance for the SRM round robin data.}
}

\references{

Kenny, D. A., Kashy, D. A., & Cook, W. L. (2006). \emph{Dyadic data analysis}. New York: Guilford Press.

}

\author{
David Leiva <dleivaur@ub.edu>, Antonio Solanas <antonio.solanas@ub.edu> & Lluís Salafranca <lsalafranca@ub.edu>.
}

\seealso{ \code{\link{SRM.effects}},\code{\link{SRM.relative.variances}}
}

\keyword{htest}

