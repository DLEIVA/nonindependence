\name{block.SRM.relative.variances}

\alias{block.SRM.relative.variances}

\title{Estimating SRM relative variances in block designs}

\description{
A function for obtaining SRM relative variances in block designs.
}

\usage{
   block.SRM.relative.variances(X)
}

\arguments{
  \item{X}{Original block sociomatrix.}
}

\details{
\code{block.SRM.relative.variances} estimates relative actor,partner and relationship variances for SRM block designs.
}

\value{
\code{block.SRM.relative.variances} returns a list containing the following components:

\item{relative.actor.variance}{ Estimation of the relative actor variance for the SRM block data.}
\item{relative.partner.variance}{ Estimation of the relative partner variance for the SRM block data.}
\item{relative.relationship.variance}{ Estimation of the relative relationship variance for the SRM block data.}
}

\references{

Kenny, D. A., Kashy, D. A., & Cook, W. L. (2006). \emph{Dyadic data analysis}. New York: Guilford Press.

}

\author{
David Leiva <dleivaur@ub.edu>, Antonio Solanas <antonio.solanas@ub.edu> & Lluís Salafranca <lsalafranca@ub.edu>.
}

\seealso{ \code{\link{block.SRM.effects}},\code{\link{block.SRM.variances}}
}

\keyword{htest}

