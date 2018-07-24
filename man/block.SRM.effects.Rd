\name{block.SRM.effects}

\alias{block.SRM.effects}

\title{Estimating SRM effects in block designs}

\description{
A function for obtaining SRM effects in block designs.
}

\usage{
   block.SRM.effects(X)
}

\arguments{
  \item{X}{Original block sociomatrix.}
}

\details{
\code{block.SRM.effects} estimates actor,partner and relationship effects for SRM block designs.
}

\value{
\code{block.SRM.effects} returns a list containing the following components:

\item{actor.effects}{ Estimation of the actor effects for the SRM block data.}
\item{partner.effects}{ Estimation of the partner effects for the SRM block data.}
\item{relationship.effects}{ Estimation of the relationship effects for the SRM block data.}
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
