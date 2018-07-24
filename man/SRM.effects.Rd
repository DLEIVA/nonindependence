\name{SRM.effects}

\alias{SRM.effects}

\title{Estimating SRM effects in round robin designs}

\description{
A function for obtaining SRM effects in round robin designs.
}

\usage{
   SRM.effects(X)
}

\arguments{
  \item{X}{Original round robin sociomatrix.}
}

\details{
\code{SRM.effects} estimates actor,partner and relationship effects for SRM round robin designs.
}

\value{
\code{SRM.effects} returns a list containing the following components:

\item{actor.effects}{ Estimation of the actor effects for the SRM round robin data.}
\item{partner.effects}{ Estimation of the partner effects for the SRM round robin data.}
\item{relationship.effects}{ Estimation of the relationship effects for the SRM round robin data.}
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
