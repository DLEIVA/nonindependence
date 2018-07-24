\name{SRM.relative.variances}

\alias{SRM.relative.variances}

\title{Estimating SRM relative variances in round robin designs}

\description{
A function for obtaining SRM relative variances in round robin designs.
}

\usage{
   SRM.relative.variances(X)
}

\arguments{
  \item{X}{Original round robin sociomatrix.}
}

\details{
\code{SRM.relative.variances} estimates relative actor,partner and relationship relative variances for SRM round robin designs.
}

\value{
\code{SRM.relative.variances} returns a list containing the following components:

\item{relative.actor.variance}{ Estimation of the relative actor variance for the SRM round robin data.}
\item{relative.partner.variance}{ Estimation of the relative partner variance for the SRM round robin data.}
\item{relative.relationship.variance}{ Estimation of the relative relationship variance for the SRM round robin data.}
}

\references{

Kenny, D. A., Kashy, D. A., & Cook, W. L. (2006). \emph{Dyadic data analysis}. New York: Guilford Press.

}

\author{
David Leiva <dleivaur@ub.edu>, Antonio Solanas <antonio.solanas@ub.edu> & Lluís Salafranca <lsalafranca@ub.edu>.
}

\seealso{ \code{\link{SRM.effects}},\code{\link{SRM.variances}}
}

\keyword{htest}


