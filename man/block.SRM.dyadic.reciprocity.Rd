\name{block.SRM.dyadic.reciprocity}

\alias{block.SRM.dyadic.reciprocity}

\title{Computing dyadic reciprocity in SRM block designs}

\description{
A function for computing dyadic reciprocity in SRM block designs.
}


\usage{
   block.SRM.dyadic.reciprocity(X,symmetric=FALSE)
}

\arguments{
  \item{X}{Original block sociomatrix.}
}

\details{
\code{block.SRM.dyadic.reciprocity} estimates dyadic reciprocity for symmetric and non-symmetric SRM block designs.
}

\value{
\item{dyadic.reciprocity}{Dyadic reciprocity measure for SRM block designs.}
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
