\name{quantile}
\Rdversion{1.3.6}
\alias{quantile}
\alias{quantile,distancematrix-method}
\title{
Quantile for upper-triangular values in distance matrix
}
\description{
Extend the stats quantile function for handling distancematrix objects.
}
\usage{
quantile(x, ...)
}

\arguments{
  \item{x}{
A distancematrix object.
}
  \item{\dots}{
Additional arguments, passed to stats:::quantile.
}
}

\details{
The upper.triangular values of the distance matrix object are passed to the quantile funtion.
}

\value{
numeric vector of quantiles corresponding to the given probabilities
}

\author{
Cole Beck
}

\examples{
plainmatrix<-as.matrix(dist(sample(1:25, 8, replace=TRUE)))
mdm<-distancematrix(plainmatrix)
quantile(mdm, probs=c(0.0, 0.25, 0.50, 0.75, 1.00))
}
