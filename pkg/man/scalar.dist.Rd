\name{scalar.dist}
\Rdversion{1.4.0}
\alias{scalar.dist}
\alias{scalar.dist,vector-method}
\title{
Calculate scalar distance
}
\description{
Calculate the scalar distance between elements of a matrix.
}
\usage{
scalar.dist(x, ...)
}

\arguments{
  \item{x}{
A vector of numeric values.
}
  \item{\dots}{
Additional arguments, not used at this time.
}
}

\details{
Take the absolute difference between all elements in a vector, and return a matrix of the distances.
}

\value{
a matrix object
}

\author{
Cole Beck
}

\examples{
scalar.dist(1:10)
}
