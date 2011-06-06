\name{nbpMatching-package}
\Rdversion{1.0}
\alias{nbpMatching-package}
\alias{nbpMatching}
\docType{package}
\title{
Nonbipartite Matching
}
\description{
This package will take an input distance matrix and generate the set of pairwise matches that minimizes the sum of distances between the pairs by running nonbimatch.
}

\author{
Bo Lu, Robert Greevy, Cole Beck

Maintainer: Cole Beck \email{cole.beck@vanderbilt.edu}
}

\references{
Lu, Greevy, Xu, and Beck. Optimal Nonbipartite Matching and its Statistical Applications. (under review)
}

\examples{
plainmatrix<-as.matrix(dist(sample(1:25, 8, replace=TRUE)))
diag(plainmatrix) <- 99999  # setting diagonal to an infinite distance for
                            # pedagogical reasons (the diagonal may be left
                            # as zero)
mdm<-distancematrix(plainmatrix)
res<-nonbimatch(mdm)

\dontrun{
try a large matrix
nonbimatch(distancematrix(as.matrix(dist(sample(1:10^8, 5000, replace=TRUE)))))
}
}
\keyword{ package }
\keyword{ cluster }
\keyword{ array }
