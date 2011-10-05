\name{nbpMatching-package}
\Rdversion{1.3.2}
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
# create a covariate matrix
df <- data.frame(id=LETTERS[1:25], val1=rnorm(25), val2=rnorm(25))
# create distances
df.dist <- gendistance(df, idcol=1)
# create distancematrix object
df.mdm <- distancematrix(df.dist)
# create matches
df.match <- nonbimatch(df.mdm)
# review quality of matches
df.qom <- qom(df.dist$cov, df.match$matches)

# some helper functions are available
# runner -- start with the covariate, run through the entire process
df.1 <- runner(df, idcol=1)
# full.qom -- start with the covariate, generate a full quality of match report
df.2 <- full.qom(df)

\dontrun{
try a large matrix
nonbimatch(distancematrix(as.matrix(dist(sample(1:10^8, 5000, replace=TRUE)))))
}
}
\keyword{ package }
\keyword{ cluster }
\keyword{ array }
