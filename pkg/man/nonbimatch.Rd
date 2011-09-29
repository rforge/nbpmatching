\name{nonbimatch}
\Rdversion{1.3.1}
\alias{nonbimatch}
\alias{nonbimatch,distancematrix-method}
\title{
Nonbipartite Matching
}
\description{
The nonbinmatch function creates the set of pairwise matches that minimizes the sum of distances between the pairs.  
}
\usage{
nonbimatch(mdm, threshold = NA, precision = 6, ...)
}

\arguments{
  \item{mdm}{
A distancematrix object.  See the distancematrix function.
}
  \item{threshold}{
An integer value, indicating the distance needed to create chameleon matches.
}
  \item{precision}{                                                                                                                        
The largest value in the matrix will have at most this many digits.  The default value is six.                                             
}
  \item{\dots}{
Additional arguments, these are not used.
}
}

\details{
The nonbinmatch function calls the Fortran code (Derigs) and set of pairwise matches that minimizes the sum of distances between the pairs.
}

\value{
    \item{matches }{data.frame containing matches}
    \item{halves }{data.frame containing each match}
    \item{total }{sum of the distances across all pairs}
    \item{mean }{mean distance for each pair}
}

\author{
Cole Beck
}

\seealso{
\code{\link{distancematrix}}
}

\examples{
plainmatrix<-as.matrix(dist(sample(1:25, 8, replace=TRUE)))
diag(plainmatrix) <- 99999  # setting diagonal to an infinite distance for
                            # pedagogical reasons (the diagonal may be left
                            # as zero)
mdm<-distancematrix(plainmatrix)
res<-nonbimatch(mdm)
}
