\name{qom}
\Rdversion{1.3.2}
\alias{qom}
\alias{qom,data.frame-method}
\title{
Quality of Match
}
\description{
Quality of matches show how well matched pairs differ.  For each variable the average distance is generated.  Each item in a pair is assigned a group and after several iterations the quantile of these average distances is returned.
}
\usage{
qom(covariate, matches, iterations=10000, probs=NA, seed=101, all.vals=FALSE, ...)
}

\arguments{
  \item{covariate}{
A data.frame object.
}
  \item{matches}{
A data.frame object.  Contains information on how to match the covariate data set.
}
  \item{iterations}{
An integer.  Number of iterations to run, defaults to 10,000.
}
  \item{probs}{
A numeric vector.  Probabilities to pass to the quantile function.
}
  \item{seed}{
Seed provided for random-number generation.  Default value of 101.
}
  \item{all.vals}{
A logical value.  Determines if false matches should be included in comparison.  Default value of FALSE.
}
  \item{\dots}{
Additional arguments, not used at the moment.
}
}

\details{
This fuction is useful for determining the effectiveness of your weights (when generating a distance matrix).  Weighting a variable more will lower the average distance, but it could penalize the distance of the other variables.
}

\value{
data.frame with quality of match information
}

\author{
Cole Beck
}

\examples{
df <- data.frame(id=LETTERS[1:25], val1=rnorm(25), val2=rnorm(25))
df.dist <- gendistance(df, idcol=1)
df.mdm <- distancematrix(df.dist)
df.match <- nonbimatch(df.mdm)
qom(df.dist$cov, df.match$matches)
}
