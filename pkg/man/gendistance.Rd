\name{gendistance}
\Rdversion{1.3.3}
\alias{gendistance}
\alias{gendistance,data.frame-method}
\title{
Generate a Distance Matrix
}
\description{
The gendistance function creates a distance matrix from a covariates matrix.
}
\usage{
gendistance(covariate, idcol=NULL, weights=NULL, prevent=NULL, force=NULL, rankcols=NULL, missing.weight=0.1, ndiscard=0, ...)
}

\arguments{
  \item{covariate}{
A data.frame object, containing the covariates of the data set.
}
  \item{idcol}{
An integer, providing the index of the column containing row ID's.
}
  \item{weights}{
A numeric vector, the length should match the number of columns.  This value determines how much weight is given to each column when generating the distance matrix.
}
  \item{prevent}{
An integer vector, providing the index of columns that should be used to prevent matches.  When generating the distance matrix, elements that match on these columns are given a maximum distance.
}
  \item{force}{
A integer, providing the index of the column contain information used to force pairs to match.
}
  \item{rankcols}{
An integer vector, providing the index of columns that should have the rank function applied to them before generating the distance matrix.
}
  \item{missing.weight}{
A numeric value, used to generate the weight of missingness indicator columns.  Missingness indicator columns are created if there is missing data within the data set.  Defaults to 0.1.
}
  \item{ndiscard}{
An integer, providing the number of elements that should be allowed to match phantom values.  The default value is 0.
}
  \item{\dots}{
Additional arguments, not used at this time.
}
}

\details{
Given a data.frame of covariates, generate a distance matrix.  Missing values are imputed with fill.missing.  For each column with missing data, a missingness indicator column will be added.  Phantoms are fake elements that perfectly match all elements.  They can be used to discard a certain number of elements.
}

\value{
a list object with several elements
    \item{dist }{generated distance matrix}
    \item{cov }{covariate matrix used to generate distances}
    \item{ignored }{ignored columns from original covariate matrix}
    \item{weights }{weights applied to each column in covariate matrix}
    \item{prevent }{columns used to prevent matches}
    \item{mates }{index of rows that should be forced to match}
    \item{rankcols }{index of columns that should use rank}
    \item{missing.weight }{weight to apply to missingness indicator columns}
    \item{ndiscard }{number of elements that will match phantoms}
}

\author{
Cole Beck
}

\seealso{
\code{\link{distancematrix}}
}

\examples{
df <- data.frame(id=LETTERS[1:25], val1=rnorm(25), val2=rnorm(25))
df[sample(seq_len(nrow(df)), ceiling(nrow(df)*0.1)), 2] <- NA
df.dist <- gendistance(df, idcol=1, ndiscard=2)
}
