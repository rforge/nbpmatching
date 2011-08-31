\name{fill.missing}
\Rdversion{1.3}
\alias{fill.missing}
\alias{fill.missing,data.frame-method}
\title{
Data Imputation
}
\description{
The fill.missing function uses the transcan function from the Hmisc package to impute values for the given data.frame.
}
\usage{
fill.missing(x, seed = 101, ...)
}

\arguments{
  \item{x}{
A data.frame object.  It should have missing values.
}
  \item{seed}{                                                                                                                        
Seed provided for random-number generation.  Default value of 101.
}
  \item{\dots}{
Additional arguments, potentially passed to transcan.
}
}

\details{
The fill.missing function with fill the missing values within a data.frame with the values imputed with the transcan function.  In addition for every column that contains missing data, a new column will be attached to the data.frame containing an indicator of missingness.  A "1" indicates that the value was missing and has been imputed.
}

\value{
data.frame with imputed values
}

\author{
Cole Beck
}

\seealso{
\code{\link{transcan}}
}

\examples{
df <- data.frame(id=LETTERS[1:25], val1=rnorm(25), val2=rnorm(25))
df[sample(seq_len(nrow(df)), ceiling(nrow(df)*0.1)), 2] <- NA
df <- fill.missing(df)
}
