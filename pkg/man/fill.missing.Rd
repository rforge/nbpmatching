\name{fill.missing}
\Rdversion{1.4.0}
\alias{fill.missing}
\alias{fill.missing,data.frame-method}
\title{
Data Imputation
}
\description{
The fill.missing function uses the \sQuote{transcan} function from the \pkg{Hmisc} package to impute values for the given data.frame.
}
\usage{
fill.missing(x, seed = 101, simplify=TRUE, idcol="id", ...)
}

\arguments{
  \item{x}{
A data.frame object.  It should have missing values.
}
  \item{seed}{                                                                                                                        
Seed provided for random-number generation.  Default value of 101.
}
  \item{simplify}{
A boolean value.  When TRUE, duplicate missingness columns will be removed.
}
  \item{idcol}{
An integer value or character string.  Indicates the column containing IDs, specified as column index or column name.  Defaults to "id", or NA, when not found.
}
  \item{\dots}{
Additional arguments, potentially passed to transcan.
}
}

\details{
The fill.missing function with fill the missing values within a data.frame with the values imputed with the transcan function.  An idcol may be specified to prevent including the use of IDs in the imputation.  In addition for every column that contains missing data, a new column will be attached to the data.frame containing an indicator of missingness.  A "1" indicates that the value was missing and has been imputed.
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
set.seed(1)
df <- data.frame(id=LETTERS[1:25], val1=rnorm(25), val2=rnorm(25))
df[sample(seq_len(nrow(df)), ceiling(nrow(df)*0.1)), 2] <- NA
df <- fill.missing(df)
}
