% File cmf-kernels.R
\name{normalize_gram}
\alias{normalize_gram}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Normalization of the Gram matrix
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Normalization of the Gram matrix
}
\usage{
normalize_gram(gram)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gram}{
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as

function(gram) {
  ndim <- dim(gram)[1]
  for (irow in 1:(ndim-1)) {
    for (icol in (irow+1):ndim) {
      gram[irow,icol] <- gram[irow,icol] / (sqrt(gram[irow,irow]) * sqrt(gram[icol,icol]))
      gram[icol,irow] <- gram[irow,icol]
    }
  }
  for (i in 1:ndim) {
    gram[i,i] <- 1.0
  }
  gram
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
