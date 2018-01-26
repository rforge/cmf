% File cinf-regression.R
\name{regr_param_ex}
\alias{regr_param_ex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Computes extended set of regression parameters: R2ex
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Computes extended set of regression parameters: R2ex
}
\usage{
regr_param_ex(predval, expval, trainval)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{predval}{
%%     ~~Describe \code{x} here~~
}
\item{expval}{
%%     ~~Describe \code{x} here~~
}
\item{trainval}{
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
regr_param_ex <- function(predval, expval, trainval) {
  n_p <- length(predval)
  rss <- 0
  for (i in 1:n_p) {
    rss <- rss + (predval[i] - expval[i]) * (predval[i] - expval[i])
  }
  r2ex <- 1 - ((rss/(n_p-1)) / var(trainval))
  r2ex
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
