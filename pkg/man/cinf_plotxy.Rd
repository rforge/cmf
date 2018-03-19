% File cinf-plots.R
\name{cinf_plotxy}
\alias{cinf_plotxy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Scatter plot
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Scatter plot
}
\usage{
cinf_plotxy(x, y, xlab="Predicted", ylab="Experimental", margin=0.5, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
\item{y}{
%%     ~~Describe \code{x} here~~
}
\item{xlab}{"Predicted"
%%     ~~Describe \code{x} here~~
}
\item{ylab}{"Experimental"
%%     ~~Describe \code{x} here~~
}
\item{margin}{0,5
%%     ~~Describe \code{x} here~~
}
\item{\dots}{
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
cinf_plotxy <- function(x, y, xlab="Predicted", ylab="Experimental", margin=0.5, ...) {
  xmin <- min(x)
  xmax <- max(x)
  ymin <- min(y)
  ymax <- max(y)
  xymin <- min(xmin, ymin) - margin
  xymax <- max(xmax, ymax) + margin
  plot(x, y, xlim=c(xymin, xymax), ylim=c(xymin, xymax), xlab=xlab, ylab=ylab, ...)
  abline(coef=c(0,1))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
