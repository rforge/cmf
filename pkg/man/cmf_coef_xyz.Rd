% File cmf-coef.R 
\name{cmf_coef_xyz}
\alias{cmf_coef_xyz}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Computes coefficient  at point (x,y,z)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Computes coefficient  at point (x,y,z)
}
\usage{
cmf_coef_xyz(mdb, a, alpha, x, y, z, field)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mdb}{
%%     ~~Describe \code{x} here~~
}
\item{a}{
%%     ~~Describe \code{x} here~~
}
\item{alpha}{
%%     ~~Describe \code{x} here~~
}
\item{x}{
%%     ~~Describe \code{x} here~~
}
\item{y}{
%%     ~~Describe \code{x} here~~
}
\item{z}{
%%     ~~Describe \code{x} here~~
}
\item{field}{
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
cmf_coef_xyz <- function(mdb, a, alpha, x, y, z, field) {
  nmols <- length(a)
  ai <- 1:nmols
  cmf_coef_xyz_cv(mdb, a, ai, alpha, x, y, z, field)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
