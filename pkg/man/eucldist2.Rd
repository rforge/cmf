% File cmf-kernels.R
\name{eucldist2}
\alias{eucldist2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Calculation of the square of the Euclidean distance between two atoms
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Calculation of the square of the Euclidean distance between two atoms
}
\usage{
eucldist2(atom1, atom2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{atom1}{
%%     ~~Describe \code{x} here~~
}
  \item{atom2}{
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
function(atom1, atom2) {
  dist2 <- (atom1$x - atom2$x)^2 + (atom1$y - atom2$y)^2 + (atom1$z - atom2$z)^2
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
