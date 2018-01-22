% File cmf-allign.R
\name{transform_xyz}
\alias{transform_xyz}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
transform_xyz
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
transform_xyz
}
\usage{
transform_xyz(xyz, R, T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xyz}{
%%     ~~Describe \code{xyz} here~~
}
  \item{R}{
%%     ~~Describe \code{R} here~~
}
  \item{T}{
%%     ~~Describe \code{T} here~~
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
function (xyz, R, T) 
{
    n <- dim(xyz)[2]
    ones <- double(n) + 1
    R \%*\% xyz + T \%*\% ones
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
