% File cmf-allign.R 
\name{rmse4xyz}
\alias{rmse4xyz}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
rmse4xyz
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
rmse4xyz
}
\usage{
rmse4xyz(xyz1, xyz2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xyz1}{
%%     ~~Describe \code{xyz1} here~~
}
  \item{xyz2}{
%%     ~~Describe \code{xyz2} here~~
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
function (xyz1, xyz2) 
{
    xyz_d2 <- (xyz1 - xyz2)^2
    rss <- sum(colSums(xyz_d2))
    N <- dim(xyz1)[2]
    sqrt(rss/N)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
