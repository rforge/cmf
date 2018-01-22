% File cmf-allign.R
\name{rmse4mol}
\alias{rmse4mol}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
rmse4mol
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
rmse4mol
}
\usage{
rmse4mol(mol1, mol2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mol1}{
%%     ~~Describe \code{mol1} here~~
}
  \item{mol2}{
%%     ~~Describe \code{mol2} here~~
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
function (mol1, mol2) 
{
    xyz1 <- mol2xyz(mol1)
    xyz2 <- mol2xyz(mol2)
    rmse4xyz(xyz1, xyz2)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
