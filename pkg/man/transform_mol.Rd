% File cmf-allign.R
\name{transform_mol}
\alias{transform_mol}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
transform_mol
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
transform_mol
}
\usage{
transform_mol(mol, R, T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mol}{
%%     ~~Describe \code{mol} here~~
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
function (mol, R, T) 
{
    xyz <- mol2xyz(mol)
    xyz1 <- transform_xyz(xyz, R, T)
    xyz2mol(mol, xyz1)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
