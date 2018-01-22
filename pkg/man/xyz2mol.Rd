% File cmf-allign.R
\name{xyz2mol}
\alias{xyz2mol}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
xyz2mol
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
xyz2mol
}
\usage{
xyz2mol(mol, xyz)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mol}{
%%     ~~Describe \code{mol} here~~
}
  \item{xyz}{
%%     ~~Describe \code{xyz} here~~
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
function (mol, xyz) 
{
    mol1 <- mol
    natoms <- length(mol$atoms)
    for (iatom in 1:natoms) {
        mol1$atoms[[iatom]]$x <- xyz[1, iatom]
        mol1$atoms[[iatom]]$y <- xyz[2, iatom]
        mol1$atoms[[iatom]]$z <- xyz[3, iatom]
    }
    mol1
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
