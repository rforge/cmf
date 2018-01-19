% File cmf-grid.R 
\name{cmf_moldim}
\alias{cmf_moldim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Returns dimensions of molecule
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Returns dimensions of molecule
}
\usage{
cmf_moldim(mol)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mol}{
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
cmf_moldim <- function(mol) {
  moldim <- list()
  natoms <- length(mol$atoms)
  if (natoms > 0) {
    atom <- mol$atoms[[1]]
    moldim$xmin <- atom$x
    moldim$xmax <- atom$x
    moldim$ymin <- atom$y
    moldim$ymax <- atom$y
    moldim$zmin <- atom$z
    moldim$zmax <- atom$z
    if (natoms > 1) {
      for (iatom in 2:natoms) {
        atom <- mol$atoms[[iatom]]
        if (atom$x < moldim$xmin) moldim$xmin <- atom$x
        if (atom$x > moldim$xmax) moldim$xmax <- atom$x
        if (atom$y < moldim$ymin) moldim$ymin <- atom$y
        if (atom$y > moldim$ymax) moldim$ymax <- atom$y
        if (atom$z < moldim$zmin) moldim$zmin <- atom$z
        if (atom$z > moldim$zmax) moldim$zmax <- atom$z
      }
    }
  }
  moldim
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
