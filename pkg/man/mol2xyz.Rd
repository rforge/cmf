% File cmf-allign.R
\name{mol2xyz}
\alias{mol2xyz}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
mol2xyz
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
mol2xyz
}
\usage{
mol2xyz(mol)
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
mol2xyz <- function(mol) {
  natoms <- length(mol$atoms)
  xyz <- matrix(0, 3, natoms)
  for (iatom in 1:natoms) {
    atom <- mol$atoms[[iatom]]
	xyz[1, iatom] <- atom$x
	xyz[2, iatom] <- atom$y
	xyz[3, iatom] <- atom$z
  }
  xyz
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
