% File cmf-mopac.R

\name{write_mopac_input_file}
\alias{write_mopac_input_file}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Writes MOPAC inpout file with Cartesian coordinates
}
\description{
Writes MOPAC inpout file with Cartesian coordinates
}
\usage{
write_mopac_input_file(mol, fname) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mol}{ }
  \item{fname}{ }
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
# Writes MOPAC inpout file with Cartesian coordinates
write_mopac_input_file <- function(mol, fname) {
  charge <- get_mol_charge(mol)
  natoms <- length(mol$atoms)
  of <- file(fname, "w")
  cat(sprintf("AUX LARGE CHARGE=\%g GEO_REF=\"SELF\" SUPER\n", charge), file=of)
  cat("\n", file=of)
  cat("\n", file=of)
  for (iatom in 1:natoms) {
    atom <- mol$atoms[[iatom]]
    cat(sprintf("\%-2s\%9.5f 1\%9.5f 1\%9.5f 1\n", atom$el, atom$x, atom$y, atom$z), file=of)
  }
  close(of)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
