% File cmf-kernels.R
\name{make_atom_lists}
\alias{make_atom_lists}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Creation of atom lists
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Creation of atom lists
}
\usage{

make_atom_lists(ft, mdb)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ft}{
%%     ~~Describe \code{x} here~~
}
  \item{mdb}{
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
make_atom_lists <- function(ft, mdb) {
  atomlists <- list()
  ncomp <- length(mdb)
  for (imol in 1:ncomp) {
    mol <- mdb[[imol]]
	natoms <- length(mol$atoms)
	atomlist <- 1:natoms
	if (ft %in% heavy_atom_fields) {
	  isheavy <- logical(natoms)
	  for (i in 1:natoms) isheavy[i] <- mol$atoms[[i]]$el != "H"
      atomlist <- atomlist[isheavy]	  
	}

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
