% File cmf-seal.R
\name{align_mdb_seal}
\alias{align_mdb_seal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Aligns molecular database mdb using template mol_t and algorithm SEAL}
\description{
Aligns molecular database mdb using template mol_t and algorithm SEAL}
\usage{
align_mdb_seal(mdb, mol_t, verbose=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mdb}{ }
  \item{mol_t}{ }
  \item{verbose}{ }

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

# The function is currently defined as
 function(mdb, mol_t, verbose=TRUE) {
  ncomp <- length(mdb)
  mdb_a <- list()
  for (imol in 1:ncomp) {
    mol_m <- mdb[[imol]]
	cat(sprintf("imol=\%d", imol))
	res <- superpose_mol_seal(mol_m, mol_t, verbose=FALSE)
	cat(sprintf(" af=\%g\n", res$af))
	mdb_a[[imol]] <- res$mol
  }
  mdb_a
}
  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
