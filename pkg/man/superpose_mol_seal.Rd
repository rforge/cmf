% File cmf-seal.R
\name{superpose_mol_seal}
\alias{superpose_mol_seal}
\title{superpose_mol_seal}
\description{superpose_mol_seal}
\usage{
superpose_mol_seal(mol_m, mol_t, verbose=TRUE, maxit=100, ntrials=5)
}
\arguments{
  \item{mol_m}{ }
\item{mol_t}{ }
\item{verbose}{TRUE}
\item{maxit}{100}
\item{ntrials}{5}
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
superpose_mol_seal <- function(mol_m, mol_t, verbose=TRUE, maxit=100, ntrials=5) {

  calc_mol2mol_transvec <- function(mol_m, mol_r) {
    xyz_m <- mol2xyz(mol_m)
    xyz_r <- mol2xyz(mol_r)
    mean_m <- rowMeans(xyz_m)
    mean_r <- rowMeans(xyz_r)
    mean_r - mean_m
  }
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
