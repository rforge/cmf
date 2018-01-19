% File cmf-kernels.R
\name{cmf_kernel_matrix_tt}
\alias{cmf_kernel_matrix_tt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Computation of the kernel (Gram) matrix for the training set
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Computation of the kernel (Gram) matrix for the training set
}
\usage{
cmf_kernel_matrix_tt(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
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
cmf_kernel_matrix_tt <- function(ft, mdb, alpha, verbose=1) {
  nmol <- length(mdb)
  gram <- matrix(0, nmol, nmol)
  atomlists <- make_atom_lists(ft, mdb)
  for (imol in 1:nmol) {
    if (verbose) {cat("."); flush.console()}
    mol <- mdb[[imol]]
    gram[imol,imol] <- cmf_kernel_al(ft, mol, mol, alpha, atomlists[[imol]], atomlists[[imol]])
  }
  if (verbose) {cat("\n"); flush.console()}
  for (imol1 in 1:(nmol-1)) {
    mol1 <- mdb[[imol1]]
	if (verbose) {cat("."); flush.console()}
    for (imol2 in (imol1+1):nmol) {
      mol2 <- mdb[[imol2]]
      gram[imol1,imol2] <- cmf_kernel_al(ft, mol1, mol2, alpha, atomlists[[imol1]], atomlists[[imol2]])
      gram[imol2,imol1] <- gram[imol1,imol2]
    }
  }
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
