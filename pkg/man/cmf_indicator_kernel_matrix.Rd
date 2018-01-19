% File cmf-kernels.R
\name{cmf_indicator_kernel_matrix}
\alias{cmf_indicator_kernel_matrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Indicator kernel matrix
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Indicator kernel matrix
}
\usage{
cmf_indicator_kernel_matrix(mdb, alpha, syb_type, verbose) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mdb}{
%%     ~~Describe \code{x} here~~
}
  \item{alpha}{
%%     ~~Describe \code{x} here~~
}
  \item{syb_type}{
%%     ~~Describe \code{x} here~~
}
  \item{verbose}{1
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
cmf_indicator_kernel_matrix <- function(mdb, alpha, syb_type, verbose=1) {
  nmol <- length(mdb)
  gram <- matrix(0, nmol, nmol)
  for (imol in 1:nmol) {
    if (verbose) {cat("."); flush.console()}
    mol <- mdb[[imol]]
    gram[imol,imol] <- cmf_indicator_kernel(mol, mol, alpha, syb_type)
  }
  if (verbose) {cat("\n"); flush.console()}
  for (imol1 in 1:(nmol-1)) {
    mol1 <- mdb[[imol1]]
	if (verbose) {cat("."); flush.console()}
    for (imol2 in (imol1+1):nmol) {
      mol2 <- mdb[[imol2]]
      gram[imol1,imol2] <- cmf_indicator_kernel(mol1, mol2, alpha, syb_type)
      gram[imol2,imol1] <- gram[imol1,imol2]
    }
  }
  if (verbose) {cat("\n"); flush.console()}
  gram
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
