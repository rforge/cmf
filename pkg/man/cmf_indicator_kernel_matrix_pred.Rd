% file cmf-kernels.R
\name{cmf_indicator_kernel_matrix_pred}
\alias{cmf_indicator_kernel_matrix_pred}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Indicator kernel matrix predictor
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Indicator kernel matrix predictor
}
\usage{
cmf_indicator_kernel_matrix_pred((mdb1, mdb2, alpha, syb_type, verbose=1) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mb1}{
%%     ~~Describe \code{x} here~~
}
  \item{mb2}{
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
cmf_indicator_kernel_matrix_pred <- function(mdb1, mdb2, alpha, syb_type, verbose=1) {
  nmol1 <- length(mdb1)
  nmol2 <- length(mdb2)
  gram <- matrix(0, nmol1, nmol2)
  for (imol1 in 1:nmol1) {
    mol1 <- mdb1[[imol1]]
    for (imol2 in 1:nmol2) {
	  mol2 <- mdb2[[imol2]]
	  gram[imol1,imol2] <- cmf_indicator_kernel(mol1, mol2, alpha, syb_type)
	}
    if (verbose) {cat("."); flush.console()}
  }
  if (verbose) {cat("\n"); flush.console()}
  gram
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
