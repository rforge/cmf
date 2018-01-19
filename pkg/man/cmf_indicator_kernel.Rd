% File cmf-kernels.R
\name{cmf_indicator_kernel}
\alias{cmf_indicator_kernel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Indicator kernel
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Indicator kernel
}
\usage{
cmf_indicator_kernel(mol1, mol2, alpha, syb_type) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mol1}{
%%     ~~Describe \code{x} here~~
}
  \item{mol2}{
%%     ~~Describe \code{x} here~~
}
  \item{alpha}{
%%     ~~Describe \code{x} here~~
}
  \item{syb_type}{
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
cmf_indicator_kernel <- function(mol1, mol2, alpha, syb_type) {
  res <- 0.0
  natoms1 <- length(mol1$atoms)
  natoms2 <- length(mol2$atoms)
  for (iatom1 in 1:natoms1) {
    atom1 <- mol1$atoms[[iatom1]]
    if(atom1$syb != syb_type) next
    for (iatom2 in 1:natoms2) {
      atom2 <- mol2$atoms[[iatom2]]
      if(atom2$syb != syb_type) next
      dist2 <- eucldist2(atom1, atom2)
      res <- res + exp(- alpha * dist2 / 4.0)
    }
  }
  coef <- sqrt(pi^3 / alpha^3)
  res <- coef * res
  res
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
