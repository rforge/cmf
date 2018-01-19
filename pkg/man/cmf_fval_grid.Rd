% File cmf-fields.R
\name{cmf_fval_grid}
\alias{cmf_fval_grid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Computes field values for grid
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Computes field values for grid
}
\usage{
cmf_fval_grid(ft, mol, alpha, grid, verbose=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ft}{
%%     ~~Describe \code{x} here~~
}
  \item{mol}{
%%     ~~Describe \code{x} here~~
}
  \item{alpha}{
%%     ~~Describe \code{x} here~~
}
  \item{grid}{
%%     ~~Describe \code{x} here~~
}
  \item{verbose}{1
%%     ~~Describe \code{x} here~~
}
  \item{x}{
%%     ~~Describe \code{x} here~~
}
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
cmf_fval_grid <- function(ft, mol, alpha, grid, verbose=1) {
  for (igridx in 1:grid$ngridx) {
    if (verbose) {cat("."); flush.console()}
    x <- grid$gridx[igridx]
    for (igridy in 1:grid$ngridy) {
      y <- grid$gridy[igridy]
      for (igridz in 1:grid$ngridz) {
        z <- grid$gridz[igridz]
        grid$val[igridx,igridy,igridz] <- cmf_fval_xyz(ft, mol,alpha,x,y,z)
      }
    }
  }
  if (verbose) {cat("\n"); flush.console()}
  grid
}
  
  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
