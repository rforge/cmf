% File cmf-coef.R
\name{cmf_coef_grid_sv}
\alias{cmf_coef_grid_sv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Computes coefficients for grid using support vectors
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Computes coefficients for grid using support vectors
}
\usage{
cmf_coef_grid_sv(mdb, a, ai, alpha, grid, field)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mdb}{
%%     ~~Describe \code{x} here~~
}
  \item{a}{
%%     ~~Describe \code{x} here~~
}
  \item{ai}{
%%     ~~Describe \code{x} here~~
}
  \item{alpha}{
%%     ~~Describe \code{x} here~~
}
  \item{grid}{
%%     ~~Describe \code{x} here~~
}
  \item{field}{
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

  cmf_coef_grid_sv <- function(mdb, a, ai, alpha, grid, field) {
  for (igridx in 1:grid$ngridx) {
  cat(sprintf("x: \%d of \%d\\n", igridx, grid$ngridx)); flush.console() 
    flush.console()
    x <- grid$gridx[igridx]
    for (igridy in 1:grid$ngridy) {
      y <- grid$gridy[igridy]
      for (igridz in 1:grid$ngridz) {
        z <- grid$gridz[igridz]
        grid$val[igridx,igridy,igridz] <- cmf_coef_xyz_sv(mdb,a,ai,alpha,x,y,z,field)
      }
    }
  }
  grid
}
  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
