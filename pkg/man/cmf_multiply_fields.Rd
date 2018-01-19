% File cmf-fields.R
\name{cmf_multiply_fields}
\alias{cmf_multiply_fields}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Get product of two fields
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
cmf_multiply_fields(grid1, grid2)
Get product of two fields
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{grid1}{
%%     ~~Describe \code{x} here~~
}
  \item{grid2}{
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
cmf_multiply_fields <- function(grid1, grid2) {
  grid <- list()
  grid$ngridx <- grid1$ngridx
  grid$ngridy <- grid1$ngridy
  grid$ngridz <- grid1$ngridz
  grid$gridx <- grid1$gridx
  grid$gridy <- grid1$gridy
  grid$gridz <- grid1$gridz
  grid$val <- grid1$val * grid2$val
  grid
}  
  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
