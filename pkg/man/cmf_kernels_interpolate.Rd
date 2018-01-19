% File cmf-kernels.R
\name{cmf_kernels_interpolate}
\alias{cmf_kernels_interpolate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Linear interpolation of kernel values
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Linear interpolation of kernel values
}
\usage{
cmf_kernels_interpolate(kernels_a, alpha, alphas)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{kernels_a}{
%%     ~~Describe \code{x} here~~
}
  \item{alpha}{
%%     ~~Describe \code{x} here~~
}
  \item{alphas}{
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
cmf_kernels_interpolate <- function(kernels_a, alpha, alphas)
{
  nalphas <- length(alphas)
  if (alpha >= alphas[1] && alpha <= alphas[nalphas]) {
    iamin <- 1
    iamax <- 2
    for (ia in 1:(nalphas-1)) {
      amin <- alphas[ia]
      amax <- alphas[ia+1]
      if ((alpha>=amin) && (alpha<=amax)) {
        iamin <- ia
        iamax <- ia + 1			
        break
      }
    }
    c1 <- (amax - alpha) / (amax - amin)
    c2 <- (alpha - amin) / (amax - amin)
    res <- c1 * kernels_a[[iamin]] + c2 * kernels_a[[iamax]]
  } else if (alpha < alphas[1]) {
    res <- kernels_a[[1]]
  } else if (alpha > alphas[nalphas]) {
    res <- kernels_a[[nalphas]]
  }  
  return(res)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
