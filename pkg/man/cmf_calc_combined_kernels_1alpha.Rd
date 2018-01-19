% file cmf-kernels.R
\name{cmf_calc_combined_kernels_1alpha}
\alias{cmf_calc_combined_kernels_1alpha}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Kernel combination with linear interpolation for fields with the same value of alpha
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Kernel combination with linear interpolation for fields with the same value of alpha
}
\usage{
cmf_calc_combined_kernels_1alpha(kernels, h, alpha, alphas) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{kernels}{
%%     ~~Describe \code{x} here~~
}
  \item{h}{
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
cmf_calc_combined_kernels_1alpha <- function(kernels, h, alpha, alphas) 
{
  mfields <- names(h)
  nfields <- length(mfields)
  alpha_f <- list()
  for (f in 1:nfields) {
    alpha_f[[ mfields[f] ]] <- alpha
  }
  cmf_calc_combined_kernels(kernels, h, alpha_f, alphas)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
