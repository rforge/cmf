% File cmf-krr.R 
\name{cmf_permute_kernels}
\alias{cmf_permute_kernels}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
To permute kernel matrix using given permutation
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
To permute kernel matrix using given permutation
}
\usage{
cmf_permute_kernels(kernels, permutation, mfields, for_pred = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{kernels}{
%%     ~~Describe \code{kernels} here~~
}
  \item{permutation}{
%%     ~~Describe \code{permutation} here~~
}
  \item{mfields}{
%%     ~~Describe \code{mfields} here~~
}
  \item{for_pred}{
%%     ~~Describe \code{for_pred} here~~
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
function (kernels, permutation, mfields, for_pred = FALSE) 
{
    alphas <- kernels$alphas
    nfields <- length(mfields)
    permuted_kernels <- list()
    permuted_kernels$alphas <- alphas
    for (f in 1:nfields) {
        field <- mfields[f]
        permuted_kernels[[field]] <- list()
        for (ialpha in 1:length(alphas)) {
            if (for_pred) 
                permuted_kernels[[field]][[ialpha]] <- kernels[[field]][[ialpha]][, 
                  permutation]
            else permuted_kernels[[field]][[ialpha]] <- kernels[[field]][[ialpha]][permutation, 
                permutation]
        }
    }
    permuted_kernels
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
