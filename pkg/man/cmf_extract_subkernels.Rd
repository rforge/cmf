% file cmf-krr.R
\name{cmf_extract_subkernels}
\alias{cmf_extract_subkernels}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
To extract subkernels
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
To extract subkernels
}
\usage{
cmf_extract_subkernels(kernels, ind_rows, ind_cols, mfields)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{kernels}{
%%     ~~Describe \code{kernels} here~~
}
  \item{ind_rows}{
%%     ~~Describe \code{ind_rows} here~~
}
  \item{ind_cols}{
%%     ~~Describe \code{ind_cols} here~~
}
  \item{mfields}{
%%     ~~Describe \code{mfields} here~~
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
function (kernels, ind_rows, ind_cols, mfields) 
{
    alphas <- kernels$alphas
    nfields <- length(mfields)
    subkernels <- list()
    subkernels$alphas <- alphas
    for (f in 1:nfields) {
        field <- mfields[f]
        subkernels[[field]] <- list()
        for (ialpha in 1:length(alphas)) {
            subkernels[[field]][[ialpha]] <- kernels[[field]][[ialpha]][ind_rows, 
                ind_cols]
        }
    }
    subkernels
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
