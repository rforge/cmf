% File cmf-krr.R
\name{cmf_krr_ecv}
\alias{cmf_krr_ecv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
External n-fold cross-validation
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
External n-fold cross-validation
}
\usage{
cmf_krr_ecv(nfolds = 5, act_all_fname = "activity-all.txt", act_colnum = 2, sep = ",", kernels_all_fname = "ligands-kernels-all.RData", ecv_fname = "ligands-ecv.RData", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nfolds}{
%%     ~~Describe \code{nfolds} here~~
}
  \item{act_all_fname}{
%%     ~~Describe \code{act_all_fname} here~~
}
  \item{act_colnum}{
%%     ~~Describe \code{act_colnum} here~~
}
  \item{sep}{
%%     ~~Describe \code{sep} here~~
}
  \item{kernels_all_fname}{
%%     ~~Describe \code{kernels_all_fname} here~~
}
  \item{ecv_fname}{
%%     ~~Describe \code{ecv_fname} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
function (nfolds = 5, act_all_fname = "activity-all.txt", act_colnum = 2, 
    sep = ",", kernels_all_fname = "ligands-kernels-all.RData", 
    ecv_fname = "ligands-ecv.RData", ...) 
{
    act <- read.table(act_all_fname, header = TRUE, sep = sep)
    y <- act[, act_colnum]
    load(kernels_all_fname)
    ecv <- cmf_krr_ecv_mem(nfolds = nfolds, y = y, kernels = kernels, 
        ...)
    save(ecv, file = ecv_fname)
    ecv
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
