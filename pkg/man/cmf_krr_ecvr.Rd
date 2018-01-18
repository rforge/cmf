% File cmf-krr.R
\name{cmf_krr_ecvr}
\alias{cmf_krr_ecvr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
External n-fold cross-validation with reshuffings
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
External n-fold cross-validation with reshuffings
}
\usage{
cmf_krr_ecvr(act_all_fname = "activity-all.txt", act_colnum = 2, sep = ",", kernels_all_fname = "ligands-kernels-all.RData", ecvr_fname = "ligands-ecvr.RData", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
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
  \item{ecvr_fname}{
%%     ~~Describe \code{ecvr_fname} here~~
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
function (act_all_fname = "activity-all.txt", act_colnum = 2, 
    sep = ",", kernels_all_fname = "ligands-kernels-all.RData", 
    ecvr_fname = "ligands-ecvr.RData", ...) 
{
    act <- read.table(act_all_fname, header = TRUE, sep = sep)
    y <- act[, act_colnum]
    load(kernels_all_fname)
    ecvr <- cmf_krr_ecvr_mem(y = y, kernels = kernels, ...)
    save(ecvr, file = ecvr_fname)
    ecvr
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
