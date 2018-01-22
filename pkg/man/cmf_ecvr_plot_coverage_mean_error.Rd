% File cmf-add.R
\name{cmf_ecvr_plot_coverage_mean_error}
\alias{cmf_ecvr_plot_coverage_mean_error}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Producing plot "coverage - mean error" for applicability domain studies
Prediction variance is taken as distance to model
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Producing plot "coverage - mean error" for applicability domain studies
Prediction variance is taken as distance to model
}
\usage{
cmf_ecvr_plot_coverage_mean_error(ecvr_fname = "ligands-ecvr.RData", cmep_type = "in", smoothing = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ecvr_fname}{
%%     ~~Describe \code{ecvr_fname} here~~
}
  \item{cmep_type}{
%%     ~~Describe \code{cmep_type} here~~
}
  \item{smoothing}{
%%     ~~Describe \code{smoothing} here~~
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
function (ecvr_fname = "ligands-ecvr.RData", cmep_type = "in", 
    smoothing = TRUE, ...) 
{
    load(ecvr_fname)
    y_err <- abs(ecvr$y_exp - ecvr$y_pred_mean)
    d2m <- ecvr$y_pred_sd
    cinf_plot_coverage_mean_error(y_err, d2m, cmep_type = cmep_type, 
        smoothing = smoothing, ...)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
