% File cmf-add.R
\name{cinf_plot_coverage_mean_error}
\alias{cinf_plot_coverage_mean_error}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Produce coverage - mean error plot for applicability domain studies
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Produce coverage - mean error plot for applicability domain studies
}
\usage{
cinf_plot_coverage_mean_error(
  y_err,             
  d2m,               
  cmep_type = "in",  
  smoothing = TRUE,  
  ...
  )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y_err}{
%%     ~~Describe \code{x} here~~
}
  \item{d2m}{
%%     ~~Describe \code{x} here~~
}

  \item{cmep_type}{
%%     ~~Describe \code{x} here~~
}
  \item{smoothing}{
%%     ~~Describe \code{x} here~~
}
  \item{\dots}{
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
cinf_plot_coverage_mean_error <- function
(
  y_err,             # Prediction errors 
  d2m,               # Distance to model
  cmep_type = "in",  # Coverage - mean error plot type ("in", "out", "diff")
  smoothing = TRUE,  # Whether to draw smoothing curve
  ...
)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
