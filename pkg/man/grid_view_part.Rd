% File cmf-gridviewlwvwl.R
\name{grid_view_part}
\alias{grid_view_part}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
grid_view_part
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
grid_view_part
}
\usage{
grid_view_part(grid, part = 0.01, alpha = 1, color_p = PT.Color[["O"]], color_n = PT.Color[["N"]], ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{grid}{
%%     ~~Describe \code{grid} here~~
}
  \item{part}{
%%     ~~Describe \code{part} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{color_p}{
%%     ~~Describe \code{color_p} here~~
}
  \item{color_n}{
%%     ~~Describe \code{color_n} here~~
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
function (grid, part = 0.01, alpha = 1, color_p = PT.Color[["O"]], 
    color_n = PT.Color[["N"]], ...) 
{
    npoints <- grid$ngridx * grid$ngridy * grid$ngridz
    %values <- rep(0, npoints)
    ivalue <- 1
    for (igridx in 1:grid$ngridx) {
        for (igridy in 1:grid$ngridy) {
            for (igridz in 1:grid$ngridz) {
                values[ivalue] <- abs(grid$val[igridx, igridy, 
                  igridz])
                ivalue <- ivalue + 1
            }
        }
    }
    svalues <- sort(values, decreasing = TRUE)
    level <- svalues[ceiling(npoints * part)]
    cat(sprintf("level=\%g\n", level))
    grid_view_level(grid, level, alpha, color_p, color_n, ...)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
