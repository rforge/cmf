% File cinf-add.R
\name{cinf_plot_d2m_error}
\alias{cinf_plot_d2m_error}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Produce distance to model - error plot for applicability domain studies
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Produce distance to model - error plot for applicability domain studies
}
\usage{
cinf_plot_d2m_error(y_err, d2m, mean_parts = TRUE, nparts = 3, color_parts = "red", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y_err}{
%%     ~~Describe \code{y_err} here~~
}
  \item{d2m}{
%%     ~~Describe \code{d2m} here~~
}
  \item{mean_parts}{
%%     ~~Describe \code{mean_parts} here~~
}
  \item{nparts}{
%%     ~~Describe \code{nparts} here~~
}
  \item{color_parts}{
%%     ~~Describe \code{color_parts} here~~
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
function (y_err, d2m, mean_parts = TRUE, nparts = 3, color_parts = "red", 
    ...) 
{
    npoints <- length(y_err)
    oldnum <- order(d2m)
    y_err <- y_err[oldnum]
    d2m <- d2m[oldnum]
    plot(d2m, y_err, main = "Distance to Model - Error Plot", 
        xlab = "Distance to model", ylab = "Error")
    if (mean_parts) {
        p_first <- integer(nparts)
        p_last <- integer(nparts)
        p_size <- floor(npoints/nparts)
        for (ip in 1:nparts) p_first[ip] <- (ip - 1) * p_size + 
            1
        if (nparts > 1) 
            for (ip in 1:(nparts - 1)) p_last[ip] <- p_first[ip] + 
                p_size
        p_last[nparts] <- npoints
        p_mean <- integer(nparts)
        for (ip in 1:nparts) p_mean[ip] <- mean(y_err[p_first[ip]:p_last[ip]])
        for (ip in 1:nparts) {
            x <- c(d2m[p_first[ip]], d2m[p_last[ip]])
            y <- c(p_mean[ip], p_mean[ip])
            lines(list(x = x, y = y), col = color_parts)
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
