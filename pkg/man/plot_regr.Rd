% File cinf-regression.R
\name{plot_regr}
\alias{plot_regr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Plots regression
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Plots regression
}
\usage{
plot_regr(x, y, xlab = "Predicted", ylab = "Experimental", rparams = TRUE)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{xlab}{
%%     ~~Describe \code{xlab} here~~
}
  \item{ylab}{
%%     ~~Describe \code{ylab} here~~
}
  \item{rparams}{
%%     ~~Describe \code{rparams} here~~
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
function (x, y, xlab = "Predicted", ylab = "Experimental", rparams = TRUE) 
{
    plot(x, y, xlab = xlab, ylab = ylab)
    abline(coef = c(0, 1))
    if (rparams) {
        rp <- regr_param(x, y)
        text <- sprintf("R2=\%6.4f RMSE=\%.3g MAE=\%.3g", rp$R2, 
            rp$RMSE, rp$MAE)
        title(sub = text)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
