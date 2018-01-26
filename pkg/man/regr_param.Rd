% File cinf-regression.R
\name{regr_param}
\alias{regr_param}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Computes regression parameters: R2, RMSE, RMSE_pc, MAE
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Computes regression parameters: R2, RMSE, RMSE_pc, MAE
}
\usage{
regr_param(predval, expval)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{predval}{
%%     ~~Describe \code{predval} here~~
}
  \item{expval}{
%%     ~~Describe \code{expval} here~~
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
function (predval, expval) 
{
    n <- length(predval)
    expval_mean <- mean(expval)
    ss <- var(expval) * (n - 1)
    rss <- 0
    rsa <- 0
    for (i in 1:n) {
        rss <- rss + (predval[i] - expval[i]) * (predval[i] - 
            expval[i])
        rsa <- rsa + abs(predval[i] - expval[i])
    }
    r2 <- (ss - rss)/ss
    rmse <- sqrt(rss/n)
    mae <- rsa/n
    list(R2 = r2, RMSE = rmse, MAE = mae)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
