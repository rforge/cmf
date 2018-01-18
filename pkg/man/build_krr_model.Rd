% File cmf-krr.R
\name{build_krr_model}
\alias{build_krr_model}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
To build krr model
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
To build krr model
}
\usage{
build_krr_model(K_train, y_train, gamma, set_b_0 = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{K_train}{
%%     ~~Describe \code{K_train} here~~
}
  \item{y_train}{
%%     ~~Describe \code{y_train} here~~
}
  \item{gamma}{
%%     ~~Describe \code{gamma} here~~
}
  \item{set_b_0}{
%%     ~~Describe \code{set_b_0} here~~
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
function (K_train, y_train, gamma, set_b_0 = FALSE) 
{
    ntrain <- length(y_train)
    if (set_b_0) {
        if (rcond(K_train) < 1e-15) 
            return(NULL)
        a <- solve(K_train + gamma * diag(ntrain), y_train)
        b <- 0
    }
    else {
        K1 <- K_train + gamma * diag(ntrain)
        K2 <- cbind(K1, matrix(1, nrow = ntrain))
        K3 <- rbind(K2, matrix(1, ncol = ntrain + 1))
        K3[ntrain + 1, ntrain + 1] <- 0
        if (rcond(K3) < 1e-15) 
            return(NULL)
        y0 <- c(y_train, 0)
        ab <- solve(K3, y0)
        a <- ab[1:ntrain]
        b <- ab[ntrain + 1]
    }
    list(a = a, b = b)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
