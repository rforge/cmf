% Filw cinf-regression.R
\name{mlr}
\alias{mlr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Build multiple linear regression model
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Build multiple linear regression model
}
\usage{
mlr(df, stepwise = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{
%%     ~~Describe \code{df} here~~
}
  \item{stepwise}{
%%     ~~Describe \code{stepwise} here~~
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
function (df, stepwise = TRUE) 
{
    names <- names(df)
    ncols <- length(names)
    propname <- names[ncols]
    formula <- paste(propname, "~ .")
    model <- lm(formula, df)
    if (stepwise) 
        model <- step(model)
    plot_regr(fitted.values(model), fitted.values(model) + residuals(model))
    model
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
