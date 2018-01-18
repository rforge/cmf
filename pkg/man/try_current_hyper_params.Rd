% File cmf-krr.R
\name{try_current_hyper_params}
\alias{try_current_hyper_params}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Try current set of hyperparameters
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Try current set of hyperparameters
}
\usage{
try_current_hyper_params()
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{}{
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
function () 
{
  m <- build_krr_model(Km, y, gamma, set_b_0)
      if (is.null(m)) return()
	  y_pred <- Km %*% m$a + m$b
      regr <- regr_param(y_pred, y)
      RMSE <- regr$RMSE
      R2 <- regr$R2
      cv <- cv_krr(10, Km, y, gamma)
      RMSEcv <- cv$RMSE
      Q2 <- cv$R2
      y_pred_cv <- cv$y_pred_cv
      minQ2R2 <- min(Q2, R2)
      if (minQ2R2 > minQ2R2_best) {
        minQ2R2_best <<- minQ2R2
        RMSE_best <<- RMSE
        R2_best <<- R2
        RMSEcv_best <<- RMSEcv
        Q2_best <<- Q2
		if (alpha_grid_search) {
          ialpha_best <<- ialpha
		}
		alpha_best <<- alpha
        gamma_best <<- gamma
        a_best <<- m$a
        b_best <<- m$b
        y_pred_best <<- y_pred
        y_pred_cv_best <<- y_pred_cv
	  }
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
