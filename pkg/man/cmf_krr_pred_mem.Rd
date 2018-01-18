% File cmf-krr.R 
\name{cmf_krr_pred_mem}
\alias{cmf_krr_pred_mem}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
To meke predictions in memory
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
To meke predictions in memory
}
\usage{
cmf_krr_pred_mem(model, kernels_pred, y_exp, print_pred = TRUE, plot_pred = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
%%     ~~Describe \code{model} here~~
}
  \item{kernels_pred}{
%%     ~~Describe \code{kernels_pred} here~~
}
  \item{y_exp}{
%%     ~~Describe \code{y_exp} here~~
}
  \item{print_pred}{
%%     ~~Describe \code{print_pred} here~~
}
  \item{plot_pred}{
%%     ~~Describe \code{plot_pred} here~~
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
function (model, kernels_pred, y_exp, print_pred = TRUE, plot_pred = TRUE, 
    ...) 
{
    alphas_pred <- kernels_pred$alphas
    y_train <- model$y_exp
    K_pred <- cmf_calc_combined_kernels(kernels_pred, model$h, 
        model$alpha, alphas_pred)
    y_pred <- K_pred \%*\% model$a + model$b
    if (!is.na(y_exp[1])) {
        if (print_pred) {
            regr <- regr_param(y_pred, y_exp)
            r2ex <- regr_param_ex(y_pred, y_exp, model$y_exp)
            cat(sprintf("R2pred=\%g RMSEpred=\%g (\%g\%\%) R2pred_ex=\%g\n", 
                regr$R2, regr$RMSE, regr$RMSE_pc, r2ex))
            flush.console()
        }
        if (plot_pred) {
            cinf_plotxy(y_pred, y_exp, xlab = "Predicted", ylab = "Experiment", 
                main = "Scatter Plot for External Prediction")
            abline(coef = c(0, 1))
        }
    }
    y_pred
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
