% File cmf-krr.R
\name{cmf_krr_ecvr_pred_mem}
\alias{cmf_krr_ecvr_pred_mem}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
To make predictions using ecvr results in memory
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
To make predictions using ecvr results in memory
}
\usage{
cmf_krr_ecvr_pred_mem(ecvr, kernels, kernels_pred, y_exp, print_pred = TRUE, plot_pred = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ecvr}{
%%     ~~Describe \code{ecvr} here~~
}
  \item{kernels}{
%%     ~~Describe \code{kernels} here~~
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
function (ecvr, kernels, kernels_pred, y_exp, print_pred = TRUE, 
    plot_pred = TRUE, ...) 
{
    ntrain <- dim(kernels_pred[[2]][[1]])[2]
    npred <- dim(kernels_pred[[2]][[1]])[1]
    nmodels <- ecvr$nmodels
    models <- ecvr$models
    permutations <- ecvr$permutations
    indexes <- ecvr$indexes
    mfields <- ecvr$mfields
    YPred <- matrix(0, npred, nmodels)
    for (imod in 1:nmodels) {
        model <- models[[imod]]
        permutation <- permutations[[imod]]
        ind_train <- indexes[[imod]]
        kernels_perm <- cmf_permute_kernels(kernels, permutation, 
            mfields)
        kernels_pred_perm <- cmf_permute_kernels(kernels_pred, 
            permutation, mfields, for_pred = TRUE)
        kernels_perm_ind <- cmf_extract_subkernels(kernels_perm, 
            ind_train, ind_train, mfields)
        kernels_pred_perm_ind <- cmf_extract_subkernels(kernels_pred_perm, 
            1:npred, ind_train, mfields)
        y_pred <- cmf_krr_pred_mem(model = model, kernels = kernels_perm_ind, 
            kernels_pred = kernels_pred_perm_ind, y_exp = y_exp, 
            print_pred = FALSE, plot_pred = FALSE, ...)
        YPred[, imod] <- y_pred
    }
    YPredMean <- rowMeans(YPred)
    YPredSD <- double(npred)
    for (i in 1:npred) YPredSD[i] <- sd(YPred[i, ])
    if (!is.na(y_exp[1])) {
        if (print_pred) {
            dif <- y_exp - YPredMean
            cat("No.   Prediction   Experiment  Difference\n")
            for (ipred in 1:npred) {
                cat(sprintf("\%3d  \%.3f +- \%.3f   \%.3f    \%6.3f\n", 
                  ipred, YPredMean[ipred], YPredSD[ipred] * 2, 
                  y_exp[ipred], dif[ipred]))
            }
            if (npred > 1) {
                regr <- regr_param(YPredMean, y_exp)
                cat(sprintf("R2pred=\%g RMSEpred=\%g (\%g\%\%)\n", 
                  regr$R2, regr$RMSE, regr$RMSE_pc))
            }
            flush.console()
        }
        if (plot_pred) {
            cinf_plotxy(YPredMean, y_exp, xlab = "Predicted", 
                ylab = "Experiment", main = "Scatter Plot for Prediction")
            abline(coef = c(0, 1))
        }
    }
    if (plot_pred) {
        cinf_plotxy(YPredMean, y_exp, xlab = "Predicted", ylab = "Experiment", 
            main = "Scatter Plot for Prediction")
        abline(coef = c(0, 1))
    }
    list(YPred = YPred, YPredMean = YPredMean, YPredSD = YPredSD, 
        y_exp = y_exp)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
