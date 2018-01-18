% File cmf-krr.R
\name{cmf_krr_ecv_mem}
\alias{cmf_krr_ecv_mem}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
External n-fold cross-validation in memory
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
External n-fold cross-validation in memory
}
\usage{
cmf_krr_ecv_mem(nfolds = 5, y, kernels, mfields = c("q", "vdw", "logp", "abra", "abrb"), print_ecv = TRUE, plot_ecv = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nfolds}{
%%     ~~Describe \code{nfolds} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{kernels}{
%%     ~~Describe \code{kernels} here~~
}
  \item{mfields}{
%%     ~~Describe \code{mfields} here~~
}
  \item{print_ecv}{
%%     ~~Describe \code{print_ecv} here~~
}
  \item{plot_ecv}{
%%     ~~Describe \code{plot_ecv} here~~
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
function (nfolds = 5, y, kernels, mfields = c("q", "vdw", "logp", 
    "abra", "abrb"), print_ecv = TRUE, plot_ecv = TRUE, ...) 
{
    y_pred_all <- double()
    y_exp_all <- double()
    new_ind <- double()
    models <- list()
    indexes <- list()
    for (ifold in 1:nfolds) {
        if (print_ecv) {
            cat(sprintf("fold = \%d of \%d\n", ifold, nfolds))
            flush.console()
        }
        ind <- 1:length(y)
        ind_train <- ind[ind\%\%nfolds != ifold - 1]
        ind_test <- ind[ind\%\%nfolds == ifold - 1]
        kernels_train <- cmf_extract_subkernels(kernels, ind_train, 
            ind_train, mfields)
        kernels_test <- cmf_extract_subkernels(kernels, ind_test, 
            ind_train, mfields)
        y_train <- y[ind_train]
        y_test <- y[ind_test]
        y_exp_all <- c(y_exp_all, y_test)
        model <- cmf_krr_train_mem(y = y_train, kernels = kernels_train, 
            mfields = mfields, ...)
        y_pred <- cmf_krr_pred_mem(model = model, kernels = kernels_train, 
            kernels_pred = kernels_test, y_exp = y_test, ...)
        y_pred_all <- c(y_pred_all, y_pred)
        new_ind <- c(new_ind, ind_test)
        models[[ifold]] <- model
        indexes[[ifold]] <- ind_train
    }
    old_ind <- double(length(y))
    for (i in 1:length(y)) old_ind[new_ind[i]] <- i
    regr <- regr_param(y_pred_all[old_ind], y_exp_all[old_ind])
    Q2ecv <- regr$R2
    RMSEecv <- regr$RMSE
    RMSEecv_pc <- regr$RMSE_pc
    if (print_ecv) {
        cat(sprintf("Q2ecv=\%.6f RMSEecv=\%.6f (\%g\%\%)\n", Q2ecv, 
            RMSEecv, RMSEecv_pc))
        flush.console()
    }
    if (plot_ecv) {
        cinf_plotxy(y_pred_all[old_ind], y_exp_all[old_ind], 
            xlab = "Predicted", ylab = "Experiment", main = "Scatter Plot for External Cross-Validation")
        abline(coef = c(0, 1))
    }
    list(Q2ecv = Q2ecv, RMSEecv = RMSEecv, RMSEecv_pc = RMSEecv_pc, 
        y_pred_ecv = y_pred_all[old_ind], y_exp = y_exp_all[old_ind], 
        models = models, indexes = indexes)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
