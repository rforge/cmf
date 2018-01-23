% File cmf-anal-contrib.R
\name{cmf_pred_anal}
\alias{cmf_pred_anal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Making predictions with analysis
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Making predictions with analysis
}
\usage{
cmf_pred_anal(model_fname = "ligands-model-pred.RData", kernels_pred_fname = "ligands-kernels-pred.RData", act_colnum = 2, sep = ",", act_pred_fname = "activity-pred.txt", is_train = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model_fname}{
%%     ~~Describe \code{model_fname} here~~
}
  \item{kernels_pred_fname}{
%%     ~~Describe \code{kernels_pred_fname} here~~
}
  \item{act_colnum}{
%%     ~~Describe \code{act_colnum} here~~
}
  \item{sep}{
%%     ~~Describe \code{sep} here~~
}
  \item{act_pred_fname}{
%%     ~~Describe \code{act_pred_fname} here~~
}
  \item{is_train}{
%%     ~~Describe \code{is_train} here~~
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
function (model_fname = "ligands-model-pred.RData", kernels_pred_fname = "ligands-kernels-pred.RData", 
    act_colnum = 2, sep = ",", act_pred_fname = "activity-pred.txt", 
    is_train = FALSE, ...) 
{
    iprop <- act_colnum
    load(kernels_pred_fname)
    load(model_fname)
    if (is_train) 
        kernels_pred <- kernels
    alphas_pred <- kernels_pred$alphas
    if (iprop > 0) {
        act <- read.table(act_pred_fname, header = TRUE, sep = sep)
        y_exp <- act[, iprop]
    }
    else {
        y_exp <- NA
    }
    mfields <- names(model$h)
    nfields <- length(mfields)
    K_pred <- cmf_calc_combined_kernels(kernels_pred, model$h, 
        model$alpha, alphas_pred)
    npred <- dim(K_pred)[1]
    ntrain <- dim(K_pred)[2]
    y_pred <- K_pred \%*\% model$a + model$b
    if (iprop > 0) {
        regr <- regr_param(y_pred, y_exp)
        cat(sprintf("R2=\%g RMSE=\%g\n", regr$R2, regr$RMSE))
        flush.console()
        plot(y_pred, y_exp, xlab = "Predicted", ylab = "Experiment")
        abline(coef = c(0, 1))
    }
    contrib <- array(0, c(nfields, npred, ntrain))
    for (f in 1:nfields) {
        fname <- mfields[f]
        kernels_interp <- cmf_kernels_interpolate(kernels_pred[[fname]], 
            model$alpha[[fname]], alphas_pred)
        for (p in 1:npred) {
            for (t in 1:ntrain) {
                contrib[f, p, t] <- model$h[[fname]] * model$a[t] * 
                  kernels_interp[p, t]
            }
        }
    }
    anal <- list()
    anal$contrib <- contrib
    anal$fields <- mfields
    anal$fld_contrib_av <- numeric(nfields)
    anal$fld_contrib <- array(0, c(npred, nfields))
    for (f in 1:nfields) {
        anal$fld_contrib_av[f] <- sum(contrib[f, , ])/npred
        for (p in 1:npred) {
            anal$fld_contrib[p, f] <- sum(contrib[f, p, ])
        }
    }
    anal$tp_contrib_av <- numeric(ntrain)
    anal$tp_contrib <- array(0, c(npred, ntrain))
    for (t in 1:ntrain) {
        anal$tp_contrib_av[t] <- sum(contrib[, , t])/npred
        for (p in 1:npred) {
            anal$tp_contrib[p, t] <- sum(contrib[, p, t])
        }
    }
    anal
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
