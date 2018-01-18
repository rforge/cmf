% File cmf-krr.R
\name{cmf_krr_ecvr_mem}
\alias{cmf_krr_ecvr_mem}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
External n-fold cross-validation with reshuffling in memory
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
External n-fold cross-validation with reshuffling in memory
}
\usage{
cmf_krr_ecvr_mem(nreshuffles = 10, y, kernels, mfields = c("q", "vdw", "logp", "abra", "abrb"), print_ecvr = TRUE, plot_ecvr = TRUE, seed = -1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nreshuffles}{
%%     ~~Describe \code{nreshuffles} here~~
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
  \item{print_ecvr}{
%%     ~~Describe \code{print_ecvr} here~~
}
  \item{plot_ecvr}{
%%     ~~Describe \code{plot_ecvr} here~~
}
  \item{seed}{
%%     ~~Describe \code{seed} here~~
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
function (nreshuffles = 10, y, kernels, mfields = c("q", "vdw", 
    "logp", "abra", "abrb"), print_ecvr = TRUE, plot_ecvr = TRUE, 
    seed = -1, ...) 
{
    if (seed >= 0) 
        set.seed(seed)
    y_init <- y
    kernels_init <- kernels
    ncomp <- length(y)
    oldnum <- integer(10)
    YPred <- array(dim = c(ncomp, nreshuffles))
    Q2ecv_array <- double(nreshuffles)
    RMSEecv_array <- double(nreshuffles)
    permutations <- list()
    models <- list()
    indexes <- list()
    for (p in 1:nreshuffles) {
        if (print_ecvr) {
            cat(sprintf("reshuffle \%d of \%d\n", p, nreshuffles))
            flush.console()
        }
        permutation <- sample(1:ncomp, ncomp)
        for (i in 1:ncomp) oldnum[permutation[i]] <- i
        y_perm <- y_init[permutation]
        kernels_perm <- cmf_permute_kernels(kernels_init, permutation, 
            mfields)
        res_perm <- cmf_krr_ecv_mem(y = y_perm, kernels = kernels_perm, 
            mfields = mfields, ...)
        YPred[, p] <- res_perm$y_pred_ecv[oldnum]
        Q2ecv_array[p] <- res_perm$Q2ecv
        RMSEecv_array[p] <- res_perm$RMSEecv
        models <- c(models, res_perm$models)
        indexes <- c(indexes, res_perm$indexes)
        last <- length(permutations)
        for (ifold in 1:nfolds) permutations[[last + ifold]] <- permutation
    }
    y_pred_mean <- rowMeans(YPred)
    y_pred_sd <- double(ncomp)
    for (i in 1:ncomp) y_pred_sd[i] <- sd(YPred[i, ])
    Q2ecv_mean <- mean(Q2ecv_array)
    Q2ecv_sd <- sd(Q2ecv_array)
    RMSEecv_mean <- mean(RMSEecv_array)
    RMSEecv_sd <- sd(RMSEecv_array)
    regr <- regr_param(y_pred_mean, y)
    Q2ecv_aggr <- regr$R2
    RMSEecv_aggr <- regr$RMSE
    RMSEecv_aggr_pc <- regr$RMSE_pc
    nmodels <- length(models)
    if (print_ecvr) {
        cat(sprintf("Q2ecv_aggr=\%.6f RMSEecv_aggr=\%.6f (\%g\%\%)\n", 
            Q2ecv_aggr, RMSEecv_aggr, RMSEecv_aggr_pc))
        cat(sprintf("Q2ecv_mean=\%.6f Q2ecv_sd=\%.6f RMSEecv_mean=\%.6f RMSEecv_sd=\%.6f\n", 
            Q2ecv_mean, Q2ecv_sd, RMSEecv_mean, RMSEecv_sd))
        flush.console()
    }
    if (plot_ecvr) {
        cinf_plotxy(y_pred_mean, y, xlab = "Predicted", ylab = "Experiment", 
            main = "Scatter Plot for External Cross-Validations with Reshuffles")
        abline(coef = c(0, 1))
    }
    list(Q2ecv_aggr = Q2ecv_aggr, RMSEecv_aggr = RMSEecv_aggr, 
        RMSEecv_aggr_pc = RMSEecv_aggr_pc, YPred = YPred, y_exp = y, 
        y_pred_mean = y_pred_mean, y_pred_sd = y_pred_sd, Q2ecv_array = Q2ecv_array, 
        Q2ecv_mean = Q2ecv_mean, Q2ecv_sd = Q2ecv_sd, RMSEecv_array = RMSEecv_array, 
        RMSEecv_mean = RMSEecv_mean, RMSEecv_sd = RMSEecv_sd, 
        nmodels = nmodels, permutations = permutations, models = models, 
        indexes = indexes, mfields = mfields)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
