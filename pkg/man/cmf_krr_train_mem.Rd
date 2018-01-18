%File cmf-krr.R
\name{cmf_krr_train_mem}
\alias{cmf_krr_train_mem}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
To build model in memory
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
To build model in memory
}
\usage{
cmf_krr_train_mem(y, kernels, alpha_grid_search = TRUE, gamma_grid_search = FALSE, conic_kernel_combination = FALSE, optimize_h = FALSE, mfields = c("q", "vdw", "logp", "abra", "abrb"), set_b_0 = FALSE, print_interm_icv = TRUE, plot_interm_icv = TRUE, print_final_icv = TRUE, plot_final_icv = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{kernels}{
%%     ~~Describe \code{kernels} here~~
}
  \item{alpha_grid_search}{
%%     ~~Describe \code{alpha_grid_search} here~~
}
  \item{gamma_grid_search}{
%%     ~~Describe \code{gamma_grid_search} here~~
}
  \item{conic_kernel_combination}{
%%     ~~Describe \code{conic_kernel_combination} here~~
}
  \item{optimize_h}{
%%     ~~Describe \code{optimize_h} here~~
}
  \item{mfields}{
%%     ~~Describe \code{mfields} here~~
}
  \item{set_b_0}{
%%     ~~Describe \code{set_b_0} here~~
}
  \item{print_interm_icv}{
%%     ~~Describe \code{print_interm_icv} here~~
}
  \item{plot_interm_icv}{
%%     ~~Describe \code{plot_interm_icv} here~~
}
  \item{print_final_icv}{
%%     ~~Describe \code{print_final_icv} here~~
}
  \item{plot_final_icv}{
%%     ~~Describe \code{plot_final_icv} here~~
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
function (y, kernels, alpha_grid_search = TRUE, gamma_grid_search = FALSE, 
    conic_kernel_combination = FALSE, optimize_h = FALSE, mfields = c("q", 
        "vdw", "logp", "abra", "abrb"), set_b_0 = FALSE, print_interm_icv = TRUE, 
    plot_interm_icv = TRUE, print_final_icv = TRUE, plot_final_icv = TRUE, 
    ...) 
{
    var_y <- var(y)
    ncomp <- length(y)
    alphas <- kernels$alphas
    nalphas <- length(alphas)
    nfields <- length(mfields)
    Q2_best_of_best <- -1000
    model <- list()
    fr <- function(par_list) {
        try_current_hyper_params <- function() {
            m <- build_krr_model(Km, y, gamma, set_b_0)
            if (is.null(m)) 
                return()
            y_pred <- Km \%*\% m$a + m$b
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
        R2_best <- -1000
        RMSE_best <- -1
        Q2_best <- -1000
        RMSEcv_best <- -1
        minQ2R2_best <- -1000
        alpha_best <- -1
        gamma_best <- -1
        a_best <- NULL
        b_best <- NULL
        y_pred_best <- double()
        y_pred_cv_best <- double()
        h <- list()
        pos <- 1
        if (optimize_h) {
            if (conic_kernel_combination) {
                for (f in 1:nfields) h[[mfields[f]]] <- abs(par_list[f])
            }
            else {
                for (f in 1:nfields) h[[mfields[f]]] <- par_list[f]
            }
            pos <- pos + nfields
            if (!alpha_grid_search) {
                alpha <- par_list[pos]
                pos <- pos + 1
            }
            if (!gamma_grid_search) 
                gamma <- par_list[pos]
        }
        else {
            for (f in 1:nfields) h[[mfields[f]]] <- 1
            if (!alpha_grid_search) {
                alpha <- par_list[pos]
                pos <- pos + 1
            }
            if (!gamma_grid_search) 
                gamma <- par_list[pos]
        }
        if (alpha_grid_search) {
            for (ialpha in 1:length(alphas)) {
                alpha <- alphas[[ialpha]]
                Km <<- matrix(0, nrow = ncomp, ncol = ncomp)
                for (f in 1:nfields) {
                  Km <<- Km + h[[mfields[f]]] * kernels[[mfields[f]]][[ialpha]]
                }
                if (gamma_grid_search) {
                  for (gamma in gamma_list) {
                    try_current_hyper_params()
                  }
                }
                else {
                  try_current_hyper_params()
                }
            }
            alpha_best <- alphas[ialpha_best]
        }
        else {
            Km <<- cmf_calc_combined_kernels_1alpha(kernels, 
                h, alpha, alphas)
            if (gamma_grid_search) {
                for (gamma in gamma_list) {
                  try_current_hyper_params()
                }
            }
            else {
                try_current_hyper_params()
            }
        }
        if (Q2_best > Q2_best_of_best) {
            Q2_best_of_best <<- Q2_best
            if (print_interm_icv) {
                for (f in 1:nfields) cat(sprintf("h_\%s=\%g ", 
                  mfields[f], h[[mfields[f]]]))
                cat(sprintf("\n"))
                cat(sprintf("best: alpha=\%g gamma=\%g RMSE=\%g R2=\%g RMSEcv=\%g Q2=\%g \n", 
                  alpha_best, gamma_best, RMSE_best, R2_best, 
                  RMSEcv_best, Q2_best))
                flush.console()
            }
            if (plot_interm_icv) {
                cinf_plotxy(y_pred_cv_best, y, xlab = "Predicted", 
                  ylab = "Experiment", main = "Scatter Plot for Cross-Validation (Internal)")
                abline(coef = c(0, 1))
            }
            model$gamma <<- gamma_best
            for (f in 1:nfields) {
                model$h[[mfields[f]]] <<- h[[mfields[f]]]
                model$alpha[[mfields[f]]] <<- alpha_best
                if (alpha_best < alphas[1]) 
                  model$alpha[[mfields[f]]] <<- alphas[1]
                if (alpha_best > alphas[nalphas]) 
                  model$alpha[[mfields[f]]] <<- alphas[nalphas]
            }
            model$R2 <<- R2_best
            model$RMSE <<- RMSE_best
            model$y_pred <<- y_pred_best
            model$y_exp <<- y
            model$Q2 <<- Q2_best
            model$RMSEcv <<- RMSEcv_best
            model$y_pred_cv <<- y_pred_cv_best
            model$a <<- a_best
            model$b <<- b_best
        }
        RMSEcv_best
    }
    par_list <- list()
    if (optimize_h) 
        par_list <- c(par_list, rep(1, nfields))
    if (!alpha_grid_search) 
        par_list <- c(par_list, 0.25)
    if (!gamma_grid_search) 
        par_list <- c(par_list, 5)
    npars <- length(par_list)
    if (npars > 1) {
        res <- optim(par_list, fr)
    }
    else if (npars == 1) {
        res <- optimize(fr, c(0.01, 20))
    }
    else {
        res <- fr()
    }
    model$set_b_0 <- set_b_0
    if (print_final_icv) {
        for (f in 1:nfields) cat(sprintf("h_\%s=\%g ", mfields[f], 
            model$h[[mfields[f]]]))
        cat(sprintf("\n"))
        cat(sprintf("final: alpha=\%g gamma=\%g RMSE=\%g R2=\%g RMSEcv=\%g Q2=\%g \n", 
            model$alpha[1], model$gamma, model$RMSE, model$R2, 
            model$RMSEcv, model$Q2))
        flush.console()
    }
    if (plot_final_icv) {
        cinf_plotxy(model$y_pred_cv, y, xlab = "Predicted", ylab = "Experiment", 
            main = "Scatter Plot for Cross-Validation (Internal)")
        abline(coef = c(0, 1))
    }
    model
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
