\name{cv_krr}
\alias{cv_krr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
KRR with cross-validation
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
KRR with cross-validation
}
\usage{
cv_krr(nfolds, K, y, gamma)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nfolds}{
%%     ~~Describe \code{nfolds} here~~
}
  \item{K}{
%%     ~~Describe \code{K} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{gamma}{
%%     ~~Describe \code{gamma} here~~
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
function (nfolds, K, y, gamma) 
{
    y_pred_all <- double()
    y_exp_all <- double()
    new_ind <- double()
    for (ifold in 1:nfolds) {
        ind <- 1:length(y)
        ind_train <- ind[ind\%\%nfolds != ifold - 1]
        ind_test <- ind[ind\%\%nfolds == ifold - 1]
        K_train <- K[ind_train, ind_train]
        K_test <- K[ind_test, ind_train]
        y_train <- y[ind_train]
        y_test <- y[ind_test]
        m <- build_krr_model(K_train, y_train, gamma, set_b_0)
        if (is.null(m)) 
            next
        y_exp_all <- c(y_exp_all, y_test)
        y_pred <- K_test \%*\% m$a + m$b
        y_pred_all <- c(y_pred_all, y_pred)
        new_ind <- c(new_ind, ind_test)
    }
    old_ind <- double(length(y))
    for (i in 1:length(y)) old_ind[new_ind[i]] <- i
    regr <- regr_param(y_pred_all, y_exp_all)
    list(R2 = regr$R2, RMSE = regr$RMSE, y_pred_cv = y_pred_all[old_ind])
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
