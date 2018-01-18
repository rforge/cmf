% File cmf-krr.R
\name{cmf_krr_pred}
\alias{cmf_krr_pred}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
To make predictions
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
To make predictiokns
}
\usage{
cmf_krr_pred(model_fname = "ligands-model.RData", kernels_train_fname = "ligands-kernels-train.RData", kernels_pred_fname = "ligands-kernels-pred.RData", act_colnum = 2, sep = ",", act_pred_fname = "activity-pred.txt", pred_fname = "ligands-pred.RData", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model_fname}{
%%     ~~Describe \code{model_fname} here~~
}
  \item{kernels_train_fname}{
%%     ~~Describe \code{kernels_train_fname} here~~
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
  \item{pred_fname}{
%%     ~~Describe \code{pred_fname} here~~
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
function (model_fname = "ligands-model.RData", kernels_train_fname = "ligands-kernels-train.RData", 
    kernels_pred_fname = "ligands-kernels-pred.RData", act_colnum = 2, 
    sep = ",", act_pred_fname = "activity-pred.txt", pred_fname = "ligands-pred.RData", 
    ...) 
{
    load(model_fname)
    load(kernels_train_fname)
    load(kernels_pred_fname)
    iprop <- act_colnum
    if (iprop > 0) {
        act <- read.table(act_pred_fname, header = TRUE, sep = sep)
        y_exp <- act[, iprop]
    }
    else {
        y_exp <- NA
    }
    y_pred <- cmf_krr_pred_mem(model = model, kernels = kernels, 
        kernels_pred = kernels_pred, y_exp = y_exp, ...)
    pred <- list(y_pred = y_pred, y_exp = y_exp)
    save(pred, file = pred_fname)
    pred
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
