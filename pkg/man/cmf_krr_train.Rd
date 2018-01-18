% File cmf-krr.R
\name{cmf_krr_train}
\alias{cmf_krr_train}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
To build model
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
To build model
}
\usage{
cmf_krr_train(act_train_fname = "activity-train.txt", act_colnum = 2, sep = ",", kernels_train_fname = "ligands-kernels-train.RData", model_fname = "ligands-model.RData", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{act_train_fname}{"activity-train.txt"
%%     ~~Describe \code{act_train_fname} here~~
}
  \item{act_colnum}{2
%%     ~~Describe \code{act_colnum} here~~
}
  \item{sep}{","
%%     ~~Describe \code{sep} here~~
}
  \item{kernels_train_fname}{"ligands-kernels-train.RData"
%%     ~~Describe \code{kernels_train_fname} here~~
}
  \item{model_fname}{"ligands-model.RData"
%%     ~~Describe \code{model_fname} here~~
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
function (act_train_fname = "activity-train.txt", act_colnum = 2, 
    sep = ",", kernels_train_fname = "ligands-kernels-train.RData", 
    model_fname = "ligands-model.RData", ...) 
{
    act <- read.table(act_train_fname, header = TRUE, sep = sep)
    y <- act[, act_colnum]
    load(kernels_train_fname)
    model <- cmf_krr_train_mem(y = y, kernels = kernels, ...)
    save(model, file = model_fname)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
