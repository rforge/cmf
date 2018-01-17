% File cmf-mopac.R

\name{cmf_comp_mopac_kernels_pred}
\alias{cmf_comp_mopac_kernels_pred}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Computes MOPAC kernel matrices for prediction and saves to file
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Computes MOPAC kernel matrices for prediction and saves to file
}
\usage{
cmf_comp_mopac_kernels_pred(
  train_fname = "ligands-train.mol2",                            # Training set file name
  train_mopac_res_fname = "ligands-mopac-res-train.RData",       # File with MOPAC results for the training set
  pred_fname = "ligands-pred.mol2",                              # Prediction set file name
  pred_mopac_res_fname = "ligands-mopac-res-pred.RData",         # File with MOPAC results for the test set
  mopac_kernels_pred_fname = "ligands-mopac-kernels-pred.RData", # Computed kernels file name
  print_comp_kernels = TRUE                                     # Verbose computation of kernels
)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ train_fname}{"ligands-train.mol2"}
  \item{train_mopac_res_fname}{"ligands-mopac-res-train.RData"}
  \item{pred_fname}{"ligands-pred.mol2"}
  \item{pred_mopac_res_fname}{"ligands-mopac-res-pred.RData"}
  \item{mopac_kernels_pred_fname}{"ligands-mopac-kernels-pred.RData",}
  \item{print_comp_kernels}{TRUE}
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
cmf_comp_mopac_kernels_pred <- function
(
  train_fname = "ligands-train.mol2",                            # Training set file name
  train_mopac_res_fname = "ligands-mopac-res-train.RData",       # File with MOPAC results for the training set
  pred_fname = "ligands-pred.mol2",                              # Prediction set file name
  pred_mopac_res_fname = "ligands-mopac-res-pred.RData",         # File with MOPAC results for the test set
  mopac_kernels_pred_fname = "ligands-mopac-kernels-pred.RData", # Computed kernels file name
  print_comp_kernels = TRUE,                                     # Verbose computation of kernels
  ...
)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
