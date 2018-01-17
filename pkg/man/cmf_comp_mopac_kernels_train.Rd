% File cmf-mopac.R

\name{cmf_comp_mopac_kernels_train}
\alias{cmf_comp_mopac_kernels_train}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Computes MOPAC kernel matrices for the training set and saves it to file
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Computes MOPAC kernel matrices for the training set and saves it to file
}
\usage{
cmf_comp_mopac_kernels_train
(
  train_fname = "ligands-train.mol2",                              # Training set file name
  mopac_res_fname = "ligands-mopac-res-train.RData",               # File with MOPAC results
  mopac_kernels_train_fname = "ligands-mopac-kernels-train.RData", # Computed MOPAC kernels file name
  print_comp_kernels = TRUE,                                       # Verbose computation of kernels
  
)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{train_fname}{"ligands-train.mol2"}
  \item{mopac_res_fname}{"ligands-mopac-res-train.RData"}
  \item{mopac_kernels_train_fname}{"ligands-mopac-kernels-train.RData"}
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
 
cmf_comp_mopac_kernels_train <- function
(
  train_fname = "ligands-train.mol2",                              # Training set file name
  mopac_res_fname = "ligands-mopac-res-train.RData",               # File with MOPAC results
  mopac_kernels_train_fname = "ligands-mopac-kernels-train.RData", # Computed MOPAC kernels file name
  print_comp_kernels = TRUE                                       # Verbose computation of kernels
  
)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
