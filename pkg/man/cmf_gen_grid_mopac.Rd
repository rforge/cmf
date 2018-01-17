% File cmf-mopac.R

\name{cmf_gen_grid_mopac }
\alias{cmf_gen_grid_mopac }
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Generation of grid for MOPAC molecular co-fields
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Generation of grid for MOPAC molecular co-fields
}
\usage{
cmf_gen_grid_mopac (
train_fname = "ligands-train.mol2",
train_mopac_res_fname = "ligands-mopac-res-train.RData",
kernels_fname = "ligands-kernels.RData", 
grid_fname = "ligands-grid-krr.RData",
verbose = TRUE
)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{train_fname}{"ligands-train.mol2"}
  \item{train_mopac_res_fname}{"ligands-mopac-res-train.RData"}  
  \item{kernels_fname}{ligands-kernels.RData"}
  \item{grid_fname}{"ligands-grid-krr.RData"}
  \item{verbose}{TRUE}
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
cmf_gen_grid_mopac <- function
(
  train_fname = "ligands-train.mol2",                      # training set file name
  train_mopac_res_fname = "ligands-mopac-res-train.RData", # File with MOPAC results for the training set
  kernels_fname = "ligands-kernels.RData",                 # Computed kernels file name
  model_fname = "ligands-model.RData",                     # Model file name
  grid_fname = "ligands-grid-krr.RData",                   # Grid with regression coefficients - file name 
  verbose = TRUE,                                          # Verbose output
  ...
)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
