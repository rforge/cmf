% File cmf-mopac.R

\name{cmf_comp_mopac_kernels_all}
\alias{cmf_comp_mopac_kernels_all}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Computes CMF kernel matrices for the combined set of molecules and mopac fields
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Computes CMF kernel matrices for the combined set of molecules and mopac fields
}
\usage{
cmf_comp_mopac_kernels_all(
all_fname = "ligands-all.mol2",
mopac_res_fname = "ligands-mopac-res-all.RData", 
mopac_kernels_all_fname = "ligands-mopac-kernels-all.RData"
)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{all_fname}{"ligands-all.mol2"}
  \item{mopac_res_fname}{"ligands-mopac-res-all.RData"}
  \item{mopac_kernels_all_fname}{"ligands-mopac-kernels-all.RData"}
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
cmf_comp_mopac_kernels_all <- function
(
  all_fname = "ligands-all.mol2",                              # The name of the file containing all molecules
  mopac_res_fname = "ligands-mopac-res-all.RData",             # File with MOPAC results
  mopac_kernels_all_fname = "ligands-mopac-kernels-all.RData", # The name of the files containing kernels for all molecules
  ...
)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
