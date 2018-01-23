% File cmf-comp-kernels.R 
\name{comp_ind_kernels_all}
\alias{comp_ind_kernels_all}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Computes CMF kernel matrices for the combined set of molecules and indicator fields
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Computes CMF kernel matrices for the combined set of molecules and indicator fields
}
\usage{
comp_ind_kernels_all(all_fname = "ligands-all.mol2", ind_kernels_all_fname = "ligands-ind-kernels-all.RData", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{all_fname}{
%%     ~~Describe \code{all_fname} here~~
}
  \item{ind_kernels_all_fname}{
%%     ~~Describe \code{ind_kernels_all_fname} here~~
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
function (all_fname = "ligands-all.mol2", ind_kernels_all_fname = "ligands-ind-kernels-all.RData", 
    ...) 
{
    comp_ind_kernels_train(train_fname = all_fname, ind_kernels_train_fname = ind_kernels_all_fname, 
        ...)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
