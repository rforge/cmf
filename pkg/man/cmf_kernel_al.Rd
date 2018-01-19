% File cmf-kernels.R
\name{cmf_kernel_al}
\alias{cmf_kernel_al}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Computation of the kernel that compares fields of two molecules
with specified atom lists
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Computation of the kernel that compares fields of two molecules
with specified atom lists
}
\usage{
cmf_kernel_al(ft, mol1, mol2, alpha, atomlist1, atomlist2) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ft}{
%%     ~~Describe \code{x} here~~
}
  \item{mol1}{
%%     ~~Describe \code{x} here~~
}
  \item{mol2}{
%%     ~~Describe \code{x} here~~
}
  \item{alpha}{
%%     ~~Describe \code{x} here~~
}
  \item{aromlist1}{
%%     ~~Describe \code{x} here~~
}
  \item{atomlist2}{
%%     ~~Describe \code{x} here~~
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
cmf_kernel_al <- function(ft, mol1, mol2, alpha, atomlist1, atomlist2) {
  res <- 0.0
  for (iatom1 in atomlist1) {
    atom1 <- mol1$atoms[[iatom1]]
    for (iatom2 in atomlist2) {
      atom2 <- mol2$atoms[[iatom2]]
	  res <- res + cmf_aa_kernel(ft, atom1, atom2, alpha)
    }
  }
  res
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
