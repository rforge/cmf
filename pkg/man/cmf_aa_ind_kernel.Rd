% File cmf-kernels.R
\name{cmf_aa_ind_kernel}
\alias{cmf_aa_ind_kernel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Indicator fields
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Indicator fields
}
\usage{
cmf_aa_ind_kernel(field, atom1, atom2, alpha)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{field}{
%%     ~~Describe \code{x} here~~
}
  \item{atom1}{
%%     ~~Describe \code{x} here~~
}
  \item{atom2}{
%%     ~~Describe \code{x} here~~
}
  \item{alpha}{
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
cmf_aa_ind_kernel <- function(field, atom1, atom2, alpha) {
  val <- 0.0
  if (field == atom1$syb && field == atom2$syb) {
    dist2 <- eucldist2(atom1, atom2)
    val <- exp(- alpha * dist2 / 4.0)
    val <- val * sqrt(pi^3 / alpha^3)
  }
  val
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
