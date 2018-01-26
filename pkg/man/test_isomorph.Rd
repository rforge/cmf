% File cinf-isomorph.R
\name{test_isomorph}
\alias{test_isomorph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Finds mol2 in mol1
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Finds mol2 in mol1
}
\usage{
test_isomorph()
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{}{
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
test_isomorph <- function() {
  # mol1
  mol1_lab <- c("C", "C", "C", "C", "C", "N", "O")
  mol1_ct <- matrix(0, nrow=7, ncol=7)
  mol1_ct[1,2] <- mol1_ct[2,1] <- 1
  mol1_ct[2,3] <- mol1_ct[3,2] <- 1
  mol1_ct[3,4] <- mol1_ct[4,3] <- 1
  mol1_ct[4,5] <- mol1_ct[5,4] <- 1
  mol1_ct[5,6] <- mol1_ct[6,5] <- 2
  mol1_ct[1,6] <- mol1_ct[6,1] <- 1
  mol1_ct[4,7] <- mol1_ct[7,4] <- 2

  # mol2
  mol2_lab <- c("O", "C", "C", "N")
  mol2_ct <- matrix(0, nrow=4, ncol=4)
  mol2_ct[1,2] <- mol2_ct[2,1] <- 2
  mol2_ct[2,3] <- mol2_ct[3,2] <- 1
  mol2_ct[3,4] <- mol2_ct[4,3] <- 2

  isomorph <- find_substr_isomorph(mol2_lab, mol2_ct, mol1_lab, mol1_ct)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
