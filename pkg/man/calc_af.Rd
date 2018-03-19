% File cmf-seal.R
\name{calc_af}
\alias{calc_af}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{calc_af}
\description{calc_af}
\usage{
calc_af(xyz1, xyz2, wij, alpha=0.29)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xyz1}{ }
\item{xyz2}{ }
\item{wij}{ }
\item{alpha}{0.29}
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
calc_af <- function(xyz1, xyz2, wij, alpha=0.29) {
  size1 <- dim(xyz1)[2]
  size2 <- dim(xyz2)[2]
  af <- 0
  for (i in 1:size1) {
    for (j in 1:size2) {
	  rij2 <- 0
	  for (k in 1:3) {
	    rij2 <- rij2 + (xyz1[k,i]-xyz2[k,j])^2
	  }
	  af <- af - wij[i,j] * exp(-alpha * rij2)
	}
  }
  af
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
