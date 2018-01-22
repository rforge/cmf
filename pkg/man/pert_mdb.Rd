% File cmf-allign.R 
\name{pert_mdb}
\alias{pert_mdb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Perturbate molecular database
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Perturbate molecular database
}
\usage{
pert_mdb(mdb, tcoef = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mdb}{
%%     ~~Describe \code{mdb} here~~
}
  \item{tcoef}{
%%     ~~Describe \code{tcoef} here~~
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
function (mdb, tcoef = 10) 
{
    ncomp <- length(mdb)
    mdb_p <- list()
    for (imol in 1:ncomp) {
        mol <- mdb[[imol]]
        mol_p <- pert_mol(mol, tcoef)
        mdb_p[[imol]] <- mol_p
    }
    mdb_p
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
