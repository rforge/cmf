% File cmf-atomlists.R
\name{gen_atomlists_generic}
\alias{gen_atomlists_generic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Generate atom lists for generic molecular field
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Generate atom lists for generic molecular field
}
\usage{
gen_atomlists_generic(mdb)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mdb}{
%%     ~~Describe \code{mdb} here~~
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
function (mdb) 
{
    nmols <- length(mdb)
    atomlists <- list()
    for (imol in 1:nmols) {
        mol <- mdb[[imol]]
        natoms <- length(mol$atoms)
        atomlists[[imol]] <- 1:natoms
    }
    atomlists
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
