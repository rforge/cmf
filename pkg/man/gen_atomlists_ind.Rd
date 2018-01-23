% File cmf-atomlists.R
\name{gen_atomlists_ind}
\alias{gen_atomlists_ind}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Generate atom lists for continuous indicator field
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Generate atom lists for continuous indicator field
}
\usage{
gen_atomlists_ind(mdb, field)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mdb}{
%%     ~~Describe \code{mdb} here~~
}
  \item{field}{
%%     ~~Describe \code{field} here~~
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
function (mdb, field) 
{
    nmols <- length(mdb)
    atomlists <- list()
    for (imol in 1:nmols) {
        mol <- mdb[[imol]]
        natoms <- length(mol$atoms)
        al <- list()
        for (iatom in 1:natoms) {
            atom <- mol$atoms[[iatom]]
            if (atom$syb == field) 
                al <- c(al, iatom)
        }
        atomlists[[imol]] <- al
    }
    atomlists
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
