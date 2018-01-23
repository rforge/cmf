% File cmf-atomlists.R
\name{gen_atomlists}
\alias{gen_atomlists}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Generate atom lists
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Generate atom lists
}
\usage{
gen_atomlists(mdb, field, field_family)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mdb}{
%%     ~~Describe \code{mdb} here~~
}
  \item{field}{
%%     ~~Describe \code{field} here~~
}
  \item{field_family}{
%%     ~~Describe \code{field_family} here~~
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
function (mdb, field, field_family) 
{
    atomlists <- list()
    if (field_family == "PHCH") {
        atomlists <- gen_atomlists_phch(mdb, field)
    }
    else if (field_family == "IND") {
        atomlists <- gen_atomlists_ind(mdb, field)
    }
    else if (field_family == "MOPAC") {
        atomlists <- gen_atomlists_mopac(mdb, field)
    }
    atomlists
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
