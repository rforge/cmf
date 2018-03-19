% File cmf-triposff.R
\name{get_syb_types_list}
\alias{get_syb_types_list}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extracts list of sybyl types of atoms
}
\description{
Extracts list of sybyl types of atoms
}
\usage{
get_syb_types_list()
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mdb}{ }
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
# Extracts list of sybyl types of atoms
get_syb_types_list <- function(mdb)
{
  syb_types <- NULL
  nmol <- length(mdb)
  for (imol in 1:nmol) {
    mol <- mdb[[imol]]
    natoms <- length(mol$atoms)
	for (iatom in 1:natoms) {
	  atom <- mol$atoms[[iatom]]
	  if (!(atom$syb \%in\% syb_types)) {
	    syb_types <- c(syb_types, atom$syb)
	  }
	}
  }
  syb_types
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line


