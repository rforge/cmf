% File cmf-mopac.R 
 
\name{read_strings}
\alias{read_strings}

\title{
read ncount string tokens
}
\description{
read ncount string tokens
}
\usage{
read_strings(count)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{count}{ }
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
 # read ncount string tokens
  read_strings <- function(count) {
	res <- character(count)
    icount <- 0
	while (icount < count) {
	  next_line()
	  tokens <- strsplit(aline, " +", perl=TRUE)[[1]]
	  ntokens <- length(tokens)
	  for (itoken in 2:ntokens) {
	    token <- tokens[itoken]
		res[icount+itoken-1] <- token
	  }
	  icount <- icount + ntokens - 1
	}
	res
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
