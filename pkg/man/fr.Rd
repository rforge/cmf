% File cmf-seal.R
\name{fr}
\alias{fr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{fr}
\description{fr}
\usage{
fr(par_list)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{par_list}{ }
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
fr <- function(par_list) {
	xyz_a <- apply_par_list(xyz_m, par_list)
    af <- calc_af(xyz_a, xyz_t, wij)
    if (verbose) cat(sprintf("af=\%g q1=\%g q2=\%g q3=\%g q4=\%g t1=\%g t2=\%g t3=\%g\n", 
      af, par_list[1], par_list[2], par_list[3], par_list[4], par_list[5], par_list[6], 
	  par_list[7])); flush.console()
    af
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
