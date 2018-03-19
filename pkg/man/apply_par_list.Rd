% File cmf-seal.R
\name{apply_par_list}
\alias{apply_par_list}

\title{apply_par_list}
\description{apply_par_list}
\usage{
apply_par_list(xyz_m, par_list)
}
\arguments{
  \item{xyz_m}{ }
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
apply_par_list <- function(xyz_m, par_list) {
    q1 <- par_list[1]
    q2 <- par_list[2]
    q3 <- par_list[3]
    q4 <- par_list[4]
    t1 <- par_list[5]
    t2 <- par_list[6]
    t3 <- par_list[7]
    R <- q2r(q1, q2, q3, q4)
    Tr <- c(t1, t2, t3)
    xyz_a <- R %*% xyz_m + T
	xyz_a
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
