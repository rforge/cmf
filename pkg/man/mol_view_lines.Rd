% File cinf-molview.R 
\name{mol_view_lines}
\alias{mol_view_lines}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
mol_view_lines
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
mol_view_lines
}
\usage{
mol_view_lines(mol)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mol}{
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
mol_view_lines <- function(mol) {
  x <- double(2)
  y <- double(2)
  z <- double(2)
  nbonds <- length(mol$bonds)
  for (ibond in 1:nbonds) {
    bond <- mol$bonds[[ibond]]
    at1 <- bond$at1
    at2 <- bond$at2
    atom1 <- mol$atoms[[at1]]
    atom2 <- mol$atoms[[at2]]
    mx <- (atom1$x + atom2$x) / 2
    my <- (atom1$y + atom2$y) / 2
    mz <- (atom1$z + atom2$z) / 2
    x[1] <- atom1$x
    x[2] <- mx
    y[1] <- atom1$y
    y[2] <- my
    z[1] <- atom1$z
    z[2] <- mz
    lines3d(x, y, z, color=PT.Color[[atom1$el]])
    x[1] <- atom2$x
    y[1] <- atom2$y
    z[1] <- atom2$z
    lines3d(x, y, z, color=PT.Color[[atom2$el]])
  }
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
