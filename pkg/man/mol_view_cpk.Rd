% File cinf-molview.R
\name{mol_view_cpk}
\alias{mol_view_cpk}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
mol_view_cpk
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
mol_view_cpk
}
\usage{
mol_view_cpk(mol, reset=TRUE, alpha=1, rfactor=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{mol}{
%%     ~~Describe \code{x} here~~
}
\item{reset}{TRUE
%%     ~~Describe \code{x} here~~
}
\item{alpha}{1
%%     ~~Describe \code{x} here~~
}
\item{rfactor}{0
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
mol_view_cpk <- function(mol, reset=TRUE, alpha=1, rfactor=0) {
  natoms <- length(mol$atoms)
  x <- double(natoms)
  y <- double(natoms)
  z <- double(natoms)
  colors <- double(natoms)
  radius <- double(natoms)
  for (i in 1:natoms) {
    atom <- mol$atoms[[i]]
    x[i] <- atom$x
    y[i] <- atom$y
    z[i] <- atom$z
    colors[i] <- PT.Color[[atom$el]]
    radius[i] <- PT.AtRad[[atom$el]] * rfactor
  }
  if (reset) {
    open3d()
    bg3d("white")
  }
  spheres3d(x, y, z, color=colors, radius=radius, alpha=alpha)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
