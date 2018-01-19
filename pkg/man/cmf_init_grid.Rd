% File cmf-grod.R 
\name{cmf_init_grid}
\alias{cmf_init_grid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Initializes grid for a given step and margin around molecules from database
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Initializes grid for a given step and margin around molecules from database
}
\usage{
cmf_init_grid(mdb, step=1.0, margin=2.0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mdb}{
%%     ~~Describe \code{x} here~~
}
  \item{step}{1.0
%%     ~~Describe \code{x} here~~
}
  \item{margin}{2.0
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
cmf_init_grid <- function(mdb, step=1.0, margin=2.0) {
  grid <- list()
  nmols <- length(mdb)
  if (nmols > 0) {
    mdbdim <- cmf_moldim(mdb[[1]])
    if (nmols > 1) {
      for (imol in 2:nmols) {
        mol <- mdb[[imol]]
        moldim <- cmf_moldim(mol)
        if (moldim$xmin < mdbdim$xmin) mdbdim$xmin <- moldim$xmin
        if (moldim$xmax > mdbdim$xmax) mdbdim$xmax <- moldim$xmax
        if (moldim$ymin < mdbdim$ymin) mdbdim$ymin <- moldim$ymin
        if (moldim$ymax > mdbdim$ymax) mdbdim$ymax <- moldim$ymax
        if (moldim$zmin < mdbdim$zmin) mdbdim$zmin <- moldim$zmin
        if (moldim$zmax > mdbdim$zmax) mdbdim$zmax <- moldim$zmax
      }
    }
    mdbdim$xmin <- mdbdim$xmin - margin
    mdbdim$xmax <- mdbdim$xmax + margin
    mdbdim$ymin <- mdbdim$ymin - margin
    mdbdim$ymax <- mdbdim$ymax + margin
    mdbdim$zmin <- mdbdim$zmin - margin
    mdbdim$zmax <- mdbdim$zmax + margin

    sizex <- mdbdim$xmax - mdbdim$xmin
    grid$ngridx <- ceiling(sizex / step) + 1
    grid$gridx <- numeric(grid$ngridx)
    grid$gridx[1] <- mdbdim$xmin
    grid$gridx[grid$ngridx] <- mdbdim$xmax
    if (grid$ngridx > 2) {
      for (igrid in 1:(grid$ngridx-2)) {
        grid$gridx[igrid+1] <- mdbdim$xmin + igrid * step
      }
    }

    sizey <- mdbdim$ymax - mdbdim$ymin
    grid$ngridy <- ceiling(sizey / step) + 1
    grid$gridy <- numeric(grid$ngridy)
    grid$gridy[1] <- mdbdim$ymin
    grid$gridy[grid$ngridy] <- mdbdim$ymax
    if (grid$ngridy > 2) {
      for (igrid in 1:(grid$ngridy-2)) {
        grid$gridy[igrid+1] <- mdbdim$ymin + igrid * step
      }
    }

    sizez <- mdbdim$zmax - mdbdim$zmin
    grid$ngridz <- ceiling(sizez / step) + 1
    grid$gridz <- numeric(grid$ngridz)
    grid$gridz[1] <- mdbdim$zmin
    grid$gridz[grid$ngridz] <- mdbdim$zmax
    if (grid$ngridz > 2) {
      for (igrid in 1:(grid$ngridz-2)) {
        grid$gridz[igrid+1] <- mdbdim$zmin + igrid * step
      }
    }

    grid$val <- array(0.0, c(grid$ngridx, grid$ngridy, grid$ngridz))
  }
  grid
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
