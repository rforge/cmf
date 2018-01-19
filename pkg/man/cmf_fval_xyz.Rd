% File cmf-fields.R
\name{cmf_fval_xyz}
\alias{cmf_fval_xyz}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Computes field value at point (x,y,z)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Computes field value at point (x,y,z)
}
\usage{
cmf_fval_xyz(ft, mol, alpha, x, y, z)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ft}{
%%     ~~Describe \code{x} here~~
}
  \item{mol}{
%%     ~~Describe \code{x} here~~
}
  \item{alpha}{
%%     ~~Describe \code{x} here~~
}
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{y}{
%%     ~~Describe \code{x} here~~
}
  \item{z}{
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
cmf_fval_xyz <- function(ft, mol, alpha, x, y, z) {
  fval <- 0.0
  natoms <- length(mol$atoms)
  for (iatom in 1:natoms) {
    atom <- mol$atoms[[iatom]]
    dist2 <- (x - atom$x)^2 + (y - atom$y)^2 + (z - atom$z)^2
	if (ft == "q") {
      fval <- fval + atom$pch * exp(- alpha * dist2 / 2.0)
	} else if (ft == "vdw") {
      dist2rel <- dist2 / (tripos_Rvdw[[atom$syb]])^2
      fval <- fval + tripos_Evdw[[atom$syb]] * exp(- alpha * dist2rel / 2.0)
	} else if (ft == "logp") {
      fval <- fval + atom$hydroph * exp(- alpha * dist2 / 2.0)
	} else if (ft == "abra") {
      fval <- fval + atom$abraham_a * exp(- alpha * dist2 / 2.0)
	} else if (ft == "abrb") {
      fval <- fval + atom$abraham_b * exp(- alpha * dist2 / 2.0)
    } else if (ft == "mop_q") {
      fval <- fval + atom$mop_q * exp(- alpha * dist2 / 2.0)
    } else if (ft == "mop_dn") {
      fval <- fval + atom$mop_dn * exp(- alpha * dist2 / 2.0)
    } else if (ft == "mop_de") {
      fval <- fval + atom$mop_de * exp(- alpha * dist2 / 2.0)
    } else if (ft == "mop_pis") {
      fval <- fval + atom$mop_pis * exp(- alpha * dist2 / 2.0)
    } else if (ft == "mop_homo") {
      fval <- fval + atom$mop_homo * exp(- alpha * dist2 / 2.0)
    } else if (ft == "mop_lumo") {
      fval <- fval + atom$mop_lumo * exp(- alpha * dist2 / 2.0)
	} else if (ft == "ind") {
      fval <- fval + exp(- alpha * dist2 / 2.0)
	} else if (ft \%in\% tripos_atom_types) {
	  if (ft == atom$syb) {
        fval <- fval + exp(- alpha * dist2 / 2.0)
	  }
	}
  }
  fval
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
