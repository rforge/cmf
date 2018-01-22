% File cmf-coef.R 
\name{cmf_coef_xyz_sv}
\alias{cmf_coef_xyz_sv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Computes coefficient  at point (x,y,z) using support vectors
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Computes coefficient  at point (x,y,z) using support vectors
}
\usage{
cmf_coef_xyz_sv(mdb, a, ai, alpha, x, y, z, field)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mdb}{
%%     ~~Describe \code{x} here~~
}
  \item{a}{
%%     ~~Describe \code{x} here~~
}
  \item{ai}{
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
  \item{field}{
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
cmf_coef_xyz_sv <- function(mdb, a, ai, alpha, x, y, z, field) {
  coef <- 0.0
  nsv <- length(ai)
  for (isv in 1:nsv) {
    imol <- ai[isv]
	mol <- mdb[[imol]]
    natoms <- length(mol$atoms)
   for (iatom in 1:natoms) {
      atom <- mol$atoms[[iatom]]
      dist2 <- (x - atom$x)^2 + (y - atom$y)^2 + (z - atom$z)^2
	  if (field == "q") {
        coef <- coef + a[isv] * atom$pch * exp(- alpha * dist2 / 2.0)
	  } else if (field == "vdw") {
        dist2rel <- dist2 / (tripos_Rvdw[[atom$syb]])^2
        coef <- coef + a[isv] * tripos_Evdw[[atom$syb]] * exp(- alpha * dist2rel / 2.0)
	  } else if (field == "logp") {
        coef <- coef + a[isv] * atom$hydroph * exp(- alpha * dist2 / 2.0)
	  } else if (field == "abra") {
        coef <- coef + a[isv] * atom$abraham_a * exp(- alpha * dist2 / 2.0)
	  } else if (field == "abrb") {
        coef <- coef + a[isv] * atom$abraham_b * exp(- alpha * dist2 / 2.0)
	  } else if (field == "vdwr") {
        coef <- coef + a[isv] * tripos_Rvdw[[atom$syb]] * exp(- alpha * dist2 / 2.0)
	  } else if (field == "abrs") {
        coef <- coef + a[isv] * atom$abraham_s * exp(- alpha * dist2 / 2.0)
	  } else if (field == "abre") {
        coef <- coef + a[isv] * atom$abraham_e * exp(- alpha * dist2 / 2.0)
	  } else if (field == "mop_q") {
        coef <- coef + a[isv] * atom$mop_q * exp(- alpha * dist2 / 2.0)
	  } else if (field == "mop_dn") {
        coef <- coef + a[isv] * atom$mop_dn * exp(- alpha * dist2 / 2.0)
	  } else if (field == "mop_de") {
        coef <- coef + a[isv] * atom$mop_de * exp(- alpha * dist2 / 2.0)
	  } else if (field == "mop_pis") {
        coef <- coef + a[isv] * atom$mop_pis * exp(- alpha * dist2 / 2.0)
	  } else if (field == "mop_homo") {
        coef <- coef + a[isv] * atom$mop_homo * exp(- alpha * dist2 / 2.0)
	  } else if (field == "mop_lumo") {
        coef <- coef + a[isv] * atom$mop_lumo * exp(- alpha * dist2 / 2.0)
	  } else if (field == "ind") {
        coef <- coef + a[isv] * exp(- alpha * dist2 / 2.0)
	  } else if (field \%in\% tripos_atom_types) {
	    if (atom$syb == field) {
          coef <- coef + a[isv] * exp(- alpha * dist2 / 2.0)
		}
	  }
    }
  }
  coef
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
