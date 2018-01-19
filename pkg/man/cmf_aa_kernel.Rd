% File cmf-kernels.R
\name{cmf_aa_kernel}
\alias{cmf_aa_kernel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Computation of the kernel that compares fields of two atoms
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Computation of the kernel that compares fields of two atoms
}
\usage{
cmf_aa_kernel(ft, atom1, atom2, alpha))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ft}{
%%     ~~Describe \code{x} here~~
}
  \item{atom1}{
%%     ~~Describe \code{x} here~~
}
  \item{atom2}{
%%     ~~Describe \code{x} here~~
}
  \item{alpha}{
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
function(ft, atom1, atom2, alpha) {
  dist2 <- eucldist2(atom1, atom2)
  val <- 0.0
  if (ft == "q") {
    val <- atom1$pch * atom2$pch * exp(- alpha * dist2 / 4.0)
  } else if (ft == "vdwr") {
    val <- tripos_Rvdw[[atom1$syb]] * tripos_Rvdw[[atom2$syb]] * exp(- alpha * dist2 / 4.0)
  } else if (ft == "vdw") {
    w1 <- tripos_Evdw[[atom1$syb]]
    w2 <- tripos_Evdw[[atom2$syb]]
    a1 <- alpha / tripos_Rvdw[[atom1$syb]]^2
    a2 <- alpha / tripos_Rvdw[[atom2$syb]]^2
    coef <- (4 * sqrt(pi^3)) / ((a1+a2) * sqrt(2*a1+2*a2))
    expart <- exp( - (a1*a2*dist2) / (2 * (a1+a2)) )
    val <- w1 * w2 * coef * expart
  } else if (ft == "logp") {
    val <- atom1$hydroph * atom2$hydroph * exp(- alpha * dist2 / 4.0)
  } else if (ft == "abra") {
    val <- atom1$abraham_a * atom2$abraham_a * exp(- alpha * dist2 / 4.0)
  } else if (ft == "abrb") {
    val <- atom1$abraham_b * atom2$abraham_b * exp(- alpha * dist2 / 4.0)
  } else if (ft == "abrs") {
    val <- atom1$abraham_s * atom2$abraham_s * exp(- alpha * dist2 / 4.0)
  } else if (ft == "abre") {
    val <- atom1$abraham_e * atom2$abraham_e * exp(- alpha * dist2 / 4.0)
  } else if (ft == "mop_q") {
    val <- atom1$mop_q * atom2$mop_q * exp(- alpha * dist2 / 4.0)
  } else if (ft == "mop_dn") {
    val <- atom1$mop_dn * atom2$mop_dn * exp(- alpha * dist2 / 4.0)
  } else if (ft == "mop_de") {
    val <- atom1$mop_de * atom2$mop_de * exp(- alpha * dist2 / 4.0)
  } else if (ft == "mop_pis") {
    val <- atom1$mop_pis * atom2$mop_pis * exp(- alpha * dist2 / 4.0)
  } else if (ft == "mop_homo") {
    val <- atom1$mop_homo * atom2$mop_homo * exp(- alpha * dist2 / 4.0)
  } else if (ft == "mop_lumo") {
    val <- atom1$mop_lumo * atom2$mop_lumo * exp(- alpha * dist2 / 4.0)
  } else if (ft == "ind") {
    if (atom1$syb == atom2$syb) {
      val <- exp(- alpha * dist2 / 4.0)
	}
  }
  if (ft != "vdw") {
    val <- val * sqrt(pi^3 / alpha^3)
  }
  val
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
