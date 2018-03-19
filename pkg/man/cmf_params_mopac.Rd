% File cmf-mopac.R

\name{cmf_params_mopac}
\alias{cmf_params_mopac}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Computing kernels based on MOPAC results
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Computing kernels based on MOPAC results
}
\usage{
cmf_params_mopac(mdb, mopac_mdb_res)
}
%- maybe also 'usage' for other objects documented here.
%\arguments{
%  \item{mdb}{}
%  \item{mdb}{}
%}
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
\dontrun{
# Computing kernels based on MOPAC results

cmf_params_mopac <- function(mdb, mopac_mdb_res) {
  ncomp <- length(mdb)
  for (imol in 1:ncomp) {
    mopac_mol_res <- mopac_mdb_res[[imol]]
    mol <- mdb[[imol]]
        natoms <- length(mol$atoms)
#       num_heavy_atom <- integer(natoms)
#       for (ihatom in 1:length(mopac_mol_res$heavy_atom_index)) {
#         num_heavy_atom[mopac_mol_res$heavy_atom_index[ihatom]] <- ihatom
#       }
        for (iatom in 1:natoms) {
          atom <- mol$atoms[[iatom]]
          atom$mop_q <- mopac_mol_res$atom_charges[iatom]
          if (mopac_mol_res$num_heavy_atom[iatom]) {
            atom$mop_dn   <- mopac_mol_res$Dn[mopac_mol_res$num_heavy_atom[iatom]]
                atom$mop_de   <- mopac_mol_res$De[mopac_mol_res$num_heavy_atom[iatom]]
                atom$mop_pis  <- mopac_mol_res$piS[mopac_mol_res$num_heavy_atom[iatom]]
                atom$mop_homo <- mopac_mol_res$chomo[mopac_mol_res$num_heavy_atom[iatom]]
                atom$mop_lumo <- mopac_mol_res$clumo[mopac_mol_res$num_heavy_atom[iatom]]
          } else {
            atom$mop_dn   <- 0
                atom$mop_de   <- 0
                atom$mop_pis  <- 0
                atom$mop_homo <- 0
                atom$mop_lumo <- 0
          }
          mdb[[imol]]$atoms[[iatom]] <- atom
        }
  }
  mdb
}


}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
