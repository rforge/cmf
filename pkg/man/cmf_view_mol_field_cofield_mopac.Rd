% File cmf-mopac.R
\name{cmf_view_mol_field_cofield_mopac}
\alias{cmf_view_mol_field_cofield_mopac}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
%%  ~~function to do ... ~~
View molecule, co-field and field
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
View molecule, co-field and field
}
\usage{

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mdb_fname}{"ligands-train.mol2"} 
  \item{train_mopac_res_fname}{"ligands-mopac-res-train.RData"}
  \item{imol}{1}
  \item{ft}{"q"}
  \item{grid_fname}{"ligands-grid-krr.RData"}
  \item{alpha_from_model}{FALSE}
  \item{model_fname}{"ligands-model.RData"} 
  \item{alpha}{0.3}
  \item{rlevel}{0.5}
  \item{alpha_g}{0.7}
  \item{draw_field}{TRUE}
  \item{draw_cofield}{TRUE} 
  \item{draw_overlap}{FALSE}
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
cmf_view_mol_field_cofield_mopac <- function
(
  mdb_fname = "ligands-train.mol2",                        # File name for molecular database 
  train_mopac_res_fname = "ligands-mopac-res-train.RData", # File with MOPAC results for the training set
  imol = 1,                                                # Molecule to visualize
  ft = "q",                                                # Field type to visualize
  grid_fname = "ligands-grid-krr.RData",                   # File name for grid with co-fields
  alpha_from_model = FALSE,                                # Take alpha from model? (TRUE/FALSE, 1/0)
  model_fname = "ligands-model.RData",                     # Model file name
  alpha = 0.3,                                             # Alpha value (if not taken from model)
  rlevel = 0.5,                                            # Isosurface level
  alpha_g = 0.7,                                           # Alpha (non-transperancy) level
  draw_field = TRUE,                                       # Whether to draw field
  draw_cofield = TRUE,                                     # Whether to draw co-field
  draw_overlap = FALSE,                                    # Whether to draw overlap between fields and co-fields
  ...
)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
