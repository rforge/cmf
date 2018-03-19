% File cmf-view-mol-fields-cofields.R
\name{cmf_view_mol_field_cofield}
\alias{cmf_view_mol_field_cofield}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
View molecule, co-field and field
}
\description{
View molecule, co-field and field
}
\usage{
cmf_view_mol_field_cofield(
  mdb_fname, 
  imol,                                  
  ft,
  grid_fname, 
  alpha_from_model, 
  model_fname, 
  alpha, 
  rlevel,                             
  alpha_g,                             
  draw_field,                         
  draw_cofield,                       
  draw_overlap
)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mdb_fname}{"ligands-train.mol2"}
  \item{imol}{1}
  \item{ft}{"q"} 
  \item{grid_fname}{"ligands-grid-krr.RData"} 
  \item{alpha_from_model}{FALSE} 
  \item{model_fname}{"ligands-model.RData"} 
  \item{alpha}{0.3} 
  \item{rlevel}{0.3} 
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
if(FALSE)
{
# View overlapping fields of regression coefficients and molecular fields

#source("cmf-view-mol-field-cofield.R")

# Molecular database file name
mdb_fname <- "ligands-train.mol2"

# Molecule to visualize
imol <- 1

# Field type 
ft <- "vdw"

# Take alpha from model? (TRUE/FALSE)
alpha_from_model <- TRUE

# Model file name
model_fname <- "ligands-model.RData"


# Alpha value (if not taken from model)
alpha <- 1.0

# Isosurface level
rlevel <- 0.5

# Alpha (non-transperancy) level
alpha_g <- 0.5

cmf_view_mol_field_cofield(
  mdb_fname = mdb_fname, 
  imol = imol, 
  ft = ft, 
  grid_fname = "ligands-grid-krr.RData",
  alpha_from_model = alpha_from_model, 
  model_fname = model_fname, 
  alpha = alpha, 
  rlevel = rlevel, 
  alpha_g = alpha_g,
  draw_field = TRUE,
  draw_cofield = FALSE,
  draw_overlap = FALSE
)

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
