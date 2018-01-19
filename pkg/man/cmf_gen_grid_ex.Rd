% File cmf-gen-grid.R
\name{cmf_gen_grid_ex}
\alias{cmf_gen_grid_ex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Generation of grid for continuous molecular co-fields
Extended version with field families
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Generation of grid for continuous molecular co-fields
Extended version with field families
}
\usage{
cmf_gen_grid_ex(train_fname = "ligands-train.mol2", kernels_fname = "ligands-kernels.RData", model_fname = "ligands-model.RData", grid_fname = "ligands-grid-krr.RData", field_family = "PHCH", verbose = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{train_fname}{
%%     ~~Describe \code{train_fname} here~~
}
  \item{kernels_fname}{
%%     ~~Describe \code{kernels_fname} here~~
}
  \item{model_fname}{
%%     ~~Describe \code{model_fname} here~~
}
  \item{grid_fname}{
%%     ~~Describe \code{grid_fname} here~~
}
  \item{field_family}{
%%     ~~Describe \code{field_family} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
function (train_fname = "ligands-train.mol2", kernels_fname = "ligands-kernels.RData", 
    model_fname = "ligands-model.RData", grid_fname = "ligands-grid-krr.RData", 
    field_family = "PHCH", verbose = TRUE, ...) 
{
    load(kernels_fname)
    load(model_fname)
    mfields <- names(model$h)
    nfields <- length(mfields)
    mdb <- read_mol2(train_fname)
    if (field_family == "PHCH") 
        mdb <- cmf_params_tripos(mdb)
    grid <- cmf_init_grid(mdb)
    grids <- list()
    for (f in 1:nfields) {
        field <- mfields[f]
        if (verbose) {
            cat(sprintf("Generating grid for field \%s...\n", 
                field))
            flush.console()
        }
        grids[[field]] <- cmf_coef_grid_ex(mdb, model$a, model$alpha[[field]], 
            grid, field, field_family)
    }
    save(grids, file = grid_fname)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
