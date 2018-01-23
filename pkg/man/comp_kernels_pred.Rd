% File cmf-comp-kernels.R
\name{comp_kernels_pred}
\alias{comp_kernels_pred}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Computes CMF kernel matrices for prediction and saves them to file
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Computes CMF kernel matrices for prediction and saves them to file
}
\usage{
comp_kernels_pred(train_fname = "ligands-train.mol2", pred_fname = "ligands-pred.mol2", kernels_pred_fname = "ligands-kernels-pred.RData", mfields = c("q", "vdw", "logp", "abra", "abrb"), print_comp_kernels = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{train_fname}{
%%     ~~Describe \code{train_fname} here~~
}
  \item{pred_fname}{
%%     ~~Describe \code{pred_fname} here~~
}
  \item{kernels_pred_fname}{
%%     ~~Describe \code{kernels_pred_fname} here~~
}
  \item{mfields}{
%%     ~~Describe \code{mfields} here~~
}
  \item{print_comp_kernels}{
%%     ~~Describe \code{print_comp_kernels} here~~
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
function (train_fname = "ligands-train.mol2", pred_fname = "ligands-pred.mol2", 
    kernels_pred_fname = "ligands-kernels-pred.RData", mfields = c("q", 
        "vdw", "logp", "abra", "abrb"), print_comp_kernels = TRUE, 
    ...) 
{
    mdb0_train <- read_mol2(train_fname)
    mdb0_pred <- read_mol2(pred_fname)
    mdb_train <- cmf_params_tripos(mdb0_train)
    mdb_pred <- cmf_params_tripos(mdb0_pred)
    nfields <- length(mfields)
    syb_types <- get_syb_types_list(mdb_train)
    kernels_pred <- list()
    kernels_pred$alphas <- alphas
    for (f in 1:nfields) {
        kernels_pred[[mfields[f]]] <- list()
    }
    for (ialpha in 1:length(alphas)) {
        alpha <- alphas[ialpha]
        for (f in 1:nfields) {
            field <- mfields[f]
            if (print_comp_kernels) {
                cat(sprintf("computing kernel_\%s for alpha=\%g\n", 
                  field, alpha))
                flush.console()
            }
            if (field == "ind") {
                Km <- 0
                for (type in syb_types) {
                  if (print_comp_kernels) 
                    cat(type)
                  Km <- Km + cmf_indicator_kernel_matrix_pred(mdb_pred, 
                    mdb_train, alpha, type, verbose = print_comp_kernels)
                }
            }
            else {
                Km <- cmf_kernel_matrix_tp(field, mdb_pred, mdb_train, 
                  alpha, verbose = print_comp_kernels)
            }
            kernels_pred[[field]][[ialpha]] <- Km
        }
    }
    save(kernels_pred, file = kernels_pred_fname)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
