% File cmf-allign.R 
\name{align_mdb_template}
\alias{align_mdb_template}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Aligns molecular database mdb using template templ
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Aligns molecular database mdb using template templ
}
\usage{
align_mdb_template(mdb, templ, iimol = 1:length(mdb))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mdb}{
%%     ~~Describe \code{mdb} here~~
}
  \item{templ}{
%%     ~~Describe \code{templ} here~~
}
  \item{iimol}{
%%     ~~Describe \code{iimol} here~~
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
function (mdb, templ, iimol = 1:length(mdb)) 
{
    templ_ct <- mol_get_ct(templ)
    templ_lab <- mol_get_chelabs(templ)
    xyz_t <- mol2xyz(templ)
    mdb_a <- list()
    imol1 <- 0
    for (imol in iimol) {
        imol1 <- imol1 + 1
        mol <- mdb[[imol]]
        mol_ct <- mol_get_ct(mol)
        mol_lab <- mol_get_chelabs(mol)
        xyz_m <- mol2xyz(mol)
        isom_list <- find_substr_isomorph(templ_lab, templ_ct, 
            mol_lab, mol_ct)
        nisom <- length(isom_list)
        substr_mol <- substruct(mol, isom_list[[1]])
        xyz_ss <- mol2xyz(substr_mol)
        align <- align_arun(xyz_ss, xyz_t)
        xyz_ss_a <- transform_xyz(xyz_ss, align$R, align$T)
        rmse_a <- rmse4xyz(xyz_ss_a, xyz_t)
        align_best <- align
        rmse_a_best <- rmse_a
        if (nisom > 1) {
            for (isom in 2:nisom) {
                substr_mol <- substruct(mol, isom_list[[isom]])
                xyz_ss <- mol2xyz(substr_mol)
                align <- align_arun(xyz_ss, xyz_t)
                xyz_ss_a <- transform_xyz(xyz_ss, align$R, align$T)
                rmse_a <- rmse4xyz(xyz_ss_a, xyz_t)
                if (rmse_a < rmse_a_best) {
                  align_best <- align
                  rmse_a_best <- rmse_a
                }
            }
        }
        xyz_a <- transform_xyz(xyz_m, align_best$R, align_best$T)
        mol_a <- xyz2mol(mol, xyz_a)
        mdb_a[[imol1]] <- mol_a
        cat(sprintf("imol=\%d imol1=\%d nisom=\%d detx=\%g rmse1=\%g rmse2=\%g\n", 
            imol, imol1, nisom, align$detx, rmse4mol(substr_mol, 
                templ), rmse_a_best))
        flush.console()
    }
    mdb_a
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
