% File cmf-allign.R
\name{superpose_mol}
\alias{superpose_mol}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Superposes moving molecule mol_m on templace molecule mol_t
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Superposes moving molecule mol_m on templace molecule mol_t
}
\usage{
superpose_mol(mol_m, mol_t)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mol_m}{
%%     ~~Describe \code{mol_m} here~~
}
  \item{mol_t}{
%%     ~~Describe \code{mol_t} here~~
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
function (mol_m, mol_t) 
{
    xyz_m <- mol2xyz(mol_m)
    xyz_t <- mol2xyz(mol_t)
    align <- align_arun(xyz_m, xyz_t)
    xyz_a <- transform_xyz(xyz_m, align$R, align$T)
    xyz2mol(mol_m, xyz_a)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
