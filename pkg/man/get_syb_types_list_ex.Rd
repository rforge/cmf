% File cmf-triposff.R
\name{get_syb_types_list_ex}
\alias{get_syb_types_list_ex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extended version of get_syb_types_list
}
\description{
Extended version of get_syb_types_list
}
\usage{
get_syb_types_list_ex(mdb, min_occ=min_occ, good_list=good_list, bad_list=bad_list)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mdb}{ }
   \item{min_occ}{0}
   \item{good_list}{NULL}
   \item{bad_list}{NULL}
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
# Remove if(FALSE) block while testing, please notice file "ligands-train.mol2"
#should exist

# Training set file name
train_fname <- "data/ligands-train.mol2"

# Minimal occurrence of sybyl atom types
min_occ <- 7

# Good list of sybyl atom types
good_list <- NULL

# Bad list of sybyl atom types
bad_list <- NULL

# Types of molecular fields
mdb0 <- read_mol2(train_fname)
mdb <- cmf_params_tripos(mdb0)
mfields <- get_syb_types_list_ex(mdb, min_occ=min_occ, good_list=good_list, bad_list=bad_list)


## The function is currently defined as
get_syb_types_list_ex <- function(mdb, min_occ=0, good_list=NULL, bad_list=NULL) {
  syb_types <- get_syb_types_list(mdb)
  
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
