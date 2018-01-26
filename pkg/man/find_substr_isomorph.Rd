% File cinf-isomorph.R
\name{find_substr_isomorph}
\alias{find_substr_isomorph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Finds all substructure isomorphisms
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Finds all substructure isomorphisms
}
\usage{
find_substr_isomorph(substr_lab, substr_ct, str_lab, str_ct)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{substr_lab}{
%%     ~~Describe \code{substr_lab} here~~
}
  \item{substr_ct}{
%%     ~~Describe \code{substr_ct} here~~
}
  \item{str_lab}{
%%     ~~Describe \code{str_lab} here~~
}
  \item{str_ct}{
%%     ~~Describe \code{str_ct} here~~
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
function (substr_lab, substr_ct, str_lab, str_ct) 
{
    isomorph <- function(stp) {
        for (i in 1:str_size) if (!used[i] && substr_lab[stp] == 
            str_lab[i]) {
            cc[stp] <<- i
            if (stp > 1) 
                for (j in 1:(stp - 1)) if ((substr_ct[stp, j] != 
                  0) && (substr_ct[stp, j] != str_ct[i, cc[j]])) {
                  to_exit <- TRUE
                  next
                }
            if (to_exit) {
                to_exit <- FALSE
                next
            }
            if (stp == substr_size) {
                num_matches <<- num_matches + 1
                isom_list[[num_matches]] <<- cc
            }
            else {
                used[i] <<- TRUE
                isomorph(stp + 1)
                used[i] <<- FALSE
            }
        }
    }
    substr_size <- length(substr_lab)
    str_size <- length(str_lab)
    used <- logical(str_size)
    cc <- integer(substr_size)
    to_exit <- FALSE
    num_matches <- 0
    isom_list <- list()
    isomorph(1)
    isom_list
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
