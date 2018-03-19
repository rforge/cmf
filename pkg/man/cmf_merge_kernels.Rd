% File cmf-merge_kernels.R

\name{cmf_merge_kernels}
\alias{cmf_merge_kernels}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Merging files with computed kernels
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Merging files with computed kernels
}
\usage{
cmf_merge_kernels(kernels_train_fname_list, kernels_train_merged_fname)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{kernels_train_fname_list}{ }
  \item{kernels_train_merged_fname}{ }
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
## The function is currently defined as
\dontrun{
function(kernels_train_fname_list, kernels_train_merged_fname)
{
  kernels_merged <- list()
  for (kernels_fname in kernels_train_fname_list) {
    load(kernels_fname)
    field_names <- names(kernels)
    for (el in field_names) kernels_merged[[el]] <- kernels[[el]]
  }
  kernels <- kernels_merged
  save(kernels, file=kernels_train_merged_fname)
}

}
}
