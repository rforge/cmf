% File cmf-allign.R
\name{align_arun}
\alias{align_arun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Rigin alignment with Arun algorithm
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Rigin alignment with Arun algorithm
}
\usage{
align_arun(p_i, p1_i)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p_i}{
%%     ~~Describe \code{p_i} here~~
}
  \item{p1_i}{
%%     ~~Describe \code{p1_i} here~~
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
function (p_i, p1_i) 
{
    p <- rowMeans(p_i)
    p1 <- rowMeans(p1_i)
    q_i <- p_i - p
    q1_i <- p1_i - p1
    H <- q_i \%*\% t(q1_i)
    s <- svd(H)
    U <- s$u
    V <- s$v
    X <- V \%*\% t(U)
    detx <- det(X)
    R <- X
    T <- p1 - R \%*\% p
    list(R = R, T = T, detx = detx)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
