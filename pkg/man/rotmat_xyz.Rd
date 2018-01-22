\name{rotmat_xyz}
\alias{rotmat_xyz}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Generation of rotation matrix around x,y,z axes:
alpha_(yaw), beta_(pitch), gamma_(roll)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Generation of rotation matrix around x,y,z axes:
alpha_(yaw), beta_(pitch), gamma_(roll)
}
\usage{
rotmat_xyz(alpha_ = 0, beta_ = 0, gamma_ = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alpha_}{
%%     ~~Describe \code{alpha_} here~~
}
  \item{beta_}{
%%     ~~Describe \code{beta_} here~~
}
  \item{gamma_}{
%%     ~~Describe \code{gamma_} here~~
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
function (alpha_ = 0, beta_ = 0, gamma_ = 0) 
{
    rx <- matrix(0, 3, 3)
    rx[1, 1] <- 1
    rx[1, 2] <- 0
    rx[1, 3] <- 0
    rx[2, 1] <- 0
    rx[2, 2] <- cos(gamma_)
    rx[2, 3] <- -sin(gamma_)
    rx[3, 1] <- 0
    rx[3, 2] <- sin(gamma_)
    rx[3, 3] <- cos(gamma_)
    ry <- matrix(0, 3, 3)
    ry[1, 1] <- cos(beta_)
    ry[1, 2] <- 0
    ry[1, 3] <- sin(beta_)
    ry[2, 1] <- 0
    ry[2, 2] <- 1
    ry[2, 3] <- 0
    ry[3, 1] <- -sin(beta_)
    ry[3, 2] <- 0
    ry[3, 3] <- cos(beta_)
    rz <- matrix(0, 3, 3)
    rz[1, 1] <- cos(alpha_)
    rz[1, 2] <- -sin(alpha_)
    rz[1, 3] <- 0
    rz[2, 1] <- sin(alpha_)
    rz[2, 2] <- cos(alpha_)
    rz[2, 3] <- 0
    rz[3, 1] <- 0
    rz[3, 2] <- 0
    rz[3, 3] <- 1
    r <- rx \%*\% ry \%*\% rz
    r
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
