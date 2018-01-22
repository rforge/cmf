\name{euler_orth}
\alias{euler_orth}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Generation of orthogonal rotation matrix for Euler angles phi, teta, psi
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Generation of orthogonal rotation matrix for Euler angles phi, teta, psi
}
\usage{
euler_orth(phi = 0, theta = 0, psi = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phi}{
%%     ~~Describe \code{phi} here~~
}
  \item{theta}{
%%     ~~Describe \code{theta} here~~
}
  \item{psi}{
%%     ~~Describe \code{psi} here~~
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
function (phi = 0, theta = 0, psi = 0) 
{
    r <- matrix(0, 3, 3)
    r[1, 1] <- cos(theta) * cos(psi)
    r[1, 2] <- -cos(phi) * sin(psi) + sin(phi) * sin(theta) * 
        cos(psi)
    r[1, 3] <- sin(phi) * sin(psi) + cos(phi) * sin(theta) * 
        cos(psi)
    r[2, 1] <- cos(theta) * sin(psi)
    r[2, 2] <- cos(phi) * cos(psi) + sin(phi) * sin(theta) * 
        sin(psi)
    r[2, 3] <- -sin(phi) * cos(psi) + cos(phi) * sin(theta) * 
        sin(psi)
    r[3, 1] <- -sin(theta)
    r[3, 2] <- sin(phi) * cos(theta)
    r[3, 3] <- cos(phi) * cos(theta)
    r
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
