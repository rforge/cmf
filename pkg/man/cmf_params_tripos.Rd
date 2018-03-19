%% File cmf-params.R
\name{cmf_params_tripos}
\alias{cmf_params_tripos}

\title{
Parameters for CMF
}
\description{
Parameters for CMF
}
\usage{
cmf_params_tripos(mdb)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mdb}{ }
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
# Parameters for CMF

cmf_params_tripos <- function(mdb) {
  ncomp <- length(mdb)
  for (imol in 1:ncomp) {
    mol <- mdb[[imol]]

    # count neighbouring hydrogens
    natoms <- length(mol$atoms)
    nchydr <- integer(natoms)
    nbonds <- length(mol$bonds)
    for (ibond in 1:nbonds) {
      bond <- mol$bonds[[ibond]]
      iat1 <- bond$at1
      iat2 <- bond$at2
      el1 <- mol$atoms[[iat1]]$el
      el2 <- mol$atoms[[iat1]]$el
      if (el1 == "H") {
        nchydr[iat2] <- nchydr[iat2] + 1
      }
      if (el2 == "H") {
        nchydr[iat1] <- nchydr[iat1] + 1
      }
    }

    # assign parameters to atoms
    for (iatom in 1:natoms) {
      atom <- mol$atoms[[iatom]]

      atom$hydroph <- 0.0
      atom$abraham_a <- 0.0
      atom$abraham_b <- 0.0
      atom$abraham_e <- 0.0
      atom$abraham_s <- 0.0

      if (atom$syb == "C.1") {
        atom$hydroph <- atom$hydroph +0.041270
        atom$abraham_e <- atom$abraham_e -0.005210
        atom$abraham_s <- atom$abraham_s -0.013906
        if (nchydr[iatom] == 1) {
          atom$hydroph <- atom$hydroph -0.093156
          atom$abraham_s <- atom$abraham_s +0.045181
        }
      }

      if (atom$syb == "C.2") {
        atom$hydroph <- atom$hydroph +0.041270
        atom$abraham_e <- atom$abraham_e -0.005210
        atom$abraham_s <- atom$abraham_s -0.013906
        atom$hydroph <- atom$hydroph +0.006410
        atom$abraham_s <- atom$abraham_s +0.016508
        if (nchydr[iatom] == 0) {
          atom$hydroph <- atom$hydroph +0.007629
          atom$abraham_b <- atom$abraham_b +0.026654
          atom$abraham_e <- atom$abraham_e +0.076572
        }
        if (nchydr[iatom] == 1) {
          atom$hydroph <- atom$hydroph +0.202552
          atom$abraham_e <- atom$abraham_e +0.058299
        }
        if (nchydr[iatom] == 2) {
          atom$hydroph <- atom$hydroph +0.422564
          atom$abraham_b <- atom$abraham_b -0.041897
          atom$abraham_e <- atom$abraham_e -0.019698
          atom$abraham_s <- atom$abraham_s -0.098377
        }
      }

      if (atom$syb == "C.3") {
        atom$hydroph <- atom$hydroph +0.041270
        atom$abraham_e <- atom$abraham_e -0.005210
        atom$abraham_s <- atom$abraham_s -0.013906
        atom$abraham_a <- atom$abraham_a -0.009684
        if (nchydr[iatom] == 0) {
          atom$hydroph <- atom$hydroph +0.089007
          atom$abraham_b <- atom$abraham_b +0.096834
          atom$abraham_e <- atom$abraham_e +0.153114
        }
        if (nchydr[iatom] == 1) {
          atom$abraham_b <- atom$abraham_b +0.041221
          atom$abraham_e <- atom$abraham_e +0.091218
        }
        if (nchydr[iatom] == 2) {
          atom$hydroph <- atom$hydroph +0.208288
        }
        if (nchydr[iatom] == 3) {
          atom$hydroph <- atom$hydroph +0.369764
          atom$abraham_b <- atom$abraham_b -0.022123
          atom$abraham_e <- atom$abraham_e -0.102229
          atom$abraham_s <- atom$abraham_s -0.055649
        }
      }

      if (atom$syb == "C.ar") {
        atom$hydroph <- atom$hydroph +0.041270
        atom$abraham_e <- atom$abraham_e -0.005210
        atom$abraham_s <- atom$abraham_s -0.013906
        atom$hydroph <- atom$hydroph +0.194585
        atom$abraham_e <- atom$abraham_e +0.117666
        atom$abraham_s <- atom$abraham_s +0.070704
        if (nchydr[iatom] == 0) {
          atom$hydroph <- atom$hydroph +0.085378
          atom$abraham_b <- atom$abraham_b +0.002903
          atom$abraham_e <- atom$abraham_e +0.107527
          atom$abraham_s <- atom$abraham_s +0.046501
        }
        if (nchydr[iatom] == 1) {
          atom$abraham_e <- atom$abraham_e -0.044799
          atom$abraham_s <- atom$abraham_s -0.011864
        }
      }

      if (atom$syb == "C.cat") {
        atom$abraham_e <- atom$abraham_e -0.005210
        atom$abraham_s <- atom$abraham_s -0.013906
      }

      if (atom$syb == "N.1") {
        atom$abraham_a <- atom$abraham_a +0.077667
        atom$abraham_s <- atom$abraham_s +0.167948
        atom$abraham_b <- atom$abraham_b +0.061824
        atom$abraham_s <- atom$abraham_s +0.232745
      }

      if (atom$syb == "N.2") {
        atom$abraham_a <- atom$abraham_a +0.077667
        atom$abraham_s <- atom$abraham_s +0.167948
        if (nchydr[iatom] == 1) {
          atom$hydroph <- atom$hydroph -0.200643
        } else {
          atom$abraham_b <- atom$abraham_b -0.763354
          atom$abraham_e <- atom$abraham_e +0.682859
        }
      }

      if (atom$syb == "N.3") {
        atom$abraham_a <- atom$abraham_a +0.077667
        atom$abraham_s <- atom$abraham_s +0.167948
        atom$hydroph <- atom$hydroph -0.091021
        atom$abraham_b <- atom$abraham_b +0.367558
        atom$abraham_e <- atom$abraham_e +0.099602
        atom$abraham_s <- atom$abraham_s +0.084806
        if (nchydr[iatom] == 0) {
          atom$hydroph <- atom$hydroph -0.495588
          atom$abraham_b <- atom$abraham_b +0.150946
          atom$abraham_e <- atom$abraham_e +0.088502
        }
        if (nchydr[iatom] == 1) {
          atom$hydroph <- atom$hydroph -0.452520
          atom$abraham_a <- atom$abraham_a +0.163329
        }
        if (nchydr[iatom] == 2) {
          atom$hydroph <- atom$hydroph -0.763440
          atom$abraham_a <- atom$abraham_a +0.120801
          atom$abraham_s <- atom$abraham_s -0.069707
        }
      }

      if (atom$syb == "N.4") {
        atom$abraham_a <- atom$abraham_a +0.077667
        atom$abraham_s <- atom$abraham_s +0.167948
        atom$hydroph <- atom$hydroph -4.382696
      }

      if (atom$syb == "N.am") {
        atom$abraham_a <- atom$abraham_a +0.077667
        atom$abraham_s <- atom$abraham_s +0.167948
        atom$hydroph <- atom$hydroph -0.091021
        if (nchydr[iatom] == 0) {
          atom$hydroph <- atom$hydroph -0.495588
          atom$abraham_a <- atom$abraham_a -0.096006
        }
        if (nchydr[iatom] == 1) {
          atom$hydroph <- atom$hydroph -0.452520
        }
        if (nchydr[iatom] == 2) {
          atom$hydroph <- atom$hydroph -0.763440
        }
      }

      if (atom$syb == "N.ar") {
        atom$abraham_a <- atom$abraham_a +0.077667
        atom$abraham_s <- atom$abraham_s +0.167948
        atom$abraham_a <- atom$abraham_a -0.104421
        atom$abraham_b <- atom$abraham_b +0.282666
      }

      if (atom$syb == "N.pl3") {
        atom$abraham_a <- atom$abraham_a +0.077667
        atom$abraham_s <- atom$abraham_s +0.167948
        atom$hydroph <- atom$hydroph +0.483597
        atom$abraham_b <- atom$abraham_b -0.296874
        atom$abraham_e <- atom$abraham_e +0.060327
        atom$abraham_s <- atom$abraham_s -0.351281
      }

      if (atom$syb == "O.2") {
        atom$hydroph <- atom$hydroph -0.224286
        atom$hydroph <- atom$hydroph -0.098496
        atom$abraham_a <- atom$abraham_a +0.038390
        atom$abraham_e <- atom$abraham_e -0.011749
        atom$abraham_s <- atom$abraham_s +0.371065
      }

      if (atom$syb == "O.3") {
        atom$hydroph <- atom$hydroph -0.224286
        atom$hydroph <- atom$hydroph -0.102250
        atom$abraham_s <- atom$abraham_s +0.111253
        atom$abraham_b <- atom$abraham_b +0.147962
        if (nchydr[iatom] == 0) {
          atom$hydroph <- atom$hydroph -0.010867
          atom$abraham_a <- atom$abraham_a -0.035354
          atom$abraham_e <- atom$abraham_e -0.020649
        } else {
          atom$hydroph <- atom$hydroph -0.279404
          atom$abraham_a <- atom$abraham_a +0.445333
          atom$abraham_b <- atom$abraham_b +0.014953
          atom$abraham_e <- atom$abraham_e +0.010690
          atom$abraham_s <- atom$abraham_s +0.027754
        }
      }

      if (atom$syb == "O.co2") {
        atom$hydroph <- atom$hydroph -0.224286
      }

      if (atom$syb == "S.2") {
        atom$hydroph <- atom$hydroph -0.077480
        atom$hydroph <- atom$hydroph +0.734821
        atom$abraham_e <- atom$abraham_e +0.136585
      }

      if (atom$syb == "S.3") {
        atom$abraham_e <- atom$abraham_e +0.195277
        atom$abraham_s <- atom$abraham_s +0.082278
        if (nchydr[iatom] == 0) {
          atom$hydroph <- atom$hydroph +0.504755
          atom$abraham_b <- atom$abraham_b +0.103649
          atom$abraham_e <- atom$abraham_e +0.071103
        }
      }

      if (atom$syb == "S.o") {
        atom$hydroph <- atom$hydroph -0.077480
        atom$hydroph <- atom$hydroph -1.055756
      }

      if (atom$syb == "S.o2") {
        atom$hydroph <- atom$hydroph -0.077480
        atom$hydroph <- atom$hydroph +0.027576
      }

      if (atom$syb == "P.3") {
        atom$abraham_b <- atom$abraham_b +0.686142
        atom$abraham_e <- atom$abraham_e +0.103799
        atom$abraham_s <- atom$abraham_s +0.215383
      }

      if (atom$syb == "H") {
      }

      if (atom$syb == "F") {
        atom$hydroph <- atom$hydroph +0.269012
        atom$abraham_b <- atom$abraham_b -0.077060
        atom$abraham_e <- atom$abraham_e -0.220887
        atom$abraham_s <- atom$abraham_s -0.077487
      }

      if (atom$syb == "Cl") {
        atom$hydroph <- atom$hydroph +0.569605
        atom$abraham_b <- atom$abraham_b -0.068040
        atom$abraham_e <- atom$abraham_e -0.008902
        atom$abraham_s <- atom$abraham_s +0.029923
      }

      if (atom$syb == "Br") {
        atom$hydroph <- atom$hydroph +0.696413
        atom$abraham_b <- atom$abraham_b -0.074564
        atom$abraham_e <- atom$abraham_e +0.167323
        atom$abraham_s <- atom$abraham_s +0.131222
      }

      if (atom$syb == "I") {
        atom$hydroph <- atom$hydroph +0.323219
        atom$abraham_e <- atom$abraham_e +0.509975
        atom$abraham_s <- atom$abraham_s +0.135703
      }

      if (atom$syb == "Si") {
        atom$hydroph <- atom$hydroph +0.617895
      }

      mdb[[imol]]$atoms[[iatom]] <- atom
    }

  }
  mdb
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
