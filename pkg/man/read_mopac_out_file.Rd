% Fiel cmf-mopac.R
%read_mopac_out_file <- function(fname, mol) 
\name{read_mopac_out_file}
\alias{read_mopac_out_file}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read mopac output file
%%  ~~function to do ... ~~
}
\description{
Read mopac output file
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
read_mopac_out_file(fname, mol) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fname}{ }
  \item{x}{mol}
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
read_mopac_out_file <- function(fname, mol) {
 
  # proceed to next line
  next_line <- function() {
    iline <<- iline + 1
	aline <<- arc_lines[iline]
  }
  
  get_num_heavy_atoms <- function(mol) {
	iheavy <- 0	
	for (atom in mol$atoms) {
	  if (atom$el != "H") iheavy <- iheavy + 1
	}
	return(iheavy)
  }
  
  arc_lines <- readLines(fname)
  nlines <- length(arc_lines)
  iline <- 0
  
  # Skip to superdelocalizabilities
  while (iline < nlines) {
    next_line()
    r <- regexpr("           SUPERDELOCALIZABILITIES", aline)
	if (r > 0) break
  }
  
  # Read several scalar values
  next_line()
  next_line()
  mulliken_electronegativity <- as.numeric(substr(aline, 34, 44))
  next_line()
  parr_pople_absolute_hardness <- as.numeric(substr(aline, 34, 44))
  next_line()
  schuurmann_mo_shift_alpha <- as.numeric(substr(aline, 34, 44))
  next_line()
  next_line()
  ehomo <- as.numeric(substr(aline, 34, 44))
  next_line()
  elumo <- as.numeric(substr(aline, 34, 44))
  for (i in 1:4) next_line()
  
  # Read arrays with superdelocalizabilities
  num_heavy_atoms <- get_num_heavy_atoms(mol)
  heavy_atom_index <- integer(num_heavy_atoms)
  Dn <- double(num_heavy_atoms)
  De <- double(num_heavy_atoms)
  qZ <- double(num_heavy_atoms)
  piS <- double(num_heavy_atoms)
  chomo <- double(num_heavy_atoms)
  clumo <- double(num_heavy_atoms)
  # Read Dn(r), De(r) and q(r)-Z(r)
  for (i in 1:num_heavy_atoms) {
    next_line()
    heavy_atom_index[i] <- as.integer(substr(aline, 5, 7))
	Dn[i] <- as.numeric(substr(aline, 11, 20))
	De[i] <- as.numeric(substr(aline, 24, 33))
	qZ[i] <- as.numeric(substr(aline, 37, 46))
  }
  for (i in 1:5) next_line()
  # Read piS
  for (i in 1:num_heavy_atoms) {
    next_line()
	piS[i] <- as.numeric(substr(aline, 11, 20))
  }
  for (i in 1:5) next_line()
  # Read homo and lumo
  for (i in 1:num_heavy_atoms) {
    next_line()
	chomo[i] <- as.numeric(substr(aline, 24, 33))
	clumo[i] <- as.numeric(substr(aline, 37, 46))
  }
  
  list(
    num_heavy_atoms = num_heavy_atoms,
	heavy_atom_index = heavy_atom_index,
	Dn = Dn,
	De = De,
	qZ = qZ,
	piS = piS,
	chomo = chomo,
	clumo = clumo,
	mulliken_electronegativity = mulliken_electronegativity,
	parr_pople_absolute_hardness = parr_pople_absolute_hardness,
	schuurmann_mo_shift_alpha = schuurmann_mo_shift_alpha,
	ehomo = ehomo,
	elumo = elumo
  )
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
