% File cinf-molattribs.R
\name{add_mol_attribs}
\alias{add_mol_attribs}
\title{Add additional attributes to all atoms in all molecules in molecular database}
\description{
  Additional attributes are computed for all atoms for all molecules in molecular database.
  %\format{
    \tabular{ll}{
	  vd_ \tab vertex degree                       \cr
	  va_ \tab valence                             \cr
	  pi_ \tab the number of pi-electrons          \cr
	  ne_ \tab vector of neighbours                \cr
	  bo_ \tab vector of bond orders to neighbours 
	}
  %}
}
\usage{
add_mol_attribs(moldbase)
}
\arguments{
	\item{moldbase}{Molecular database}
}
\value{Molecular database with additional attributes appended to each atom}
\examples{
data("alkans74")
alkans74.mdb.1 <- add_mol_attribs(alkans74.mdb)
}
