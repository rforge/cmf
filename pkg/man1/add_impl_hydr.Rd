% File cinf-hydrogens.R
\name{add_impl_hydr}
\alias{add_impl_hydr}
\title{Add implicit hydrogens}
\description{
  Adds implicit hydrogens.
}
\usage{
add_impl_hydr(mol)
}
\arguments{
	\item{mol}{Molecule}
}
\value{Molecule with added implicit hydrogen atoms}
\examples{
data("alkans74")
mol <- alkans74.mdb[[10]]
mol_h <- add_impl_hydr(mol)
}
