% File cinf-hydrogens.R
\name{del_expl_hydr}
\alias{del_expl_hydr}
\title{Delete explicit hydrogens}
\description{
  Deletes explicit hydrogens.
}
\usage{
del_expl_hydr(mol, add_impl, recalc_attribs)
}
\arguments{
	\item{mol}{Molecule}
	\item{add_impl}{Whether to add implicit hydrogens. Dafault: TRUE}
	\item{recalc_attribs}{Whether to recalculate additional attributes for atoms. Dafault: FALSE}
}
\value{Molecule without explicit hydrogen atoms}
\examples{
data("alkans74")
mol <- alkans74.mdb[[10]]
mol_noh <- del_expl_hydr(mol)
}
