% File cinf-mol
\name{mol_get_ct}
\alias{mol_get_ct}
\title{Return connection table for a molecule}
\description{
  Returns connection table for a molecule.
}
\usage{
mol_get_ct(mol, bond_orders)
}
\arguments{
	\item{mol}{molecule}
	\item{bond_orders}{whether to specify bond orders in connection table (1/0, default=0)}
}
\value{Connection table}
\examples{
data("alkans74")
mol <- alkans74.mdb[[1]]
mol_get_ct(mol)
}
