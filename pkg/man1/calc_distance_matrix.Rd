% File cinf-mol.R
\name{calc_distance_matrix}
\alias{calc_distance_matrix}
\title{Compute distance matrix from connection table}
\description{
  Computes distance matrix from connection table.
}
\usage{
calc_distance_matrix(connTable)
}
\arguments{
	\item{connTable}{connection table}
}
\value{distance matrix}
\examples{
data("alkans74")
mol <- alkans74.mdb[[5]]
ct <- mol_get_ct(mol)
calc_distance_matrix(ct)
}
