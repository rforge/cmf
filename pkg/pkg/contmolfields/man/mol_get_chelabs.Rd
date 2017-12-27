% File cinf-mol
\name{mol_get_chelabs}
\alias{mol_get_chelabs}
\title{Return array of chemical element labels}
\description{
  Returns array of chemical element labels
}
\usage{
mol_get_chelabs(mol)
}
\arguments{
	\item{mol}{molecule}
}
\value{Array of chemical element labels}
\examples{
data("alkans74")
mol <- alkans74.mdb[[1]]
mol_get_chelabs(mol)
}
