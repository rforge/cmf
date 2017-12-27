% File cinf-mol
\name{substr_mask}
\alias{substr_mask}
\title{Extract substructure from a molecule using mask}
\description{
  Extracts substructure from a molecule using mask.
}
\usage{
substr_mask(mol, mask)
}
\arguments{
	\item{mol}{molecule}
	\item{mask}{mask defining substructure to be extracted}
}
\value{Substructure}
\examples{
data("alkans74")
mol <- alkans74.mdb[[2]]
substr_mask(mol, c(1,1,0))
}
