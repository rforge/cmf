% File cinf-mol2
\name{write_mol2}
\alias{write_mol2}
\title{Write Sybyl mol2-file}
\description{
  Writes file in Sybyl mol2 format.
}
\usage{
write_mol2(mdb, fname)
}
\arguments{
	\item{mdb}{molecular database (mdb)}
	\item{fname}{the name of the mol2 file to be written to disk}
}
\examples{
## Save mol2 file
data("alkans74")
\dontrun{
write_mol2(alkans74.mdb, "fname.mol2")
}
}
