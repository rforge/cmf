% File cinf-mol2
\name{read_mol2}
\alias{read_mol2}
\title{Read Sybyl mol2-file}
\description{
  Reads file in Sybyl mol2 format.
}
\usage{
read_mol2(filename)
}
\arguments{
	\item{filename}{the name of the mol2 file to be read}
}
\value{Molecular database (mdb)}
\examples{
## Load mol2 file
\dontrun{
mdb <- read_mol2("fname.mol2")
}
}
