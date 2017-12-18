% File cinf-mdb.R
\name{mdb_keep_matr_with_prop}
\alias{mdb_keep_matr_with_prop}
\title{Keep in matrix only rows and columns corresponding to compounds possessing values of certain property}
\description{
  Keep in matrix only rows and columns corresponding to compounds possessing values of certain property
}
\usage{
mdb_keep_matr_with_prop(matr, mdb, prop_to_keep)
}
\arguments{
	\item{matr}{Matrix to be changed}
	\item{mdb}{Molecular database}
	\item{prop_to_keep}{Property that should be specified in all compounds of the returned molecular database}
}
\value{Matrix after change}
\examples{
data("alkans74")
matr <- matrix(0, nrow=74, ncol=74)
mdb_keep_matr_with_prop(matr, alkans74.mdb, "mp")
}
