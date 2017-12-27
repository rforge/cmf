% File cinf-mdb.R
\name{mdb_get_prop_matr1}
\alias{mdb_get_prop_matr1}
\title{Extract property one-column matrix by name}
\description{
  Extracts property one-column matrix by name
}
\usage{
mdb_get_prop_matr1(mdb, propname)
}
\arguments{
	\item{mdb}{Molecular database}
	\item{propname}{Property name}
}
\value{One-column matrix of property values}
\examples{
data("alkans74")
mdb_get_prop_matr1(alkans74.mdb, "bp")
}
