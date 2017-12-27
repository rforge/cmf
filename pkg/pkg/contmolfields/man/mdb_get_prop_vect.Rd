% File cinf-mdb.R
\name{mdb_get_prop_vect}
\alias{mdb_get_prop_vect}
\title{Get property values}
\description{
  Returns vector of property values for specified molecular database and property name
}
\usage{
mdb_get_prop_vect(mdb, propname)
}
\arguments{
	\item{mdb}{Molecular database}
	\item{propname}{Property name}
}
\value{Vector of property values}
\examples{
data("alkans74")
mdb_get_prop_vect(alkans74.mdb, "bp")
}
