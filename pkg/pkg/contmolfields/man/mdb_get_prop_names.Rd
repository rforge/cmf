% File cinf-mdb.R
\name{mdb_get_prop_names}
\alias{mdb_get_prop_names}
\title{Get the list of property names}
\description{
  Returns the list of property names from molecular database
}
\usage{
mdb_get_prop_names(mdb)
}
\arguments{
	\item{mdb}{Molecular database}
}
\value{Array of property names}
\examples{
data("alkans74")
mdb_get_prop_names(alkans74.mdb)
}
