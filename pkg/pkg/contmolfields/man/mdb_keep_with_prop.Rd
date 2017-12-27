% File cinf-mdb.R
\name{mdb_keep_with_prop}
\alias{mdb_keep_with_prop}
\title{Keep in molecular database only compounds containing values of a certain property}
\description{
  Returns molecular databased composed of compounds possessing values of a given property
}
\usage{
mdb_keep_with_prop(mdb, prop_to_keep)
}
\arguments{
	\item{mdb}{Molecular database}
	\item{prop_to_keep}{Property that should be specified in all compounds of the returned molecular database}
}
\value{Molecular database}
\examples{
data("alkans74")
mdb_keep_with_prop(alkans74.mdb, "mp")
}
