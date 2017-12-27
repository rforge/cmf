% File cinf-getprops.R
\name{get_prop_names}
\alias{get_prop_names}
\title{Get the list of property names}
\description{
  Returns the list of property names from molecular database
}
\usage{
get_prop_names(moldbase)
}
\arguments{
	\item{moldbase}{Molecular database}
}
\value{Array of property names}
\examples{
data("alkans74")
get_prop_names(alkans74.mdb)
}
