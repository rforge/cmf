% File cinf-getprops.R
\name{get_prop_vec}
\alias{get_prop_vec}
\title{Get property values}
\description{
  Returns vector of property values for specified molecular database and property name
}
\usage{
get_prop_vec(moldbase, propname)
}
\arguments{
	\item{moldbase}{Molecular database}
	\item{propname}{Property name}
}
\value{Vector of property values}
\examples{
data("alkans74")
get_prop_vec(alkans74.mdb, "bp")
}
