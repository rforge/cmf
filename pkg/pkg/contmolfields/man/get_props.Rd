% File cinf-getprops.R
\name{get_props}
\alias{get_props}
\title{Get data frame with property values}
\description{
  Returns a data frame for the specified property names
}
\usage{
get_props(moldbase, propnames)
}
\arguments{
	\item{moldbase}{Molecular database}
	\item{propnames}{The list of property names. Default - the list of all property names}
}
\value{Data frame with property values}
\examples{
data("alkans74")
get_props(alkans74.mdb, c("bp", "mp"))
}
