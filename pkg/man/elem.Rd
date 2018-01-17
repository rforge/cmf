% File cinf-descr-elem.R
\name{elem}
\alias{elem}
\title{Compute descriptors based on elemetal composition}
\description{
  Computes descriptors based on elemental composition
  %\format{
    \tabular{ll}{
	  nC  \tab the number of carbon atoms     \cr
	  nH  \tab the number of hydrogen atoms   \cr
	  nN  \tab the number of nitrogen atoms   \cr
	  nO  \tab the number of oxygen atoms     \cr
	  nS  \tab the number of sulfur atoms     \cr
	  nP  \tab the number of phosphorus atoms \cr
	  nSe \tab the number of selenium atoms   \cr
	  nF  \tab the number of fluorine atoms   \cr
	  nCl \tab the number of chlorine atoms   \cr
	  nBr \tab the number of bromine atoms    \cr
	  nI  \tab the number of iodine atoms
	}
  %}
}
\usage{
elem(moldbase)
}
\arguments{
	\item{moldbase}{molecular database}
}
\value{data frame with computed descriptors}
\examples{
data("alkans74")
df <- elem(alkans74.mdb)
}
