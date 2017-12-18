% File cinf-mdb.R
\name{mdb_get_num_comp}
\alias{mdb_get_num_comp}
\title{Get the number of compounds in molecular database}
\description{
  Returns the number of compounds in molecular database
}
\usage{
mdb_get_num_comp(mdb)
}
\arguments{
	\item{mdb}{Molecular database}
}
\value{The number of compounds}
\examples{
data("alkans74")
mdb_get_num_comp(alkans74.mdb)
}
