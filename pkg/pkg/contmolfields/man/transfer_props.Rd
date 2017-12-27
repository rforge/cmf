% File cinf-mdb.R
\name{transfer_props}
\alias{transfer_props}
\title{Transfer properties from one mdb to another}
\description{
  Transfers properties from one mdb to another
}
\usage{
transfer_props(mdb_target, mdb_source)
}
\arguments{
	\item{mdb_target}{Target molecular database for property transfer}
	\item{mdb_source}{Source molecular database for property transfer}
}
\value{Target molecular database with transferred properties}
\examples{
\dontrun{
mdb_t <- transfer_props(mdb_t, mdb_s)
}
}
