% File cinf-sdf.R
\name{write_sdf}
\alias{write_sdf}
\title{Write sdf-file}
\description{
  Writes file in sdf format.
}
\usage{
write_sdf(mdb, filename)
}
\arguments{
	\item{mdb}{molecular database}
	\item{filename}{the name of the sdf file to be written to disk}
}
\examples{
## Save sdf file
data("alkans74")
\dontrun{
write_sdf(alkans74.mdb, "fname.sdf")
}
}
