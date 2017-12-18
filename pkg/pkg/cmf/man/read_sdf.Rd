% File cinf-sdf.R
\name{read_sdf}
\alias{read_sdf}
\title{Read sdf-file}
\description{
  Reads file in sdf format.
}
\usage{
read_sdf(filename, coord, prop, to_numeric, delete_expl_hydr)
}
\arguments{
	\item{filename}{the name of the sdf file to be read}
	\item{coord}{whether to read coordinates of atoms; default:TRUE}
	\item{prop}{whether to read property values; default:TRUE}
	\item{to_numeric}{whether to convert property values to numbers; default:FALSE}
	\item{delete_expl_hydr}{whether to delete explicit hydrogens; default:FALSE}
}
\value{Molecular database}
\examples{
## Load sdf file
\dontrun{
mdb <- read_sdf("fname.sdf")
}
}
