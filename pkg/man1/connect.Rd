% File cinf-descr-connect.R
\name{connect}
\alias{connect}
\title{Compute Kier-Hall connectivity topological indices}
\description{
  Computes Kier-Hall connectivity topological indices:
  \format{
    \tabular{ll}{
	  Chi0V   \tab Chi0V   \cr
	  Chi1V   \tab Chi1V   \cr
	  Chi2V   \tab Chi2V   \cr
	  Chi3cV  \tab Chi3cV  \cr
	  Chi3pV  \tab Chi3pV  \cr
	  Chi4pV  \tab Chi4pV  \cr
	  Chi4pcV \tab Chi4pcV \cr
	  Chi5cV  \tab Chi5cV  \cr
	  Chi5pV  \tab Chi5pV  \cr
	  Chi6pV  \tab Chi6pV
	}
  }
}
\usage{
connect(moldbase)
}
\arguments{
	\item{moldbase}{molecular database}
}
\value{data frame with computed topological indices}
\examples{
data("alkans74")
df <- connect(alkans74.mdb)
}
