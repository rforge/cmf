% File cinf-free-wilson.R
\name{process_free_wilson}
\alias{process_free_wilson}
\title{Process data table using the Free-Wilson approach}
\description{
  Processes data table using the Free-Wilson approach.
}
\usage{
process_free_wilson(t, activity_column)
}
\arguments{
	\item{t}{data frame with table containing data for Free-Wilson analysis}
	\item{activity_column}{whether the data frame contains the activity data in its last column. Default: TRUE}
}
\value{Data frame with descriptors and properties}
\examples{
\dontrun{
  t <- read.table("free-wilson-boehm.txt", header=TRUE)
  df <- process_free_wilson(t)
}
}
