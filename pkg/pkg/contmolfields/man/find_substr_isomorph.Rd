% File cinf-רûשüשךחנ.R
\name{find_substr_isomorph}
\alias{find_substr_isomorph}
\title{Find all substructure isomorphisms}
\description{
  Finds all substructure isomorphisms.
}
\usage{
find_substr_isomorph(substr_lab, substr_ct, str_lab, str_ct)
}
\arguments{
	\item{substr_lab}{labels for substructure}
	\item{substr_ct}{connection table for substructure}
	\item{str_lab}{labels for structure}
	\item{str_ct}{connection table for structure}
}
\value{list of isomorphisms}
\examples{
# Finds mol2 in mol1

# mol1
mol1_lab <- c("C", "C", "C", "C", "C", "N", "O")
mol1_ct <- matrix(0, nrow=7, ncol=7)
mol1_ct[1,2] <- mol1_ct[2,1] <- 1
mol1_ct[2,3] <- mol1_ct[3,2] <- 1
mol1_ct[3,4] <- mol1_ct[4,3] <- 1
mol1_ct[4,5] <- mol1_ct[5,4] <- 1
mol1_ct[5,6] <- mol1_ct[6,5] <- 2
mol1_ct[1,6] <- mol1_ct[6,1] <- 1
mol1_ct[4,7] <- mol1_ct[7,4] <- 2

# mol2
mol2_lab <- c("O", "C", "C", "N")
mol2_ct <- matrix(0, nrow=4, ncol=4)
mol2_ct[1,2] <- mol2_ct[2,1] <- 2
mol2_ct[2,3] <- mol2_ct[3,2] <- 1
mol2_ct[3,4] <- mol2_ct[4,3] <- 2

isomorph <- find_substr_isomorph(mol2_lab, mol2_ct, mol1_lab, mol1_ct)
}

