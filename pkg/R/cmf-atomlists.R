# Generates atom lists for different types of molecular fields

# Generate atom lists for generic molecular field
gen_atomlists_generic <- function(mdb) {
  nmols <- length(mdb)
  atomlists <- list()
  for (imol in 1:nmols) {
    mol <- mdb[[imol]]
    natoms <- length(mol$atoms)
	atomlists[[imol]] <- 1:natoms
  }
  atomlists
}

# Generate atom lists for physico-chemical molecular field
gen_atomlists_phch <- function(mdb, field) {
  gen_atomlists_generic(mdb)
}

# Generate atom lists for continuous indicator field
gen_atomlists_ind <- function(mdb, field) {
  nmols <- length(mdb)
  atomlists <- list()
  for (imol in 1:nmols) {
    mol <- mdb[[imol]]
    natoms <- length(mol$atoms)
	al <- list()
	for (iatom in 1:natoms) {
	  atom <- mol$atoms[[iatom]]
	  if (atom$syb == field) al <- c(al, iatom)
	}
	atomlists[[imol]] <- al
  }
  atomlists
}

# Generate atom lists for MOPAC molecular field
gen_atomlists_mopac <- function(mdb, field) {
  gen_atomlists_generic(mdb)
}

# Generate atom lists
gen_atomlists <- function(mdb, field, field_family) {
  atomlists <- list()
  if (field_family == "PHCH") {
    atomlists <- gen_atomlists_phch(mdb, field)
  } else if (field_family == "IND") {
    atomlists <- gen_atomlists_ind(mdb, field)
  } else if (field_family == "MOPAC") {
    atomlists <- gen_atomlists_mopac(mdb, field)
  }
  atomlists
}
