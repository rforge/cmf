# Various procedures with metals

# Delete metal from molecule by name
mol_delete_metal <- function(mol, mname) {
  natoms <- length(mol$atoms)
  mask <- rep(TRUE, natoms)
  for (iatom in 1:natoms) {
    atom <- mol$atoms[[iatom]]
	if (atom$el == mname) {mask[iatom] <- FALSE}
  }
  substr_mask(mol, mask)
}

# Delete metal from molecular database by name
mdb_delete_metal <- function(mdb, mname) {
  ncomp <- length(mdb)
  mdb_new <- list()
  for (imol in 1:ncomp) {
    mol <- mdb[[imol]]
	mol_new <- mol_delete_metal(mol, mname)
	mdb_new[[imol]] <- mol_new 
  }
  mdb_new
}


