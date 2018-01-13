# Reads TRIPOS force field

tripos_Rvdw <- list()
tripos_Evdw <- list()
tripos_atom_types <- list()

read_tripos_ff <- function(filename="data/cmf-tripos.prm") {
  lines <- readLines(filename)
  nlines <- length(lines)
  for (i in 1:nlines) {
    if (any(grep("atomtype", lines[i]))) {
      syb <- substr(lines[i], 11, 15)
      syb <- sub('[[:space:]]+$', '', syb)
      Rvdw_str <- substr(lines[i], 21, 27)
      Rvdw <- as.double(Rvdw_str)
      Evdw_str <- substring(lines[i], 29)
      Evdw <- as.double(Evdw_str)
      tripos_Rvdw[[syb]] <<- Rvdw
      tripos_Evdw[[syb]] <<- Evdw
	  tripos_atom_types <<- c(tripos_atom_types, syb)
    }
  }
}

read_tripos_ff()

# Extracts list of sybyl types of atoms
get_syb_types_list <- function(mdb)
{
  syb_types <- NULL
  nmol <- length(mdb)
  for (imol in 1:nmol) {
    mol <- mdb[[imol]]
    natoms <- length(mol$atoms)
	for (iatom in 1:natoms) {
	  atom <- mol$atoms[[iatom]]
	  if (!(atom$syb %in% syb_types)) {
	    syb_types <- c(syb_types, atom$syb)
	  }
	}
  }
  syb_types
}

# Extended version of get_syb_types_list
get_syb_types_list_ex <- function(mdb, min_occ=0, good_list=NULL, bad_list=NULL) {
  syb_types <- get_syb_types_list(mdb)
  
  # Check occurrencies
  if (min_occ) {
    # Find occurrencies
	occ_mdb <- list()
	for (st in syb_types) occ_mdb[[st]] <- 0
	occ_mol <- list()
	for (mol in mdb) {
	  for (st in syb_types) occ_mol[[st]] <- 0
	  for (atom in mol$atoms) occ_mol[[atom$syb]] <- occ_mol[[atom$syb]] + 1
	  for (st in syb_types) if (occ_mol[[st]] > 0) occ_mdb[[st]] <- occ_mdb[[st]] + 1
	}
	syb_types_1 <- NULL
	for (st in syb_types) if (occ_mdb[[st]] >= min_occ) syb_types_1 <- c(syb_types_1, st)
	syb_types <- syb_types_1
  }
  
  # Check good list
  if (length(good_list)) syb_types <- intersect(syb_types, good_list)
  
  # Check bad list
  if (length(bad_list)) syb_types <- setdiff(syb_types, bad_list)
  
  syb_types
}

